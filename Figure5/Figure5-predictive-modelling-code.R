####################
#################### Predictive modelling

### packages load
library(caret)
library(rms)
library(MASS)
library(pROC)

### data preparation
load("Figure5-rawdata-clinic-info.RData")

# random plit data by 7:3 ratio
set.seed(123)

trainIndex <- createDataPartition(clinic_all$Group, p=0.7, list=FALSE)
train_set <- clinic_all[trainIndex,]
test_set <- clinic_all[-trainIndex,]

########## logistic regression modelling

# feature selection
X <- clinic_all[, 1:14]  # features
y <- clinic_all[, 15]    # labels
control <- rfeControl(functions = caretFuncs, method = "cv", number = 10)

modelFuncs <- caretFuncs
set.seed(123)
rfeResults <- rfe(x=X, y=y, sizes=1:17, rfeControl=control, methods="glm")
print(rfeResults)

### logistic regression model construction
ddist <- datadist(train_set)
options(datadist='ddist')

lrm_model_6 <- lrm(Group ~differentiation.grade +CRC.site+Pathological.Type,data = train_set, x = TRUE, y = TRUE)

lrm_model_11 <- lrm(Group ~differentiation.grade +CRC.site+Pathological.Type+Perineural.invasion+
                      Family_cancer_history+CEA..5+Early_Onset+Vascular.invasion,data = train_set, x = TRUE, y = TRUE)

lrm_model_all <- lrm(Group ~.,data = train_set, x = TRUE, y = TRUE)

# calibrate residuals analysis
cal_data_all <- calibrate(lrm_model_all, method = "boot", B = 1000)
cal_data_11 <- calibrate(lrm_model_11, method = "boot", B = 1000)
cal_data_6 <- calibrate(lrm_model_6, method = "boot", B = 1000)

# plot calibrate residuals for model
plot(cal_data_all)
plot(cal_data_11)
plot(cal_data_6)

# Predicted in test dataset
lrm_model_all <- predict(lrm_model_all, newdata=test_set)
lrm_model_11 <- predict(lrm_model_11, newdata=test_set)
lrm_model_6 <- predict(lrm_model_6, newdata=test_set)

roc_obj1 <- roc(test_set$Group, as.numeric(lrm_model_all))
roc_obj2 <- roc(test_set$Group, as.numeric(lrm_model_11))
roc_obj3 <- roc(test_set$Group, as.numeric(lrm_model_6))

# plot ROC curve
plot.roc(roc_obj1, col="blue")
lines.roc(roc_obj2, col="red")
lines.roc(roc_obj3, col="green")
legend("bottomright", legend=c("all-features", "11-features", "6-features"), col=c("blue", "red", "green"), lty=1)

############### nomogram construction
nom <- nomogram(lrm_model_11, fun=plogis)
plot(nom)




