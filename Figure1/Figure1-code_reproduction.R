#############################
#############################
########packages loading
library(survival)
library(survminer)

######## survival data loading
load("Figure1BC-survive-data.RData")     

######## Creating survival object for OS/DFS
surv_object_OS <- Surv(survival_data_clinic$OS, survival_data_clinic$Dead)
surv_object_DFS <- Surv(survival_data_clinic$DFS, survival_data_clinic$Recurrent)

######## Creating a survival curves for OS/DFS
fit_OS <- survfit(surv_object_OS ~ Group,data=survival_data_clinic)
fit_DFS <- survfit(surv_object_DFS ~ Group,data=survival_data_clinic)

######## Getting 5 years old OS & DFS
summary(fit_OS, times = 60)     #获得5年生存率
summary(fit_DFS, times = 60)

######## Figure1B_5yrs_DFS
plot_survDFS <- ggsurvplot(fit_DFS, data = survival_data_clinic,
                           xlim = c(0, 60),
                           xlab = "Month", 
                           pval = TRUE,
                           pval.method = TRUE,
                           ylab = "Disease-Free Survival",  

                           break.time.by = 12, 
                           risk.table = "abs_pct",  
                           risk.table.y.text.col = T,
                           risk.table.y.text = FALSE,
                           ggtheme = theme_survminer())

print(plot_survDFS)

######## Figure1B_5yrs_OS
plot_survDFS <- ggsurvplot(fit_OS, data = survival_data_clinic,
                           xlim = c(0, 60),
                           xlab = "Month", 
                           pval = TRUE,
                           pval.method = TRUE,
                           ylab = "Overall survival",  
                           
                           break.time.by = 12, 
                           risk.table = "abs_pct",  
                           risk.table.y.text.col = T,
                           risk.table.y.text = FALSE,
                           ggtheme = theme_survminer())

print(plot_survDFS)

###########################
###########################

######## paired log-rank p-value
group_levels <- unique(survival_data_clinic$Group)

######## initialization
pairwise_p_values <- matrix(NA, nrow = length(group_levels), ncol = length(group_levels), 
                            dimnames = list(group_levels, group_levels))

######## Paried Log-rank test
for(i in 1:(length(group_levels) - 1)) {
  for(j in (i + 1):length(group_levels)) {
    subset_data <- survival_data_clinic[survival_data_clinic$Group %in% c(group_levels[i], group_levels[j]),]

    test <- survdiff(Surv(DFS, Recurrent) ~ Group, data = subset_data)  ## DFS
    #test <- survdiff(Surv(OS, Dead) ~ Group, data = subset_data)  ## OS
    pairwise_p_values[i, j] <- 1 - pchisq(test$chisq, length(test$n) - 1)
    pairwise_p_values[j, i] <- pairwise_p_values[i, j] 
  }
}
# Print the matrix of log-rank p-value
print(pairwise_p_values)

###########################
###########################
######## HR analysis
## transform OS/DFS datato 5 years old
survival_data_clinic$status_5yr <- ifelse(survival_data_clinic$OS >= 60, 0, survival_data_clinic$OS)
survival_data_clinic$statusDFS_5yr <- ifelse(survival_data_clinic$DFS >= 60, 0, survival_data_clinic$DFS)
surv_object_OS <- Surv(survival_data_clinic$status_5yr, survival_data_clinic$Dead)
surv_object_DFS <- Surv(survival_data_clinic$statusDFS_5yr, survival_data_clinic$Recurrent)

## set reference: pMMR-MSS
survival_data_clinic$Group <- factor(survival_data_clinic$Group,levels = c("pMMR_MSS","pMMR_MSI","dMMR_MSI"))
survival_data_clinic$Group <- relevel(survival_data_clinic$Group,ref="pMMR_MSS")

## cox proportional hazards regression model analysis
## DFS analysis
cox_model_DFS <- coxph(surv_object_DFS ~ Group, data = survival_data_clinic) 
## OS analysis
cox_model_OS <- coxph(surv_object_OS ~ Group, data = survival_data_clinic)

## get the odds ratio and 95%CI
summary(cox_model_DFS)     
summary(cox_model_OS)     

######################## 
######################## estimate the residual, and cox regression assumption analysis
## get schoenfeld_residuals
DFS_residuals <- cox.zph(cox_model_DFS)
OS_residuals <- cox.zph(cox_model_OS)

plot(DFS_residuals)
plot(OS_residuals)

## martingale residuals analysis
group_numeric <- as.numeric(as.factor(survival_data_clinic$Group))

## Assuming cox_model is your Cox regression model
DFS_mar_residuals <- residuals(cox_model_DFS, type = "martingale")
OS_mar_residuals <- residuals(cox_model_OS, type = "martingale")

# If there are NA values, you might want to exclude them
valid_indices <- !is.na(group_numeric) & !is.na(DFS_mar_residuals)
valid_indices <- !is.na(group_numeric) & !is.na(OS_mar_residuals)

group_numeric <- group_numeric[valid_indices]

DFS_mar_residuals <- DFS_mar_residuals[valid_indices]
OS_mar_residuals <- OS_mar_residuals[valid_indices]

# If Group_msi2 is a factor, consider a boxplot to see the distribution of residuals by group
boxplot(DFS_mar_residuals ~ survival_data_clinic$Group)
boxplot(OS_mar_residuals ~ survival_data_clinic$Group)










