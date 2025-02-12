# ============================================================================ #
# File: TimeAnalysis_ADNI.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Perform survival (time) analysis in the ADNI cohort.
# ============================================================================ #

# Load packages
library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(tidyverse)
library(caret)
library(patchwork)
library(ranger)
library(pROC)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("ADNI/Data/predictedScore_factors_ADNI.RData")
load("ADNI/Data/MetaData_ADNI.RData")

# Filter for midlife samples
metaData_fil <- as.data.frame(MetaData_baseline)
rownames(metaData_fil) <- metaData_fil$Basename
midlife_samples <- intersect(metaData_fil$Basename[metaData_fil$Age <= 75], 
                             rownames(predictedScore_factors))
metaData_fil <- metaData_fil[midlife_samples,]

# prepare data
predictedScore_factors_fil <- predictedScore_factors[metaData_fil$Basename,]
table(metaData_fil$DX)

# Meta data for all time points
metaData_fil_time <- MetaData_allTime[MetaData_allTime$RID %in% metaData_fil$RID,]

# Make predictions using epi-MCI (RF-RFE) model
load("Models/EMIF_Models/MRS/Fit_EMIF_MCI_RF.RData")
pred_RF <- predict(fit, predictedScore_factors_fil, type = "prob")
predictDF <- data.frame(RID = metaData_fil$RID, 
                        pred = pred_RF$MCI,
                        Age = metaData_fil$Age,
                        Sex = metaData_fil$Sex)

# Get mean prediction for same individual 
# (some individuals have more than one baseline sample available)
predictDF <- predictDF %>%
  group_by(RID) %>%
  reframe(RID = RID,
          pred= mean(pred),
          Age = Age,
          Sex = Sex)

predictDF <- unique(predictDF)


# Epi-MCI (RF-RFE) score as predictor
predictDF$RID <- as.character(predictDF$RID)
predictDF$predClass <- "Intermediate risk"
predictDF$predClass[predictDF$pred < quantile(predictDF$pred,0.33)] <- "Low risk"
predictDF$predClass[predictDF$pred > quantile(predictDF$pred,0.67)] <- "High risk"
table(predictDF$predClass)


# Histogram of predictions
predictDF$predClass <- factor(predictDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
colors <- c("#FDD0A2","#FD8D3C","#D94801")
colors <- c("#4292C6","grey","#EF3B2C")
p <- ggplot(predictDF) + 
  geom_histogram(aes(x = log(pred/(1-pred)), fill = predClass),
                 position = "identity",color = "black", bins= 35) +
  theme_classic() +
  xlab("Epi-MCI Score") +
  ylab("Count") +
  scale_fill_manual(values = colors) +
  theme(legend.position = "bottom",
        legend.title = element_blank())

# Save plot
ggsave(p, file = "ADNI/TimeAnalysis/RiskScoreDistribution_ADNI.png", width = 7, height = 5)


###############################################################################

# MMSE

###############################################################################

# Which variable
var <- "MMSE"

# Remove missing values
metaData_fil_var <- metaData_fil_time[!is.na(metaData_fil_time[,var]),]
metaData_fil_var <- metaData_fil_var[!is.na(metaData_fil_var[,paste0(var,".bl")]),]

# Set sex
metaData_fil_var$Sex <- ifelse(metaData_fil_var$PTGENDER == "Male",0,1)

# Prepare data for survival analysis
testDF <- metaData_fil_var[,c("RID", "VISCODE", var, paste0(var,".bl"), "AGE", "Sex")]

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Set cognitive impairment status
testDF <- testDF[testDF[,paste0(var,".bl")]>=26,]
testDF$Status <- ifelse((testDF[,var] < 24), "Cognitive Impaired", "Normal")

# Only keep samples with diagnosis available at last or second-to-last time point
keep_samples <- unique(testDF$RID[(testDF$Status == "Cognitive Impaired") & testDF$Time <= 48])
keep_samples <- unique(c(keep_samples, 
                         as.character(testDF$RID[(testDF$VISCODE == "m48") | (testDF$VISCODE == "m42")])))


testDF <- testDF[testDF$RID %in% keep_samples,]

# Only include cognitive impaired
testDF_imp <- testDF[testDF$Status == "Cognitive Impaired",]

# Only include up to 48 months in analysis
testDF_imp <- testDF_imp[testDF_imp$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF_imp$RID)){
  testDF_imp[testDF_imp$RID == i, "Time"] <- min(testDF_imp[testDF_imp$RID == i, "Time"])
}

testDF_imp <- testDF_imp[!duplicated(testDF_imp[,c(1,7,8)]),]
colnames(testDF_imp) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")

# Set the samples not converted to a cognitively healthy status at month 49
testDF_normal <- testDF[!(testDF$RID %in% testDF_imp$RID),]
colnames(testDF_normal) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")
testDF_normal$Time <- 49
testDF_normal$VISCODE <- "m49"
testDF_normal$Var <- NA
testDF_normal <- testDF_normal[!duplicated(testDF_normal$RID),]


testDF_all <- rbind.data.frame(testDF_imp, testDF_normal)
testDF_all$RID <- as.character(testDF_all$RID)
length(unique(testDF_all$RID)) == length(keep_samples)

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF_all[,c("RID", "Time", "Status", "Age", "Sex")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF[,c("RID", "pred", "predClass")], by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
table(kaplanDF$predClass)

colors <- c("#FDD0A2","#FD8D3C","#D94801")
colors <- c("#4292C6","grey","#EF3B2C")
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("MMSE") +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = "ADNI/TimeAnalysis/KaplanMeier_ADNI_MMSE1.jpg", width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF
kaplanDF1$predClass <- factor(ifelse(kaplanDF1$predClass == "High risk", "High risk", "Low risk"), 
                              levels = c("Low risk", "High risk"))

kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass + Age + Sex, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# RAVLT learning

###############################################################################

# Which variable
var <- "RAVLT.learning"

# Remove missing values
metaData_fil_var <- metaData_fil_time[!is.na(metaData_fil_time[,var]),]
metaData_fil_var <- metaData_fil_var[!is.na(metaData_fil_var[,paste0(var,".bl")]),]

# Set sex
metaData_fil_var$Sex <- ifelse(metaData_fil_var$PTGENDER == "Male",0,1)

# Prepare data for survival analysis
testDF <- metaData_fil_var[,c("RID", "VISCODE", var, paste0(var,".bl"), "AGE", "Sex")]

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Set cognitive impairment status
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,paste0(var,".bl")] >= mean_var - 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] < mean_var - 2*sd_var), "Cognitive Impaired", "Normal")

# Only keep samples with diagnosis available at last or second-to-last time point
keep_samples <- unique(testDF$RID[(testDF$Status == "Cognitive Impaired") & testDF$Time <= 48])
keep_samples <- unique(c(keep_samples, 
                         as.character(testDF$RID[(testDF$VISCODE == "m48") | (testDF$VISCODE == "m42")])))


testDF <- testDF[testDF$RID %in% keep_samples,]

# Only include cognitive impaired
testDF_imp <- testDF[testDF$Status == "Cognitive Impaired",]

# Only include up to 48 months in analysis
testDF_imp <- testDF_imp[testDF_imp$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF_imp$RID)){
  testDF_imp[testDF_imp$RID == i, "Time"] <- min(testDF_imp[testDF_imp$RID == i, "Time"])
}

testDF_imp <- testDF_imp[!duplicated(testDF_imp[,c(1,7,8)]),]
colnames(testDF_imp) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")

# Set the samples not converted to a cognitively healthy status at month 49
testDF_normal <- testDF[!(testDF$RID %in% testDF_imp$RID),]
colnames(testDF_normal) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")
testDF_normal$Time <- 49
testDF_normal$VISCODE <- "m49"
testDF_normal$Var <- NA
testDF_normal <- testDF_normal[!duplicated(testDF_normal$RID),]


testDF_all <- rbind.data.frame(testDF_imp, testDF_normal)
testDF_all$RID <- as.character(testDF_all$RID)
length(unique(testDF_all$RID)) == length(keep_samples)

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF_all[,c("RID", "Time", "Status", "Age", "Sex")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF[,c("RID", "pred", "predClass")], by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
table(kaplanDF$predClass)

colors <- c("#FDD0A2","#FD8D3C","#D94801")
colors <- c("#4292C6","grey","#EF3B2C")
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("RAVLT (learning)") +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,"1.jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF
kaplanDF1$predClass <- factor(ifelse(kaplanDF1$predClass == "High risk", "High risk", "Low risk"), 
                              levels = c("Low risk", "High risk"))

kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass + Age + Sex, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# RAVLT forgetting

###############################################################################

# Which variable
var <- "RAVLT.forgetting"

# Remove missing values
metaData_fil_var <- metaData_fil_time[!is.na(metaData_fil_time[,var]),]
metaData_fil_var <- metaData_fil_var[!is.na(metaData_fil_var[,paste0(var,".bl")]),]

# Set sex
metaData_fil_var$Sex <- ifelse(metaData_fil_var$PTGENDER == "Male",0,1)

# Prepare data for survival analysis
testDF <- metaData_fil_var[,c("RID", "VISCODE", var, paste0(var,".bl"), "AGE", "Sex")]

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Set cognitive impairment status
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,paste0(var,".bl")] <= mean_var + 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

# Only keep samples with diagnosis available at last or second-to-last time point
keep_samples <- unique(testDF$RID[(testDF$Status == "Cognitive Impaired") & testDF$Time <= 48])
keep_samples <- unique(c(keep_samples, 
                         as.character(testDF$RID[(testDF$VISCODE == "m48") | (testDF$VISCODE == "m42")])))


testDF <- testDF[testDF$RID %in% keep_samples,]

# Only include cognitive impaired
testDF_imp <- testDF[testDF$Status == "Cognitive Impaired",]

# Only include up to 48 months in analysis
testDF_imp <- testDF_imp[testDF_imp$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF_imp$RID)){
  testDF_imp[testDF_imp$RID == i, "Time"] <- min(testDF_imp[testDF_imp$RID == i, "Time"])
}

testDF_imp <- testDF_imp[!duplicated(testDF_imp[,c(1,7,8)]),]
colnames(testDF_imp) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")

# Set the samples not converted to a cognitively healthy status at month 49
testDF_normal <- testDF[!(testDF$RID %in% testDF_imp$RID),]
colnames(testDF_normal) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")
testDF_normal$Time <- 49
testDF_normal$VISCODE <- "m49"
testDF_normal$Var <- NA
testDF_normal <- testDF_normal[!duplicated(testDF_normal$RID),]


testDF_all <- rbind.data.frame(testDF_imp, testDF_normal)
testDF_all$RID <- as.character(testDF_all$RID)
length(unique(testDF_all$RID)) == length(keep_samples)

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF_all[,c("RID", "Time", "Status", "Age", "Sex")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF[,c("RID", "pred", "predClass")], by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
table(kaplanDF$predClass)

colors <- c("#FDD0A2","#FD8D3C","#D94801")
colors <- c("#FDD0A2","#F16913","#8C2D04")
colors <- c("#4292C6","grey","#EF3B2C")
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("RAVLT (forgetting)") +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,"1.jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass + Age + Sex, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# RAVLT percent forgetting

###############################################################################

# Which variable
var <- "RAVLT.perc.forgetting"

# Remove missing values
metaData_fil_var <- metaData_fil_time[!is.na(metaData_fil_time[,var]),]
metaData_fil_var <- metaData_fil_var[!is.na(metaData_fil_var[,paste0(var,".bl")]),]

# Set sex
metaData_fil_var$Sex <- ifelse(metaData_fil_var$PTGENDER == "Male",0,1)

# Prepare data for survival analysis
testDF <- metaData_fil_var[,c("RID", "VISCODE", var, paste0(var,".bl"), "AGE", "Sex")]

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Set cognitive impairment status
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,paste0(var,".bl")] <= mean_var + 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

# Only keep samples with diagnosis available at last or second-to-last time point
keep_samples <- unique(testDF$RID[(testDF$Status == "Cognitive Impaired") & testDF$Time <= 48])
keep_samples <- unique(c(keep_samples, 
                         as.character(testDF$RID[(testDF$VISCODE == "m48") | (testDF$VISCODE == "m42")])))


testDF <- testDF[testDF$RID %in% keep_samples,]

# Only include cognitive impaired
testDF_imp <- testDF[testDF$Status == "Cognitive Impaired",]

# Only include up to 48 months in analysis
testDF_imp <- testDF_imp[testDF_imp$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF_imp$RID)){
  testDF_imp[testDF_imp$RID == i, "Time"] <- min(testDF_imp[testDF_imp$RID == i, "Time"])
}

testDF_imp <- testDF_imp[!duplicated(testDF_imp[,c(1,7,8)]),]
colnames(testDF_imp) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")

# Set the samples not converted to a cognitively healthy status at month 49
testDF_normal <- testDF[!(testDF$RID %in% testDF_imp$RID),]
colnames(testDF_normal) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")
testDF_normal$Time <- 49
testDF_normal$VISCODE <- "m49"
testDF_normal$Var <- NA
testDF_normal <- testDF_normal[!duplicated(testDF_normal$RID),]


testDF_all <- rbind.data.frame(testDF_imp, testDF_normal)
testDF_all$RID <- as.character(testDF_all$RID)
length(unique(testDF_all$RID)) == length(keep_samples)

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF_all[,c("RID", "Time", "Status", "Age", "Sex")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF[,c("RID", "pred", "predClass")], by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
table(kaplanDF$predClass)

colors <- c("#FDD0A2","#FD8D3C","#D94801")
colors <- c("#FDD0A2","#F16913","#8C2D04")
colors <- c("#4292C6","grey","#EF3B2C")

p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("RAVLT (percent forgetting)") +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,"1.jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass + Age + Sex, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)

###############################################################################

# RAVLT immediate

###############################################################################

# Which variable
var <- "RAVLT.immediate"

# Remove missing values
metaData_fil_var <- metaData_fil_time[!is.na(metaData_fil_time[,var]),]
metaData_fil_var <- metaData_fil_var[!is.na(metaData_fil_var[,paste0(var,".bl")]),]

# Set sex
metaData_fil_var$Sex <- ifelse(metaData_fil_var$PTGENDER == "Male",0,1)

# Prepare data for survival analysis
testDF <- metaData_fil_var[,c("RID", "VISCODE", var, paste0(var,".bl"), "AGE", "Sex")]

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Set cognitive impairment status
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,paste0(var,".bl")] >= mean_var - 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] < mean_var - 2*sd_var), "Cognitive Impaired", "Normal")

# Only keep samples with diagnosis available at last or second-to-last time point
keep_samples <- unique(testDF$RID[(testDF$Status == "Cognitive Impaired") & testDF$Time <= 48])
keep_samples <- unique(c(keep_samples, 
                         as.character(testDF$RID[(testDF$VISCODE == "m48") | (testDF$VISCODE == "m42")])))


testDF <- testDF[testDF$RID %in% keep_samples,]

# Only include cognitive impaired
testDF_imp <- testDF[testDF$Status == "Cognitive Impaired",]

# Only include up to 48 months in analysis
testDF_imp <- testDF_imp[testDF_imp$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF_imp$RID)){
  testDF_imp[testDF_imp$RID == i, "Time"] <- min(testDF_imp[testDF_imp$RID == i, "Time"])
}

testDF_imp <- testDF_imp[!duplicated(testDF_imp[,c(1,7,8)]),]
colnames(testDF_imp) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")

# Set the samples not converted to a cognitively healthy status at month 49
testDF_normal <- testDF[!(testDF$RID %in% testDF_imp$RID),]
colnames(testDF_normal) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")
testDF_normal$Time <- 49
testDF_normal$VISCODE <- "m49"
testDF_normal$Var <- NA
testDF_normal <- testDF_normal[!duplicated(testDF_normal$RID),]


testDF_all <- rbind.data.frame(testDF_imp, testDF_normal)
testDF_all$RID <- as.character(testDF_all$RID)
length(unique(testDF_all$RID)) == length(keep_samples)

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF_all[,c("RID", "Time", "Status", "Age", "Sex")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF[,c("RID", "pred", "predClass")], by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
table(kaplanDF$predClass)

colors <- c("#FDD0A2","#FD8D3C","#D94801")
colors <- c("#FDD0A2","#F16913","#8C2D04")
colors <- c("#4292C6","grey","#EF3B2C")
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("RAVLT (immediate)") +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,"1.jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass + Age + Sex, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# ADASQ4

###############################################################################

# Which variable
var <- "ADASQ4"

# Remove missing values
metaData_fil_var <- metaData_fil_time[!is.na(metaData_fil_time[,var]),]
metaData_fil_var <- metaData_fil_var[!is.na(metaData_fil_var[,paste0(var,".bl")]),]

# Set sex
metaData_fil_var$Sex <- ifelse(metaData_fil_var$PTGENDER == "Male",0,1)

# Prepare data for survival analysis
testDF <- metaData_fil_var[,c("RID", "VISCODE", var, paste0(var,".bl"), "AGE", "Sex")]

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Set cognitive impairment status
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,paste0(var,".bl")] <= mean_var + 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

# Only keep samples with diagnosis available at last or second-to-last time point
keep_samples <- unique(testDF$RID[(testDF$Status == "Cognitive Impaired") & testDF$Time <= 48])
keep_samples <- unique(c(keep_samples, 
                         as.character(testDF$RID[(testDF$VISCODE == "m48") | (testDF$VISCODE == "m42")])))


testDF <- testDF[testDF$RID %in% keep_samples,]

# Only include cognitive impaired
testDF_imp <- testDF[testDF$Status == "Cognitive Impaired",]

# Only include up to 48 months in analysis
testDF_imp <- testDF_imp[testDF_imp$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF_imp$RID)){
  testDF_imp[testDF_imp$RID == i, "Time"] <- min(testDF_imp[testDF_imp$RID == i, "Time"])
}

testDF_imp <- testDF_imp[!duplicated(testDF_imp[,c(1,7,8)]),]
colnames(testDF_imp) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")

# Set the samples not converted to a cognitively healthy status at month 49
testDF_normal <- testDF[!(testDF$RID %in% testDF_imp$RID),]
colnames(testDF_normal) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")
testDF_normal$Time <- 49
testDF_normal$VISCODE <- "m49"
testDF_normal$Var <- NA
testDF_normal <- testDF_normal[!duplicated(testDF_normal$RID),]


testDF_all <- rbind.data.frame(testDF_imp, testDF_normal)
testDF_all$RID <- as.character(testDF_all$RID)
length(unique(testDF_all$RID)) == length(keep_samples)

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF_all[,c("RID", "Time", "Status", "Age", "Sex")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF[,c("RID", "pred", "predClass")], by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
table(kaplanDF$predClass)

colors <- c("#FDD0A2","#FD8D3C","#D94801")
colors <- c("#FDD0A2","#F16913","#8C2D04")
colors <- c("#4292C6","grey","#EF3B2C")
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("ADAS-Q4") +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,"1.jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass + Age + Sex, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# ADAS13

###############################################################################

# Which variable
var <- "ADAS13"

# Remove missing values
metaData_fil_var <- metaData_fil_time[!is.na(metaData_fil_time[,var]),]
metaData_fil_var <- metaData_fil_var[!is.na(metaData_fil_var[,paste0(var,".bl")]),]

# Set sex
metaData_fil_var$Sex <- ifelse(metaData_fil_var$PTGENDER == "Male",0,1)

# Prepare data for survival analysis
testDF <- metaData_fil_var[,c("RID", "VISCODE", var, paste0(var,".bl"), "AGE", "Sex")]

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Set cognitive impairment status
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,paste0(var,".bl")] <= mean_var + 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

# Only keep samples with diagnosis available at last or second-to-last time point
keep_samples <- unique(testDF$RID[(testDF$Status == "Cognitive Impaired") & testDF$Time <= 48])
keep_samples <- unique(c(keep_samples, 
                         as.character(testDF$RID[(testDF$VISCODE == "m48") | (testDF$VISCODE == "m42")])))


testDF <- testDF[testDF$RID %in% keep_samples,]

# Only include cognitive impaired
testDF_imp <- testDF[testDF$Status == "Cognitive Impaired",]

# Only include up to 48 months in analysis
testDF_imp <- testDF_imp[testDF_imp$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF_imp$RID)){
  testDF_imp[testDF_imp$RID == i, "Time"] <- min(testDF_imp[testDF_imp$RID == i, "Time"])
}

testDF_imp <- testDF_imp[!duplicated(testDF_imp[,c(1,7,8)]),]
colnames(testDF_imp) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")

# Set the samples not converted to a cognitively healthy status at month 49
testDF_normal <- testDF[!(testDF$RID %in% testDF_imp$RID),]
colnames(testDF_normal) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")
testDF_normal$Time <- 49
testDF_normal$VISCODE <- "m49"
testDF_normal$Var <- NA
testDF_normal <- testDF_normal[!duplicated(testDF_normal$RID),]


testDF_all <- rbind.data.frame(testDF_imp, testDF_normal)
testDF_all$RID <- as.character(testDF_all$RID)
length(unique(testDF_all$RID)) == length(keep_samples)

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF_all[,c("RID", "Time", "Status", "Age", "Sex")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF[,c("RID", "pred", "predClass")], by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
table(kaplanDF$predClass)

colors <- c("#FDD0A2","#FD8D3C","#D94801")
colors <- c("#FDD0A2","#F16913","#8C2D04")
colors <- c("#4292C6","grey","#EF3B2C")
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle(var) +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,"1.jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass + Age + Sex, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# ADAS11

###############################################################################

# Which variable
var <- "ADAS11"

# Remove missing values
metaData_fil_var <- metaData_fil_time[!is.na(metaData_fil_time[,var]),]
metaData_fil_var <- metaData_fil_var[!is.na(metaData_fil_var[,paste0(var,".bl")]),]

# Set sex
metaData_fil_var$Sex <- ifelse(metaData_fil_var$PTGENDER == "Male",0,1)

# Prepare data for survival analysis
testDF <- metaData_fil_var[,c("RID", "VISCODE", var, paste0(var,".bl"), "AGE", "Sex")]

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Set cognitive impairment status
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,paste0(var,".bl")] <= mean_var + 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

# Only keep samples with diagnosis available at last or second-to-last time point
keep_samples <- unique(testDF$RID[(testDF$Status == "Cognitive Impaired") & testDF$Time <= 48])
keep_samples <- unique(c(keep_samples, 
                         as.character(testDF$RID[(testDF$VISCODE == "m48") | (testDF$VISCODE == "m42")])))


testDF <- testDF[testDF$RID %in% keep_samples,]

# Only include cognitive impaired
testDF_imp <- testDF[testDF$Status == "Cognitive Impaired",]

# Only include up to 48 months in analysis
testDF_imp <- testDF_imp[testDF_imp$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF_imp$RID)){
  testDF_imp[testDF_imp$RID == i, "Time"] <- min(testDF_imp[testDF_imp$RID == i, "Time"])
}

testDF_imp <- testDF_imp[!duplicated(testDF_imp[,c(1,7,8)]),]
colnames(testDF_imp) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")

# Set the samples not converted to a cognitively healthy status at month 49
testDF_normal <- testDF[!(testDF$RID %in% testDF_imp$RID),]
colnames(testDF_normal) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")
testDF_normal$Time <- 49
testDF_normal$VISCODE <- "m49"
testDF_normal$Var <- NA
testDF_normal <- testDF_normal[!duplicated(testDF_normal$RID),]


testDF_all <- rbind.data.frame(testDF_imp, testDF_normal)
testDF_all$RID <- as.character(testDF_all$RID)
length(unique(testDF_all$RID)) == length(keep_samples)

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF_all[,c("RID", "Time", "Status", "Age", "Sex")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF[,c("RID", "pred", "predClass")], by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
table(kaplanDF$predClass)

colors <- c("#FDD0A2","#FD8D3C","#D94801")
colors <- c("#FDD0A2","#F16913","#8C2D04")
colors <- c("#4292C6","grey","#EF3B2C")
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle(var) +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,"1.jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass + Age + Sex, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# Log rank test
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)


###############################################################################

# TRABSCOR: Trails B

###############################################################################

# Which variable
var <- "TRABSCOR"

# Remove missing values
metaData_fil_var <- metaData_fil_time[!is.na(metaData_fil_time[,var]),]
metaData_fil_var <- metaData_fil_var[!is.na(metaData_fil_var[,paste0(var,".bl")]),]

# Set sex
metaData_fil_var$Sex <- ifelse(metaData_fil_var$PTGENDER == "Male",0,1)

# Prepare data for survival analysis
testDF <- metaData_fil_var[,c("RID", "VISCODE", var, paste0(var,".bl"), "AGE", "Sex")]

# Make time point numeric
testDF$VISCODE[testDF$VISCODE == "bl"] <- 0
testDF$Time <- as.numeric(str_remove(testDF$VISCODE, "m"))
testDF$Time <- as.numeric(testDF$Time)

# Set cognitive impairment status
sd_var <- sd(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
mean_var <- mean(metaData_fil[metaData_fil$DX == "CN",var], na.rm = TRUE)
testDF <- testDF[testDF[,paste0(var,".bl")] <= mean_var + 2*sd_var,]
testDF$Status <- ifelse((testDF[,var] > mean_var + 2*sd_var), "Cognitive Impaired", "Normal")

# Only keep samples with diagnosis available at last or second-to-last time point
keep_samples <- unique(testDF$RID[(testDF$Status == "Cognitive Impaired") & testDF$Time <= 48])
keep_samples <- unique(c(keep_samples, 
                         as.character(testDF$RID[(testDF$VISCODE == "m48") | (testDF$VISCODE == "m42")])))


testDF <- testDF[testDF$RID %in% keep_samples,]

# Only include cognitive impaired
testDF_imp <- testDF[testDF$Status == "Cognitive Impaired",]

# Only include up to 48 months in analysis
testDF_imp <- testDF_imp[testDF_imp$Time <= 48,]

# Set time to earliest time point converted to cognitive impairments
for (i in unique(testDF_imp$RID)){
  testDF_imp[testDF_imp$RID == i, "Time"] <- min(testDF_imp[testDF_imp$RID == i, "Time"])
}

testDF_imp <- testDF_imp[!duplicated(testDF_imp[,c(1,7,8)]),]
colnames(testDF_imp) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")

# Set the samples not converted to a cognitively healthy status at month 49
testDF_normal <- testDF[!(testDF$RID %in% testDF_imp$RID),]
colnames(testDF_normal) <- c("RID", "VISCODE", "Var"," Var.bl","Age","Sex","Time","Status")
testDF_normal$Time <- 49
testDF_normal$VISCODE <- "m49"
testDF_normal$Var <- NA
testDF_normal <- testDF_normal[!duplicated(testDF_normal$RID),]


testDF_all <- rbind.data.frame(testDF_imp, testDF_normal)
testDF_all$RID <- as.character(testDF_all$RID)
length(unique(testDF_all$RID)) == length(keep_samples)

# Put data into correct format:
# 1: censored
# 2: disease
kaplanDF <- testDF_all[,c("RID", "Time", "Status", "Age", "Sex")]
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF[,c("RID", "pred", "predClass")], by = c("RID" = "RID"))
kaplanDF$Time <- as.numeric(kaplanDF$Time)/12

# Make survival curve
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low risk", "Intermediate risk", "High risk"))
table(kaplanDF$predClass)

colors <- c("#FDD0A2","#FD8D3C","#D94801")
colors <- c("#FDD0A2","#F16913","#8C2D04")
colors <- c("#4292C6","grey","#EF3B2C")
p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval(alpha = 0.15) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Follow-up time (years)") +
  #scale_x_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("Trail Making Test Part B Time") +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("ADNI/TimeAnalysis/KaplanMeier_ADNI_",var,"1.jpg"), width = 7, height = 5)

# Compare low and high risk
kaplanDF1 <- kaplanDF[kaplanDF$predClass != "Intermediate risk",]
kaplanDF1$predClass <- factor(kaplanDF1$predClass, levels = c("Low risk", "High risk"))

# Log rank test
survdiff(Surv(Time, Test) ~ predClass, 
         data = kaplanDF1)

# Cox regression
coxph(Surv(Time, Test) ~ predClass + Sex + Age, 
      data = kaplanDF1) %>% 
  tbl_regression(exp = TRUE) 

# AUROC
kaplanDF1$Status <- factor(kaplanDF1$Status, levels = c("Normal", "Cognitive Impaired"))
test_roc <- pROC::roc(kaplanDF1$Status, kaplanDF1$pred, direction = "<")
auc(test_roc)
ci(test_roc)