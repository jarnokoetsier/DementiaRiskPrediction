# ============================================================================ #
# File: ShapValues.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Calculate SHAP values of Epi-MCI model's features.
# ============================================================================ #

#######################################################################################

# SHAP values: MRSs

#######################################################################################

# Load packages
library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library("Numero")
library("DALEX") # SHAP values
library(fmsb) # Radar Chart

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF-AD/Data/X_train_EMIF.RData")
load("EMIF-AD/Data/Y_train_EMIF.RData")
load("EMIF-AD/Data/X_test_EMIF.RData")
load("EMIF-AD/Data/Y_test_EMIF.RData")


###############################################################################

# Get SHAP values

###############################################################################

# Load MCI model
load("Models/EMIF_Models/MRS/Fit_EMIF_MCI_sPLS.RData")

# Make explained object
explainer <- DALEX::explain(fit, data = X_train, y = NULL)

# Calculate SHAP values
for (i in 1:nrow(X_test)){
  shap <- predict_parts(explainer = explainer, 
                        new_observation = X_test[i,], 
                        type = "shap",
                        B = 50)
  
  # get mean contribution (i.e., SHAP value)
  shap_fil <- shap %>%
    group_by(variable_name) %>%
    reframe(
      Contribution = mean(contribution))
  
  if (i == 1){
    output <- shap_fil
  } else{
    output <- inner_join(output, shap_fil, by = c("variable_name" = "variable_name"))
  }
  
}
output <- as.data.frame(output)
rownames(output) <- output$variable_name
output <- output[,-1]
colnames(output) <- rownames(X_test)

# Save SHAP values 
save(output, file = "Models/ModelInterpretation/SHAP/output_shap_sPLS.RData")


###############################################################################

# Compare SHAP values between different models

###############################################################################

#==============================================================================#
# MCI models from EMIF
#==============================================================================#

# Models
models <- c("EN", "sPLS", "RF")

plotDF <- NULL
for (i in 1:length(models)){
  
  # Load SHAP values
  load(paste0("Models/ModelInterpretation/SHAP/output_shap_", models[i], ".RData"))
  
  # Scale SHAP values per individual
  #output_scaled <- t(t(output)/colSums(abs(output)))
  output_scaled <- output
  
  # Get average SHAP values
  temp <- data.frame(Variable = rownames(output_scaled),
                        value = rowMeans(abs(output_scaled)))
  
  if (i > 1){
    plotDF <- inner_join(plotDF, temp, by = c("Variable" = "Variable"))
  } else{
    plotDF <- temp
  }
}

# Set row names
rownames(plotDF) <- c("Alcohol Intake", "BMI", "Depression", "Type II Diabetes",
                      "Unhealthy Diet", "Low Education", "Age","HDL Chol.", "Heart Disease",
                      "Physical Inact.", "Sex", "Smoking", "Syst. Blood Pressure", "Total Chol.")
plotDF <- plotDF[,-1]

# Set column names
colnames(plotDF) <- c("ElasticNet", "sPLS", "Random Forest")

# Scale SHAP values
for (i in 1:ncol(plotDF)){
  plotDF[,i] <- plotDF[,i]/sum(plotDF[,i])
}


# Prepare data for plotting
data <- as.data.frame(t(plotDF))
data <- rbind(rep(0.3,14) , rep(0,14) , data)
data <- data[,c("Alcohol Intake", "BMI", "HDL Chol.", "Type II Diabetes",
                "Unhealthy Diet", "Low Education", "Age","Total Chol.", "Heart Disease",
                "Physical Inact.", "Sex", "Smoking", "Syst. Blood Pressure", "Depression")]

# Color vector
colors_border=rev(c("#EF3B2C","#CB181D", "#99000D") )
colors_border=rev(c("#FCBBA1","#EF3B2C", "#99000D") )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )

# plot with default options:
jpeg("Models/ModelInterpretation/SHAP/shaply_emif_new1.jpg", width = 8000, height = 7000, quality = 100)
radarchart( data  , axistype=0 , 
            #custom polygon
            pcol=colors_border, plwd=100 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,0.1,0.05), cglwd=8,
            #custom labels
            vlcex=20 
)
dev.off()



#######################################################################################

# SHAP values: MRSs + CSF biomarkers (NOT TESTED YET)

#######################################################################################

# Load packages
library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library("Numero")
library("DALEX") # SHAP values
library(fmsb) # Radar Chart

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF-AD/Data/X_train_EMIF.RData")
load("EMIF-AD/Data/Y_train_EMIF.RData")
load("EMIF-AD/Data/X_test_EMIF.RData")
load("EMIF-AD/Data/Y_test_EMIF.RData")
load("EMIF-AD/Data/metaData_fil.RData")

# Prepare data
rownames(metaData_fil) <- metaData_fil$X
CSFbio <- metaData_fil[,c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore", "Age")]
colnames(CSFbio) <- c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore", "ChrAge")
samples <- rownames(CSFbio)[(!is.na(CSFbio$Ptau_ASSAY_Zscore)) & 
                              (!is.na(CSFbio$AB_Zscore)) &
                              (!is.na(CSFbio$Ttau_ASSAY_Zscore))]

# Add CSF info
Y_train <- Y_train[intersect(samples, rownames(Y_train)),]
Y_test <- Y_test[intersect(samples, rownames(Y_test)),]
X_train <- cbind.data.frame(X_train[rownames(Y_train),], CSFbio[rownames(Y_train),])
X_test <- cbind.data.frame(X_test[rownames(Y_test),], CSFbio[rownames(Y_test),])


###############################################################################

# Get SHAP values

###############################################################################

# Load MCI model
load("Model/EMIF_Models/CSF/Fit_EMIF_MCI_sPLS_CSFbio.RData")

# Make explained object
explainer <- DALEX::explain(fit, data = X_train, y = NULL)

# Calculate SHAP values
for (i in 1:nrow(X_test)){
  shap <- predict_parts(explainer = explainer, 
                        new_observation = X_test[i,], 
                        type = "shap",
                        B = 50)
  
  # get mean contribution (i.e., SHAP value)
  shap_fil <- shap %>%
    group_by(variable_name) %>%
    reframe(
      Contribution = mean(contribution))
  
  if (i == 1){
    output <- shap_fil
  } else{
    output <- inner_join(output, shap_fil, by = c("variable_name" = "variable_name"))
  }
  
}
output <- as.data.frame(output)
rownames(output) <- output$variable_name
output <- output[,-1]
colnames(output) <- rownames(X_test)

# Save SHAP values 
save(output, file = "Models/ModelInterpretation/SHAP/output_shap_sPLS_CSFbio.RData")


###############################################################################

# Compare SHAP values between different models

###############################################################################

#==============================================================================#
# MCI models from EMIF
#==============================================================================#

# Models
models <- c("EN", "sPLS", "RF")

plotDF <- NULL
for (i in 1:length(models)){
  
  # Load SHAP values
  load(paste0("Models/ModelInterpretation/SHAP/output_shap_", models[i], "_CSFbio.RData"))
  
  # Scale SHAP values per individual
  output_scaled <- t(t(output)/colSums(abs(output)))
  
  # Get average SHAP values
  temp <- data.frame(Variable = rownames(output_scaled),
                     value = rowMeans(abs(output_scaled)))
  
  if (i > 1){
    plotDF <- inner_join(plotDF, temp, by = c("Variable" = "Variable"))
  } else{
    plotDF <- temp
  }
}
plotDF <- plotDF[plotDF$Variable != "ChrAge",]

# Set row names
rownames(plotDF) <- c("Amyloid-beta", "Age", "Alcohol Intake", "BMI", "Depression", "Type II Diabetes",
                      "Dietary Intake", "Education", "HDL Chol.", "Heart Disease",
                      "Physical Act.", "P-tau", "Sex", "Smoking", "Syst. Blood Pressure", "Total Chol.",
                       "T-tau")
plotDF <- plotDF[,-1]

# Set column names
colnames(plotDF) <- c("ElasticNet", "sPLS", "Random Forest")

# Scale SHAP values
for (i in 1:ncol(plotDF)){
  plotDF[,i] <- plotDF[,i]/sum(plotDF[,i])
}


# Prepare data for plotting
data <- as.data.frame(t(plotDF))
data <- rbind(rep(0.3,5) , rep(0,5) , data)
data <- data[,c("Amyloid-beta", "Age", "Alcohol Intake", "BMI", "Depression", "Type II Diabetes",
                "Dietary Intake", "Education", "HDL Chol.", "Heart Disease",
                "Physical Act.", "P-tau", "Sex", "Smoking", "Syst. Blood Pressure", "Total Chol.",
                "T-tau")]

# Color vector
colors_border=rev(c("#EF3B2C","#CB181D", "#99000D") )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )

# plot with default options:
jpeg("Models/ModelInterpretation/SHAP/shaply_emif_CSF.jpg", width = 8000, height = 7000, quality = 100)
radarchart( data  , axistype=0 , 
            #custom polygon
            pcol=colors_border, plwd=100 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,0.1,0.05), cglwd=8,
            #custom labels
            vlcex=20 
)
dev.off()

