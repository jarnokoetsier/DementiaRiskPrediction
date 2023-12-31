# ============================================================================ #
# File: PlotMPSs_EMIF.R
# Author: Jarno Koetsier
# Date: August 6, 2023
# Description: Plot distribution of MPSs among the diagnostic groups in the 
#              EMIF-AD cohort.
# ============================================================================ #

# Load packages
library(prospectr)
library(tidyverse)
library(caret)
library(patchwork)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF-AD/Data/metaData_fil.RData")                # Meta data
load("EMIF-AD/Data/predictedScore_factors_fil.RData")  # MRSs

# Check data
all(metaData_fil$X == rownames(predictedScore_factors_fil))
table(metaData_fil$CTR_Convert)
table(metaData_fil$Diagnosis)

# Format predicted scores
testDF <- predictedScore_factors_fil
colnames(testDF) <- c("Syst. Blood Pressure", "Total Chol.", "Low Education",
                      "Physical Inact.", "Unhealthy Diet", "Depression",
                      "Type II Diabetes", "Heart Disease", "Sex (male)", "Age",
                      "Alcohol Intake",  "BMI","HDL Chol.","Smoking")

plotDF <- gather(testDF)
plotDF$SampleID <- rep(rownames(testDF, ncol(testDF)))

# Add diagnosis
plotDF <- inner_join(plotDF, metaData_fil[,c("X", "Diagnosis")],
                     by = c("SampleID" = "X"))
plotDF$Diagnosis[plotDF$Diagnosis == "NL"] <- "Control"
plotDF$Diagnosis <- factor(plotDF$Diagnosis,
                           levels = c("Control", "SCI","MCI", "AD")) 

# For each risk factor (MPS), make boxplot
for (i in 1:ncol(testDF)){
  
  p <- ggplot(plotDF[plotDF$key == colnames(testDF)[i],]) +
    geom_boxplot(aes(x = Diagnosis, y = value, fill = Diagnosis)) +
    xlab(NULL) +
    ylab("Methylation Profile Score") +
    ggtitle(colnames(testDF)[i]) +
    scale_fill_manual(values = c("#FCBBA1","#FB6A4A","#CB181D","#99000D")) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5,
                                    face = "bold",
                                    size = 16),
          plot.subtitle = element_text(hjust = 0.5,
                                       size = 10,
                                       face = "italic"))
  
  if (i == 1){
    all <- p
  } else{
    all <- all + p
  }
  
}
finalPlot <- all +  
  plot_layout(ncol = 4, nrow = 4)

# Save plot
ggsave(finalPlot, file = "EMIF-AD/ModelPerformance/FactorVsDiagnosis.jpg", width = 12, height = 8)

# ANOVA (PARAMETRIC)
test <- aov(value ~ Diagnosis, data = plotDF[plotDF$key== "Depression",])
summary(test)
TukeyHSD(test)

# KRUSKAL-WALLIS(NON-PARAMETRIC)
data1 <- plotDF[plotDF$key== "SexMale",]
test <- kruskal.test(value ~ Diagnosis, data = data1)
test
pairwise.wilcox.test(data1$value, data1$Diagnosis,
                     p.adjust.method = "BH")

