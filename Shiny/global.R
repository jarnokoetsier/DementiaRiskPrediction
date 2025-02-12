options(timeout=3600)

# Load packages
library(caret)
library(glmnet)
library(ranger)
library(data.table)
library(DT)
library(shinyWidgets)
library(shinycssloaders)

# Install files that are not in Shiny directory
all_files <- list.files()
files <- c("finalModels.RData", 
           "all_cpgs.RData", 
           "example.csv",
           "Fit_EMIF_MCI_RF.RData")

for (f in files){
  if (!(f %in% all_files)){
    download.file(url = paste0("https://zenodo.org/record/8306113/files/",f), 
                  f, mode = "wb")
  }
}

if (!("Predictors_Shiny_by_Groups.csv" %in% all_files)){
  download.file(url = "https://zenodo.org/record/4646300/files/Predictors_Shiny_by_Groups.csv", 
                "Predictors_Shiny_by_Groups.csv", mode = "wb")
}




# Load EXTEND models
load("finalModels.RData")

# Load model's CpGs
load("all_cpgs.RData")

# Load Marioni models
cpgs <- read.csv("Predictors_Shiny_by_Groups.csv", header = T) 

# Load epi-MCI model
load("Fit_EMIF_MCI_RF.RData")


#******************************************************************************#
# Convert NAs to Mean Value for all individuals across each probe 
#******************************************************************************#
na_to_mean <-function(methyl){
  methyl[is.na(methyl)]<-mean(methyl,na.rm=T)
  return(methyl)
}

#******************************************************************************#
# Conversion to Beta Values
#******************************************************************************#
m_to_beta <- function (val){
  beta <- 2^val/(2^val + 1)
  return(beta)
}

#******************************************************************************#
# Conversion to M Values
#******************************************************************************#
beta_to_m <- function (val){
  m <- log2(val/(1-val))
  return(m)
}

#******************************************************************************#
# Run EXTEND models
#******************************************************************************#
extendModel <- function(data, finalModels, all_cpgs){
  
  # Check if Data needs to be Transposed
  if(ncol(data) < nrow(data)){
    message("It seems that individuals are columns - data will be transposed!")
    data <- t(data) 
  }
  
  # Copy data
  coef <- data
  
  # Conversion to M Values if Beta Values are likely present  
  coef<-if((range(coef,na.rm=T)<= 1)[[2]] == "TRUE") { message("Suspect that Beta Values are present. Converting to M Values");beta_to_m(coef) } else { message("Suspect that M Values are present");coef}
  
  # Identify CpGs missing from input dataframe, include them and provide values 
  # as mean methylation value at that site
  if (!all((all_cpgs$id %in% colnames(coef)))){
    missing_cpgs <- all_cpgs$id[!(all_cpgs$id %in% colnames(coef))]
    mat <- matrix(ncol=length(unique(missing_cpgs)), nrow = nrow(coef))
    colnames(mat) <- unique(missing_cpgs)
    rownames(mat) <- rownames(coef) 
    
    for (c in 1:ncol(mat)){
      mat[,c] <- rep(all_cpgs[colnames(mat)[c],"meanValue"])
    }
    coef <- cbind(coef,mat)
    
  }
  
  # Convert NAs to Mean Value for all individuals across each probe 
  coef <- apply(coef,2,function(x) na_to_mean(x))
  coef <- data.frame(coef)
  
  # Predict factors
  factors <- names(finalModels)
  predictedScore_factors <- matrix(NA, nrow = nrow(coef), ncol = length(factors))
  colnames(predictedScore_factors) <- factors
  rownames(predictedScore_factors) <- rownames(coef)
  for (f in 1:length(factors)){
    f1 <- factors[f]
    model <- finalModels[[f1]]
    
    # Get features for model fitting
    features <- colnames(model$trainingData)[-ncol(model$trainingData)]
    
    # Make predictions
    if (f <= 2){
      predictedScore_factors[,f] <- predict(model, coef[,features])
    }
    if (f > 2){
      prob <- predict(model, coef[,features], type = "prob")$Yes
      predictedScore_factors[,f] <- log(prob/(1-prob))
    }
    
  }
  predictedScore_factors <- as.data.frame(predictedScore_factors)
  return(predictedScore_factors)
}

#******************************************************************************#
# Run Marioni models
#******************************************************************************#
marioniModel <- function(data, cpgs){
  
  # Check if Data needs to be Transposed
  if(ncol(data) > nrow(data)){
    message("It seems that individuals are rows - data will be transposed!")
    data <- t(data) 
  }
  
  # Subset CpG sites to those present on list for predictors 
  coef <- data[intersect(rownames(data), cpgs$CpG_Site),]
  
  # Identify CpGs missing from input dataframe, include them and provide values 
  # as mean methylation value at that site
  
  coef <- if(nrow(coef) == 3629) { message("All sites present"); coef } else if(nrow(coef)==0){ 
    message("There Are No Necessary CpGs in The dataset - All Individuals Would Have Same Values For Predictors. Analysis Is Not Informative!")
  } else { 
    missing_cpgs = cpgs[-which(cpgs$CpG_Site %in% rownames(coef)),c("CpG_Site","Mean_Beta_Value")]
    message(paste(length(unique(missing_cpgs$CpG_Site)), "unique sites are missing - add to dataset with mean Beta Value from Training Sample", sep = " "))
    mat = matrix(nrow=length(unique(missing_cpgs$CpG_Site)),ncol = ncol(coef))
    row.names(mat) <- unique(missing_cpgs$CpG_Site)
    colnames(mat) <- colnames(coef) 
    mat[is.na(mat)] <- 1
    missing_cpgs1 <- if(length(which(duplicated(missing_cpgs$CpG_Site))) > 1) { 
      missing_cpgs[-which(duplicated(missing_cpgs$CpG_Site)),]
    } else {missing_cpgs
    }  
    ids = unique(row.names(mat))
    missing_cpgs1 = missing_cpgs1[match(ids,missing_cpgs1$CpG_Site),]
    mat=mat*missing_cpgs1$Mean_Beta_Value
    coef=rbind(coef,mat)
  } 
  

  # Convert NAs to Mean Value for all individuals across each probe 
  coef <- t(apply(coef,1,function(x) na_to_mean(x)))
  
  # Conversion to Beta Values if M Values are likely present  
  coef<-if((range(coef,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(coef) } else { message("Suspect that Beta Values are present");coef}
  
  # Calculating the predictors 
  loop = unique(cpgs$Predictor)
  out <- data.frame()
  for(i in loop){ 
    tmp=coef[intersect(row.names(coef),cpgs[cpgs$Predictor %in% i,"CpG_Site"]),]
    tmp_coef = cpgs[cpgs$Predictor %in% i, ]
    if(nrow(tmp_coef) > 1) { 
      tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_Site),]
      out[colnames(coef),i]=colSums(tmp_coef$Coefficient*tmp)
    } else {
      tmp2 = as.matrix(tmp)*tmp_coef$Coefficient 
      out[colnames(coef),i] = colSums(tmp2)
    }
  } 
  out$'Epigenetic Age (Zhang)' <- out$'Epigenetic Age (Zhang)' + 65.79295
  out$ID <- row.names(out) 
  out <- out[,c(ncol(out),1:(ncol(out)-1))] 
  colnames(out) <- c("ID", "EpiAge", "Alcohol", "BMI", "BodyFat", "HDL", "Smoking", "WHR")
  return(out)
}