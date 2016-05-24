#Ignore this. This file contains some routines I use to get the data into the package
TempDataFetch <- function() {

  loadMyData <- function(filename, doi) {
    names <- read.csv("S:/OMOP/Null/SCCSLikeTata/DrugNames.csv")
    data <- read.csv(filename)
    data$LOGRR <- data$Beta
    data$LOGLB95RR <- log(data$CI95Down)
    data$LOGUB95RR <- log(data$CI95Up)
    data$DRUG_CONCEPT_ID <- data$ATC
    data$DRUG_CONCEPT_NAME <- data$ATC
    data$CONDITION_CONCEPT_ID <- data$EventType
    data$GROUND_TRUTH <- rep(0, nrow(data))
    data$GROUND_TRUTH[data$DRUG_CONCEPT_ID == doi] <- 1
    data$SE <- (data$LOGRR - data$LOGLB95R)/qnorm(0.975)
    data$Z <- data$LOGRR/data$SE
    data$P <- 2 * pmin(pnorm(data$Z), 1 - pnorm(data$Z))  # 2-sided p-value
    data$SIGNIFICANT <- data$LOGLB95RR > 0 | data$LOGUB95RR < 0
    # data <- data[!is.na(data$Z) & data$Z != 0 & !is.infinite(data$Z),]
    data <- merge(data, names, by.x = "DRUG_CONCEPT_ID", by.y = "drug_concept_id")
    data$DRUG_CONCEPT_NAME <- data$drug_concept_name
    data
  }

  method <- "SCCS"
  doi <- 739138
  filename <- "S:/OMOP/Null/SCCSLikeTata/SCCSresults.csv"
  HOI <- 500001003
  data <- loadMyData(filename, doi)
  data <- data[!is.na(data$Z) & data$Z != 0 & !is.infinite(data$Z), ]
  data <- data[, c("DRUG_CONCEPT_NAME", "GROUND_TRUTH", "LOGRR", "SE")]
  sccs <- data
  colnames(sccs) <- c("drugName", "groundTruth", "logRr", "seLogRr")
  save(sccs, file = "data/sccs.rda")

  # logRr <- sccs$logRr seLogRr <- sccs$seLogRr

  method <- "CC"
  doi <- 739138
  filename <- "S:/OMOP/Null/SCCSLikeTata/CCresults.csv"
  HOI <- 500001003
  data <- loadMyData(filename, doi)
  data <- data[!is.na(data$Z) & data$Z != 0 & !is.infinite(data$Z), ]
  data <- data[, c("DRUG_CONCEPT_NAME", "GROUND_TRUTH", "LOGRR", "SE")]
  caseControl <- data
  colnames(caseControl) <- c("drugName", "groundTruth", "logRr", "seLogRr")
  save(caseControl, file = "data/caseControl.rda")


  loadData <- function(method, analysisSourceID, HOI) {
    source <- paste("S:/OMOP/Null/Data_", method, ".csv", sep = "")
    data <- read.csv(source)
    data$ANALYSIS_ID <- as.factor(data$ANALYSIS_ID)
    data$GROUND_TRUTH <- as.factor(data$GROUND_TRUTH)
    data$ANALYSIS_SOURCE_ID <- paste(data$ANALYSIS_ID, "-", data$SOURCE_ID)
    data <- data[data$ANALYSIS_SOURCE_ID == analysisSourceID, ]
    data <- data[data$CONDITION_CONCEPT_ID == HOI, ]
    data$SE <- (data$LOGRR - data$LOGLB95R)/qnorm(0.975)
    data$Z <- data$LOGRR/data$SE
    data$P <- 2 * pmin(pnorm(data$Z), 1 - pnorm(data$Z))  # 2-sided p-value
    data <- data[!is.na(data$Z), ]
    data$LOGUB95RR <- data$LOGRR + data$LOGRR - data$LOGLB95RR
    data$SIGNIFICANT <- data$LOGLB95RR > 0 | data$LOGUB95RR < 0
    data$DRUG_CONCEPT_NAME <- as.character(data$DRUG_CONCEPT_NAME)
    data
  }
  method <- "CM"
  doi <- 1782521
  analysisSourceID <- "21000211 - 28"  #MDCR
  HOI <- 500000301

  data <- loadData(method, analysisSourceID, HOI)
  data <- data[data$GROUND_TRUTH == 0 | data$DRUG_CONCEPT_ID == doi, ]
  # require(Hmisc)
  #data$DRUG_CONCEPT_NAME <- Hmisc::capitalize(data$DRUG_CONCEPT_NAME)
  #data$DRUG_CONCEPT_NAME[data$DRUG_CONCEPT_NAME == "Sodium Phosphate, Monobasic"] <- "Sodium Phosphate"
  #data <- data[!is.na(data$Z) & data$Z != 0 & !is.infinite(data$Z), ]
  #data <- data[, c("DRUG_CONCEPT_NAME", "GROUND_TRUTH", "LOGRR", "SE")]
  #cohortMethod <- data
  #colnames(cohortMethod) <- c("drugName", "groundTruth", "logRr", "seLogRr")
  #save(cohortMethod, file = "data/cohortMethod.rda")

}
