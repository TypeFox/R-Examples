######################
# tauToPearson

# Function: tauToPearson
# Description: If the two variables are bivariate normal distributed, this formula maps kendalls tau to bravais pearson correlation coefficient
# Input:
# Tau: Numeric vector in the interval [-1, 1]
# Output: Vector of Bravais Pearson correlation coefficient

tauToPearson <- function (Tau) {
  sin(Tau*pi/2)
}

########################
# simCompRiskDisc

#######################
# Description
# Simulates responses and covariates of discrete competing risk models
# Responses and covariates are modelled separately by gaussian copulas with given kendalls tau

#######################
# Inputs

# Input rank correlation structures
# responseCorr: Regulates correlation between event survival times. Numeric Matrix with kendalls tau version b rank correlations 
  # Each cell is restricted to be between -1 and 1
  # Diagonal elements are always 1
# covariateCorr: Regulates correlation between event covariates. Numeric Matrix with kendalls tau version b rank correlations (each cell is restricted to be between -1 and 1)
# Each cell is restricted to be between -1 and 1
# Diagonal elements are always 1
# Uses singular, value decomposition for invertation of covariance matrices

# sampleSize: Integer scalar specifying the number of rows of the data frame

# Input for covariates:
# covariateSeed: Integer scalar, specifying seed of covariate random number generation
# covariateQuantFunc: Character vector, giving the marginal quantile functions of all covariates
# covariateArgs: List of Arguments for each quantile function of the marginal distributions
# Each list element should be a named numeric vector (names giving the arguments)

# intercept: Logical vector, giving if intercept is given in true coefficients for each Event (Default all TRUE)
# trueCoef: List of numeric vectors containing the \beta of true coefficients, 
# e. g. linear predictor \eta = X %*% \beta; with X as simulated design matrix
# Each list element has all \beta coefficients of one event

# Input for responses:
# responseSeed: Integer scalar, specifying seed of response random number generation
# responseQuantFunc: Character vector, giving the marginal quantile functions of all survival events
# responseFixArgs: List of Arguments for each quantile function of the marginal distributions
# Each list element should be a named numeric vector
# responseLinPredArgs: List of lists, specifying the relationship of linear predictor and parameters of the marginal distributions
# Each list element is a list of all functional relationships between linear predictor and parameters of one marginal distribution
# Each list element is a function giving the functional relationship between linear predictor and one parameter

# Input for censoring:
# censorRN: Integer scalar, specifying seed of censoring random number generation
# censorArgs: Named numeric vector, giving all arguments of the marginal censoring distribution
# censorSeed: Integer scalar, specifying seed of censoring random number generation   

####################
# Output

# List with following components
# Data: Simulated data as data.frame
# responseCorr: Original input rank correlation structure of responses
# covariateCorr: Original input rank correlation structure of covariates
# covariateSeed: Covariate seed
# ... (all arguments specified in Input are saved other the corresponding names)

###################
# Main Code

simCompRisk <- function (responseCorr, covariateCorr, sampleSize, covariateSeed=NULL, covariateQuantFunc, covariateArgs, 
                           intercept=TRUE, trueCoef, responseSeed, responseQuantFunc, responseFixArgs, responseLinPredArgs, censorRN, censorArgs, censorSeed) {
  # Count of covariates
  CountCovariates <- nrow(covariateCorr)
  CountResponses <- nrow(responseCorr)
  RN_generator <- NULL
  
  # Transform to Pearson Correlation
  responseCorrOrg <- responseCorr
  covariateCorrOrg <- covariateCorr
  responseCorr <- apply(responseCorr, c(1,2), tauToPearson)
  covariateCorr <- apply(covariateCorr, c(1,2), tauToPearson)
  
  # Simulate variables with covariateCorr (based on Kendalls tau)
  set.seed(covariateSeed)
  RN_copula <- rmvnorm(n=sampleSize, mean = rep(0, CountCovariates), sigma = covariateCorr, method="svd")
  RN_covariate <- matrix(NA, nrow=sampleSize, ncol=CountCovariates, dimnames=list(NULL, paste("Var", 1:CountCovariates, sep="")))
  for(i in 1:CountCovariates) {
    RN_probs <- pnorm(q=RN_copula [, i], mean = 0, sd = 1)
    CodeString <- paste("RN_generator <- Curry(eval(parse(text=covariateQuantFunc [i])), ", paste(names(covariateArgs [[i]]), covariateArgs [[i]], sep="=", collapse=", "), ")", sep="")
    eval(parse(text=CodeString))
    RN_covariate [, i] <- RN_generator (RN_probs)
  }
  
  # Construct linear predictor
  LinPred <- vector("list", CountResponses)
  for(i in 1:CountResponses) {
    if(intercept [i]) {
      LinPred [[i]] <- cbind(1, RN_covariate) %*% trueCoef [[i]]
    }
    else {
      LinPred [[i]] <- RN_covariate %*% trueCoef [[i]]
    }
  }
  
  # Simulate survival times with dependence structure (based on Kendalls tau)
  # Caution: Linear predictor must be included in the parameters of the response quantile functions!
  set.seed(responseSeed)
  RN_copula <- rmvnorm(n=sampleSize, mean = rep(0, CountResponses), sigma = responseCorr)
  RN_response <- matrix(NA, nrow=sampleSize, ncol=CountResponses, dimnames=list(NULL, paste("Resp", 1:CountResponses, sep="")))
  for(i in 1:CountResponses) {
    
    # Generate dependent multivariate gaussian observations
    RN_probs <- pnorm(q=RN_copula [, i], mean = 0, sd = 1)
    TempStringFix <- paste(names(responseFixArgs [[i]]), responseFixArgs [[i]], sep="=", collapse=", ")
    
    for(j in 1:sampleSize) {
      
      # Generate survival times of specified marginal distributions with parameters dependent on covariates
      TransformLinPred <- sapply(1:length(responseLinPredArgs [[i]]), function (x) responseLinPredArgs [[i]] [[x]] (LinPred [[i]] [j]))
      TempStringLinPred <- paste(names(responseLinPredArgs [[i]]), TransformLinPred, sep="=", collapse=", ")
      CodeString <- paste("RN_generator <- Curry(eval(parse(text=responseQuantFunc [i])), ", TempStringFix, ", ", TempStringLinPred, ")", sep="")
      eval(parse(text=CodeString))
      RN_response [j, i] <- RN_generator (RN_probs [j])
      
    }
  }
  
  # Simulate censoring times independent of survival times
  CodeString <- paste("RN_generator <- Curry(eval(parse(text=censorRN)), ", paste(names(censorArgs), censorArgs, sep="=", collapse=", "), ")", sep="")
  eval(parse(text=CodeString))
  set.seed(censorSeed)
  RN_censor <- RN_generator (sampleSize)
  
  # Calculate observed times and event indicators
  RN_response <- cbind(RN_censor, RN_response)
  EventInd <- sapply(1:sampleSize, function (x) which.min(RN_response [x, ]))
  ObsTime <- sapply(1:sampleSize, function (x) RN_response [x, EventInd [x]])
  EventInd <- model.matrix(~I(factor(EventInd, levels=1:(CountResponses+1)))-1)
  dimnames(EventInd) [[2]] <- c("Censor", paste("Event", 1:CountResponses, sep=""))
  
  # Output and store input parameters
  temp <- cbind(Time=ObsTime, Event=EventInd, RN_covariate)
  rownames (temp) <- 1:sampleSize
  DataComplete <- as.data.frame(temp)
  Output <- list(Data=DataComplete, responseCorr=responseCorrOrg, covariateCorr=covariateCorrOrg, sampleSize=sampleSize, covariateSeed=covariateSeed, covariateQuantFunc=covariateQuantFunc, 
                 covariateArgs=covariateArgs, intercept=intercept, trueCoef=trueCoef, responseSeed=responseSeed, responseQuantFunc=responseQuantFunc, responseFixArgs=responseFixArgs, 
                 responseLinPredArgs=responseLinPredArgs, censorRN=censorRN, censorArgs=censorArgs, censorSeed=censorSeed)
  class(Output) <- "discSurvSimCompRisk"
  return(Output)
}

print.discSurvSimCompRisk <- function (x, ...) {
  print(x$Data)
}