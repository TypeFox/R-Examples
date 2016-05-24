#########################
# brierScore

# Description:
# Computes the Brier Score of person i, based on generalized, linear models given Survival and Censoring functions

# Steps
# 0. Extract training and test data
# 1. Convert training data to long format
# 2. Convert test data to long format
# 3. Convert response in training data to censoring variable
# 4. Fit survival model on training data in long format
# 5. Fit censoring model on training data in long format
# 6. Predict survival times on test data
# 7. Predict censoring times on test data
# 8. Estimate survival functions of survival times for each person (aggregate)
# 9. Estimate survival functions of censoring times for each person (aggregate)
# 10. Calculate brier score

# Input
# dataSet: Original data in short format. Must be of type data.frame
# trainIndices: List of integer Indices, which give the rows of *dataSet* as training data in cross validation
# survModelFormula: Formula of the survival process. First argument must be the univariate time (numeric) variable
# censModelFormula: Formula of the censoring process. First argument must be the univariate event (binary) variable 
# linkFunc: Gives the type of link used in the glm function. May be customized
# idColumn: Gives the column name of the identification numbers as character scalar. 
# Default NULL means, that each row equals one person (no repeated measurements)

# Output
# Numeric vector giving the brier scores of the persons available in all test subsets

brierScore <- function (dataSet, trainIndices, survModelFormula, linkFunc="logit", censColumn, idColumn=NULL) {
  
  # Input Checks
  if(!is.data.frame(dataSet)) {stop("Argument *dataSet* is not in the correct format! Please specify as data.frame object.")}
  if(!is.list(trainIndices)) {stop("Argument *trainIndices* is not in the correct format! Please specify a list.")}
  InputCheck1 <- all(sapply(1:length(trainIndices), function (x) is.integer(trainIndices [[x]])))
  if(!InputCheck1) {stop("Sublists of *trainIndices* are not all integer values! Please specify a list of integer Indices.")}
  if(!("formula" %in% class(survModelFormula))) {stop("*survModelFormula* is not of class formula! Please specify a valid formula, e. g. y ~ x + z.")}
  if(!any(names(dataSet)==censColumn)) {stop("Argument *censColumn* is not available in *dataSet*! Please specify the correct column name of the event indicator.")}
  
  # Help function
  B <- function(k) {
    probs <- estMargProb (LambdaSplit [[k]] [, "Lambda" ])
    probs <- probs [-length(probs)]
    
    if(length(probs [-length(probs)])!=0) {
      brierVec <- as.numeric(tail(LambdaSplit [[k]] [, censColumn], 1) * (1 - tail(probs, 1))^2 + sum (probs [-length(probs)]))
    }
    else {
      brierVec <- as.numeric(tail(LambdaSplit [[k]] [, censColumn], 1) * (1 - tail(probs, 1))^2)
    }
    return(brierVec)
  }
  
  RET <- vector("list", length(trainIndices))
  for(i in 1:length(trainIndices)) {
    
    # 0. Extract training, test data and responses
    TrainSet <- dataSet [trainIndices [[i]], ]
    if(length(trainIndices)!=1) {
      TestSet <- dataSet [-trainIndices [[i]], ]
    }
    else {
      TestSet <- TrainSet
    }
    
    # 1. Convert training data to long format
    if(!is.null(idColumn)) {
      TrainLong <- dataLongTimeDep (dataSet=TrainSet, timeColumn=as.character(survModelFormula) [2], censColumn=censColumn, idColumn=idColumn)
    }
    else {
      TrainLong <- dataLong (dataSet=TrainSet, timeColumn=as.character(survModelFormula) [2], censColumn=censColumn)
    }
    
    # 2. Convert response in training data to censoring variable
    TrainLong <- dataCensoring (dataSetLong=TrainLong, respColumn="y", idColumn="obj")
    
    # 3. Convert test data to long format
    if(!is.null(idColumn)) {
      TestLong <- dataLongTimeDep (dataSet=TestSet, timeColumn=as.character(survModelFormula) [2], censColumn=censColumn, idColumn=idColumn)
    }
    else {
      TestLong <- dataLong (dataSet=TestSet, timeColumn=as.character(survModelFormula) [2], censColumn=censColumn)
    }

    SurvnewFormula <- update(survModelFormula, y ~ timeInt + .)
    SurvFit <- glm (formula=SurvnewFormula, data=TrainLong, family=binomial(link=linkFunc), control=glm.control(maxit=2500))

    # 6. Estimate survival curves on test data
    Check <- "error" %in% class(tryCatch(predict(SurvFit, TestLong, type="response"), error= function (e) e))
    if(Check) {
      # Which columns are factors in test data?
      IndexFactor <- which(sapply(1:dim(TestLong)[2], function (x) is.factor(TestLong [, x]))==TRUE)
      # What are the levels of these factors?
      TestLevelsFactor <- sapply(IndexFactor, function (x) levels(TestLong [, x]))
      # First column does not count (censoring process)
      # What are the levels of the corresponding factors in the training data?
      TrainLevelsFactor <- sapply(IndexFactor, function (x) levels(TrainLong [, x+1]))
      # Which levels of the test data exist in the training data factors?
      InLevelsFactor <- lapply(1:length(TestLevelsFactor), function (x) which((TestLevelsFactor [[x]] %in% TrainLevelsFactor [[x]])==FALSE))
      ExcludeRows <- lapply (1:length(IndexFactor), function (j) which(TestLong [, IndexFactor [j]] %in% TestLevelsFactor [[j]] [InLevelsFactor [[j]]]))
      ExcludeRows <- do.call(c, ExcludeRows)
      TestLong <- TestLong [-ExcludeRows, ]
    }
    Lambda <- predict(SurvFit, TestLong, type="response")
    LambdaSplit <- split(cbind(Lambda=Lambda, TestLong), TestLong$obj)
    RET [[i]] <- sapply(1:length(LambdaSplit), B)
  }

  # Output
  RET <- do.call(cbind, RET)
  RET <- rowMeans(RET)
  return(RET)
}

#############################
# tprUno
#############################

##############
# Description
# Estimates the predictive true positive rate (tpr) based on cross validation and generalized, linear models

#######
# Input
# timepoint: Discrete time interval given that the false positive rate is evaluated (integer scalar)
# dataSet: Original data. Should be in format data.frame()
# trainIndices: List of Indices from original data used for training (list of integer vectors). 
# The length of the list is equal to the number of cross valdiation samples
# survModelFormula: Formula of the survival model
# censModelFormula: Formula of the censoring model. Normally this is done without covariates
# linkFunc: Link function of the generalized, linear model see glm
# idColumn: Name of the column with identification numbers of persons. 
# Default NULL means, that each row equals one person (no repeated measurements).

# Output
# data.frame with columns
# cutoff: Cut off values of the linear predictor (numeric vector)
# fpr: False positive rate (numeric vector)

tprUno <- function(timepoint, dataSet, trainIndices, survModelFormula, censModelFormula, linkFunc="logit", idColumn=NULL, timeAsFactor=TRUE) {
  
  # Input Checks
  if(length(timepoint)!=1 || !(timepoint==floor(timepoint))) {stop("Argument *timepoint* is not in the correct format! Please specify as integer scalar value.")}
  if(!is.data.frame(dataSet)) {stop("Argument *dataSet* is not in the correct format! Please specify as data.frame object.")}
  if(!is.list(trainIndices)) {stop("Argument *trainIndices* is not in the correct format! Please specify a list.")}
  InputCheck1 <- all(sapply(1:length(trainIndices), function (x) is.integer(trainIndices [[x]])))
  if(!InputCheck1) {stop("Sublists of *trainIndices* are not all integer values! Please specify a list of integer Indices.")}
  if(length(trainIndices)!=1) {
    InputCheck2 <- all(sort(as.numeric(do.call(c, lapply(trainIndices, function (x) setdiff(1:dim(dataSet) [1], x)))))==(1:dim(dataSet) [1]))
  }
  else {
    InputCheck2 <- all(trainIndices [[1]]==(1:dim(dataSet) [1]))
  }
  if(!InputCheck2) {stop("Argument *trainIndices* does not contain cross validation samples! Please ensure that the union of all test indices equals the indices of the complete data set.")}
  if(!("formula" %in% class(censModelFormula))) {stop("*censModelFormula* is not of class formula! Please specify a valid formula, e. g. yCens ~ 1")}
  if(!("formula" %in% class(survModelFormula))) {stop("*survModelFormula* is not of class formula! Please specify a valid formula, e. g. y ~ x")}
  if(!(any(names(dataSet)==idColumn) | is.null(idColumn))) {stop("Argument *idColumn* is not available in *dataSet*! Please specify the correct column name of the identification number.")}
  
  # Help function
  sens <- function(k) {
    
    sensNum <- sum((marker > k) * (newTime == timepoint) * newEvent / GT, na.rm = TRUE)
    sensDenom <- sum((newTime == timepoint) * newEvent / GT, na.rm = TRUE)
    
    if (sensDenom > 0)
      return(sensNum / sensDenom) else
        return(0)
  }
  
  # Loop across all training data sets
  markerList <- vector("list", length(trainIndices))
  ExcludeRowsCensList <- vector("list", length(trainIndices))
  ExcludeRowsDataSetList <- vector("list", length(trainIndices))
  oneMinuslambdaList <- vector("list", length(trainIndices))
  
  # Convert full sample to long format
  if(!is.null(idColumn)) {
    TrainLongFull <- dataLongTimeDep (dataSet=dataSet, timeColumn=as.character(survModelFormula) [2], 
                                      censColumn=as.character(censModelFormula) [2], 
                                      idColumn=idColumn, timeAsFactor=timeAsFactor)
  }
  else {
    TrainLongFull <- dataLong (dataSet=dataSet, timeColumn=as.character(survModelFormula) [2], 
                               censColumn=as.character(censModelFormula) [2], timeAsFactor=timeAsFactor)
  }
  
  for(i in 1:length(trainIndices)) {
    
    # 0. Extract training, test data and responses
    TrainSet <- dataSet [trainIndices [[i]], ]
    if(length(trainIndices)!=1) {
      TestSet <- dataSet [-trainIndices [[i]], ]
    }
    else {
      TestSet <- TrainSet
    }

    # 1. Convert training data to long format
    if(!is.null(idColumn)) {
      TrainLong <- dataLongTimeDep (dataSet=TrainSet, timeColumn=as.character(survModelFormula) [2], 
                                    censColumn=as.character(censModelFormula) [2], idColumn=idColumn,
                                    timeAsFactor=timeAsFactor)
    }
    else {
      TrainLong <- dataLong (dataSet=TrainSet, timeColumn=as.character(survModelFormula) [2], 
                             censColumn=as.character(censModelFormula) [2], timeAsFactor=timeAsFactor)
    }
    
    # 2. Convert response in training data to censoring variable
    TrainLong <- dataCensoring (dataSetLong=TrainLong, respColumn="y", idColumn="obj")
    
    # 3. Convert test data to long format
    if(!is.null(idColumn)) {
      TestLong <- dataLongTimeDep (dataSet=TestSet, timeColumn=as.character(survModelFormula) [2], 
                                   censColumn=as.character(censModelFormula) [2], idColumn=idColumn,
                                   timeAsFactor=timeAsFactor)
    }
    else {
      TestLong <- dataLong (dataSet=TestSet, timeColumn=as.character(survModelFormula) [2], 
                            censColumn=as.character(censModelFormula) [2], timeAsFactor=timeAsFactor)
    }

    # 4. Fit censoring model on training data in long format
    CensnewFormula <- update(censModelFormula, yCens ~ timeInt + .)
    CensFit <- glm (formula=CensnewFormula, data=TrainLong, family=binomial(link=linkFunc), control=glm.control(maxit=2500))
    
    # 5. Estimate censoring curves on test data
    # Exclude cases with new factor levels in test data in long format
    Check <- "error" %in% class(tryCatch(predict(CensFit, TestLong, type="response"), error= function (e) e))
    if(Check) {
      
      # Which columns are factors in test data?
      IndexFactor <- which(sapply(1:dim(TestLong)[2], function (x) is.factor(TestLong [, x]))==TRUE)
      
      # What are the levels of these factors?
      TestLevelsFactor <- sapply(IndexFactor, function (x) levels(TestLong [, x]))
      
      # First column does not count (censoring process)
      # What are the levels of the corresponding factors in the training data?
      TrainLevelsFactor <- sapply(IndexFactor, function (x) levels(TrainLong [, x+1]))
      
      # Which levels of the test data exists in the training data factors?
      InLevelsFactor <- lapply(1:length(TestLevelsFactor), function (x) which((TestLevelsFactor [[x]] %in% TrainLevelsFactor [[x]])==FALSE))
      ExcludeRows <- lapply (1:length(IndexFactor), function (j) which(TestLong [, IndexFactor [j]] %in% TestLevelsFactor [[j]] [InLevelsFactor [[j]]]))
      ExcludeRows <- do.call(c, ExcludeRows)
      
      # Convert Indices of left out test data to complete data set in long format
      ExcludeRowsConv <- vector("integer", length(ExcludeRows))
      for(j in 1:length(ExcludeRows)) {
        I1 <- sapply(1:dim(TrainLongFull) [1], function(x) TrainLongFull [x, -1] == TestLong [ExcludeRows [j], -1])
        ExcludeRowsConv [j] <- which(sapply(1:dim(I1) [2], function (x) all(I1 [, x]))==TRUE)
      }
      ExcludeRowsCensList [[i]] <- ExcludeRowsConv
      
      # Exclude rows of TestLong
      TestLong <- TestLong [-ExcludeRows, ]
    }
    oneMinuslambdaList [[i]] <- 1 - predict(CensFit, TestLong, type="response")
    
    # 7. Estimate marker values
    SurvnewFormula <- update(survModelFormula, y ~ timeInt + .)
    SurvFit <- glm (formula=SurvnewFormula, data=TrainLong, family=binomial(link=linkFunc), control=glm.control(maxit=2500))
    if(timeAsFactor) {
      TestSetExt <- cbind(TestSet, timeInt=factor(TestSet [, as.character(survModelFormula) [2] ]))
      TrainSetExt <- cbind(TrainSet, timeInt=factor(TrainSet [, as.character(survModelFormula) [2] ]))
    }
    else{
      TestSetExt <- cbind(TestSet, timeInt=TestSet [, as.character(survModelFormula) [2] ])
      TrainSetExt <- cbind(TrainSet, timeInt=TrainSet [, as.character(survModelFormula) [2] ])
    }
    
    # Exclude cases with new factor levels in test data in short format
    Check <- "error" %in% class(tryCatch(predict(SurvFit, TestSetExt), error= function (e) e))
    if(Check) {
      # Which columns are factors in test data?
      IndexFactor <- which(sapply(1:dim(TestSetExt)[2], function (x) is.factor(TestSetExt [, x]))==TRUE)
      
      # What are the levels of these factors?
      TestLevelsFactor <- sapply(IndexFactor, function (x) levels(TestSetExt [, x]))
      
      # First column does not count (censoring process)
      # What are the levels of the corresponding factors in the training data?
      TrainLevelsFactor <- sapply(IndexFactor, function (x) levels(TrainSetExt [, x]))
      
      # Which levels of the test data exists in the training data factors?
      InLevelsFactor <- lapply(1:length(TestLevelsFactor), function (x) which((TestLevelsFactor [[x]] %in% TrainLevelsFactor [[x]])==FALSE))
      ExcludeRows <- lapply (1:length(IndexFactor), function (j) which(TestSetExt [, IndexFactor [j]] %in% TestLevelsFactor [[j]] [InLevelsFactor [[j]] ]))
      ExcludeRows <- do.call(c, ExcludeRows)
      
      # Convert excluded rows of test data in short format to complete data in short format (necessary for newEvent and newTime)
      ExcludeRowsConvShort <- vector("integer", length(ExcludeRows))
      for(j in 1:length(ExcludeRows)) {
        I1 <- sapply(1:dim(dataSet) [1], function(x) dataSet [x, ] == TestSetExt [ExcludeRows [j], -dim(TestSetExt) [2] ])
        ExcludeRowsConvShort [j] <- which(sapply(1:dim(I1) [2], function (x) all(I1 [, x]))==TRUE)
      }
      ExcludeRowsDataSetList [[i]] <- ExcludeRowsConvShort
      
      # Exclude rows of test data in short format
      TestSetExt <- TestSetExt [-ExcludeRows, ]
    }
    markerList [[i]] <- predict(SurvFit, TestSetExt)
  }

  # 8. Estimate sensitivity
  # Estimate survival function of censoring process (complete data set)
  oneMinuslambda <- do.call(c, oneMinuslambdaList)
  # Merge excluded rows
  ExcludeRowsCens <- do.call(c, ExcludeRowsCensList)
  ExcludeRowsDataSet <- do.call(c, ExcludeRowsDataSetList)
  if(!is.null(ExcludeRowsCens)) {
    TrainLongFullExc <- TrainLongFull [-ExcludeRowsCens, ]
  }
  else {
    TrainLongFullExc <- TrainLongFull
  }
  G <- aggregate(oneMinuslambda ~ obj, FUN = cumprod, data = TrainLongFullExc, simplify = FALSE)
  if(!is.null(ExcludeRowsDataSet)) {
    newEvent <- dataSet [-ExcludeRowsDataSet, as.character(censModelFormula) [2]]
    newTime <- dataSet [-ExcludeRowsDataSet, as.character(survModelFormula) [2]]
  }
  else {
    newEvent <- dataSet [, as.character(censModelFormula) [2]]
    newTime <- dataSet [, as.character(survModelFormula) [2]]
  }
  n <- length(newEvent)
  if(is.null(idColumn)) {
    GT <- sapply(1:n, function(u){
      if (newTime[u] > 1)
        return(G[[2]] [u] [[1]] [newTime[u]-1]) else
          return(1) } )
  }
  else{
    GT <- sapply(1:n, function(u){
      if (newTime[u] > 1)
        return(G[[2]] [TrainLongFullExc [dataSet[u, idColumn], "obj"] ] [[1]] [newTime[u]-1]) else
          return(1) } )
  }
  # Merge markers
  marker <- do.call(c, markerList)
  RET <- sapply(marker, sens)
  tempDat <- data.frame(cutoff = sort(marker), tpr = RET [order(marker)])
  rownames(tempDat) <- 1:dim(tempDat) [1]
  
  RET <- list(Output=tempDat, 
              Input=list(timepoint=timepoint, dataSet=dataSet, trainIndices=trainIndices, 
                         survModelFormula=survModelFormula, censModelFormula=censModelFormula, 
                         linkFunc=linkFunc, idColumn=idColumn, Short=FALSE, timeAsFactor=timeAsFactor))
  class(RET) <- "discSurvTprUno"
  return(RET)
}

print.discSurvTprUno <- function (x, ...) {
  print(round(x$Output, 4))
}

plot.discSurvTprUno <- function (x, ...) {
  plot(x=x$Output [, "cutoff"], y=x$Output [, "tpr"], xlab="Cutoff", ylab="Tpr", las=1, type="l", main=paste("Tpr(c, t=", x$Input$timepoint, ")", sep=""), ...)
}

###########################
# tprUnoShort
###########################

# Description
# Computes tprUno given marker values for a general model without implicit estimation

# Input
# timepoint: Timepoint (integer scalar)
# marker: Linear predictor values of each person in the test data
# newTime: Discrete time of each person in the test data
# newEvent: Event indicator of each person in the test data
# trainTime: 
# trainEvent: 

# Output
# data.frame with columns:
  # cutoff:
  # tpr: True positive rate (numeric) \in [0, 1]

tprUnoShort <- function (timepoint, marker, newTime, newEvent, trainTime, trainEvent) {
  
  # Help function
  TransformLongToShort <- function (dataSetLong, idColumn) {
    SplitLongData <- split(dataSetLong, dataSetLong [, idColumn])
    NewDataSplit <- list()
    for(i in 1:length(SplitLongData)) {
      if(is.na(tail(SplitLongData [[i]], n=1) [,"yCens"])) {
        NewDataSplit [[i]] <- tail(SplitLongData [[i]], n=2) [-2,]
      }
      else {
        NewDataSplit [[i]] <- tail(SplitLongData [[i]], n=1)
      }
    }
    result <- do.call(rbind, NewDataSplit)
    return(result)
  }
  
  # Expand training data in long format with censoring variable
  dataSetLong <- dataLong (dataSet=data.frame(trainTime=trainTime, trainEvent=trainEvent), timeColumn="trainTime", censColumn="trainEvent")
  dataSetLongCens <- dataCensoring (dataSetLong=dataSetLong, respColumn="y", idColumn="obj")

  # Convert back to short format
  dataSetLongCensTrans <- TransformLongToShort (dataSetLong=dataSetLongCens, idColumn="obj")
  dataSetLongCensTrans <- na.omit(dataSetLongCensTrans)
  
  # Estimate nonparametric survival function of censoring variable 
  tempLifeTab <- lifeTable (dataSet=dataSetLongCensTrans, timeColumn="timeInt", censColumn="yCens")
  preG <- tempLifeTab [[1]] [, "S"]
  GT <- c(1, preG [-length(preG)])
  GT <- GT [newTime]
  
  # Help function
  sens <- function(k) {
    sensNum <- sum((marker > k) * (newTime == timepoint) * newEvent / GT, na.rm = TRUE)
    sensDenom <- sum((newTime == timepoint) * newEvent / GT, na.rm = TRUE)
    
    if (sensDenom > 0)
      return(sensNum / sensDenom) else
        return(0)
  }
  
  RET <- sapply(marker, sens)
  tempDat <- data.frame(cutoff = sort(marker), tpr = RET [order(marker)])
  rownames(tempDat) <- 1:dim(tempDat) [1]
  RET <- list(Output=tempDat, Input=list(timepoint=timepoint, marker=marker, newTime=newTime, newEvent=newEvent, trainTime=trainTime, trainEvent=trainEvent, Short=TRUE))
  class(RET) <- "discSurvTprUno"
  return(RET)
}

#############################
# fprUno
#############################

##############
# Description
# Estimates the predictive false positive rate (fpr) based on cross validation and generalized, linear models

#######
# Input
# timepoint: Discrete time interval given that the false positive rate is evaluated (integer scalar)
# dataSet: Original data. Should be in format data.frame()
# trainIndices: List of Indices from original data used for training (list of integer vectors). 
  # The length of the list is equal to the number of cross valdiation samples
# survModelFormula: Formula of the survival model
# censModelFormula: Formula of the censoring model. Normally this is done without covariates
# linkFunc: Link function of the generalized, linear model see glm
# idColumn: Name of the column with identification numbers of persons. 
# Default NULL means, that each row equals one person (no repeated measurements).

# Output
# data.frame with columns
  # cutoff: Cut off values of the linear predictor (numeric vector)
  # fpr: False positive rate (numeric vector)

fprUno <- function(timepoint, dataSet, trainIndices, survModelFormula, censModelFormula, linkFunc="logit", idColumn=NULL, timeAsFactor=TRUE) {
  
  # Input Checks
  if(length(timepoint)!=1 || !(timepoint==floor(timepoint))) {stop("Argument *timepoint* is not in the correct format! Please specify as integer scalar value.")}
  if(!is.data.frame(dataSet)) {stop("Argument *dataSet* is not in the correct format! Please specify as data.frame object.")}
  if(!is.list(trainIndices)) {stop("Argument *trainIndices* is not in the correct format! Please specify a list.")}
  InputCheck1 <- all(sapply(1:length(trainIndices), function (x) is.integer(trainIndices [[x]])))
  if(!InputCheck1) {stop("Sublists of *trainIndices* are not all integer values! Please specify a list of integer Indices.")}
  # Checks if union of test data indices equals all indices in the complete data
  if(length(trainIndices)!=1) {
    InputCheck2 <- all(sort(as.numeric(do.call(c, lapply(trainIndices, function (x) setdiff(1:dim(dataSet) [1], x)))))==(1:dim(dataSet) [1]))
  }
  else {
    InputCheck2 <- all(trainIndices [[1]]==(1:dim(dataSet) [1]))
  }
  if(!InputCheck2) {stop("Argument *trainIndices* does not contain cross validation samples! Please ensure that the union of all test indices equals the indices of the complete data set.")}
  # Formula checks
  if(!("formula" %in% class(censModelFormula))) {stop("*censModelFormula* is not of class formula! Please specify a valid formula, e. g. yCens ~ 1.")}
  if(!("formula" %in% class(survModelFormula))) {stop("*survModelFormula* is not of class formula! Please specify a valid formula, e. g. y ~ x")}
  if(!(any(names(dataSet)==idColumn) | is.null(idColumn))) {stop("Argument *idColumn* is not available in *dataSet*! Please specify the correct column name of the identification number.")}
  
#   # Help function
#   spec <- function(k){
#     
#     specNum <- sum((marker <= k) * (newTime > timepoint), na.rm = TRUE)
#     specDenom <- sum(newTime > timepoint, na.rm = TRUE)
#     
#     if (specDenom > 0)
#       return(specNum / specDenom) else
#         return(0)
#   }
  # Help function
  # Changed to avoid undefined specificity at the last time point.
  # Persons in the test data are likewise controls, 
  # if the new time point in the test data is equal to considered time point and 
  # the observation is censored
  spec <- function(k){
    
    specNum <- sum((marker <= k) * ((newTime > timepoint) | (newTime==timepoint & newEvent==0) ), na.rm = TRUE)
    specDenom <- sum(((newTime > timepoint) | (newTime==timepoint & newEvent==0) ), na.rm = TRUE)
    
    if (specDenom > 0)
      return(specNum / specDenom) else
        return(0)
  }
  
  # Loop across all training data sets
  RET <- vector("list", length(trainIndices))
  markerList <- vector("list", length(trainIndices))
  ExcludeRowsDataSetList <- vector("list", length(trainIndices))
  for(i in 1:length(trainIndices)) {
    
    # 0. Extract training, test data and responses
    TrainSet <- dataSet [trainIndices [[i]], ]
    if(length(trainIndices)!=1) {
      TestSet <- dataSet [-trainIndices [[i]], ]
    }
    else {
      TestSet <- TrainSet
    }
    
    # 1. Convert training data to long format
    if(!is.null(idColumn)) {
      TrainLong <- dataLongTimeDep (dataSet=TrainSet, timeColumn=as.character(survModelFormula) [2], 
                                    censColumn=as.character(censModelFormula) [2], idColumn=idColumn, 
                                    timeAsFactor=timeAsFactor)
    }
    else {
      TrainLong <- dataLong (dataSet=TrainSet, timeColumn=as.character(survModelFormula) [2], 
                             censColumn=as.character(censModelFormula) [2], timeAsFactor=timeAsFactor)
    }
    
    # 2. Convert response in training data to censoring variable
    TrainLong <- dataCensoring (dataSetLong=TrainLong, respColumn="y", idColumn="obj")
    
    # 3. Convert test data to long format
    if(!is.null(idColumn)) {
      TestLong <- dataLongTimeDep (dataSet=TestSet, timeColumn=as.character(survModelFormula) [2], 
                                   censColumn=as.character(censModelFormula) [2], idColumn=idColumn, timeAsFactor=timeAsFactor)
    }
    else {
      TestLong <- dataLong (dataSet=TestSet, timeColumn=as.character(survModelFormula) [2], 
                            censColumn=as.character(censModelFormula) [2], timeAsFactor=timeAsFactor)
    }

    # 7. Estimate marker values
    SurvnewFormula <- update(survModelFormula, y ~ timeInt + .)
    SurvFit <- glm (formula=SurvnewFormula, data=TrainLong, family=binomial(link=linkFunc), control=glm.control(maxit=2500))
    if(timeAsFactor) {
      TestSetExt <- cbind(TestSet, timeInt=factor(TestSet [, as.character(survModelFormula) [2] ]))
      TrainSetExt <- cbind(TrainSet, timeInt=factor(TrainSet [, as.character(survModelFormula) [2] ]))
    }
    else{
      TestSetExt <- cbind(TestSet, timeInt=TestSet [, as.character(survModelFormula) [2] ])
      TrainSetExt <- cbind(TrainSet, timeInt=TrainSet [, as.character(survModelFormula) [2] ])
    }

    Check <- "error" %in% class(tryCatch(predict(SurvFit, TestSetExt), error= function (e) e))
    if(Check) {
      # Which columns are factors in test data?
      IndexFactor <- which(sapply(1:dim(TestSetExt)[2], function (x) is.factor(TestSetExt [, x]))==TRUE)
      
      # What are the levels of these factors?
      TestLevelsFactor <- sapply(IndexFactor, function (x) levels(TestSetExt [, x]))
      
      # First column does not count (censoring process)
      # What are the levels of the corresponding factors in the training data?
      TrainLevelsFactor <- sapply(IndexFactor, function (x) levels(TrainSet [, x]))
      
      # Which levels of the test data exists in the training data factors?
      InLevelsFactor <- lapply(1:length(TestLevelsFactor), function (x) which((TestLevelsFactor [[x]] %in% TrainLevelsFactor [[x]])==FALSE))
      ExcludeRows <- lapply (1:length(IndexFactor), function (j) which(TestSetExt [, IndexFactor [j]] %in% TestLevelsFactor [[j]] [InLevelsFactor [[j]] ]))
      ExcludeRows <- do.call(c, ExcludeRows)
      
      # Convert excluded rows of test data in short format to complete data in short format (necessary for newEvent and newTime)
      ExcludeRowsConvShort <- vector("integer", length(ExcludeRows))
      for(j in 1:length(ExcludeRows)) {
        I1 <- sapply(1:dim(dataSet) [1], function(x) dataSet [x, ] == TestSetExt [ExcludeRows [j], -dim(TestSetExt) [2] ])
        ExcludeRowsConvShort [j] <- which(sapply(1:dim(I1) [2], function (x) all(I1 [, x]))==TRUE)
      }
      ExcludeRowsDataSetList [[i]] <- ExcludeRowsConvShort
      
      TestSetExt <- TestSetExt [-ExcludeRows, ]
    }
    markerList [[i]] <- predict(SurvFit, TestSetExt)
  }
  
  # 8. Estimate specificity
  marker <- do.call(c, markerList)
  ExcludeRowsDataSet <- do.call(c, ExcludeRowsDataSetList)
  if(!is.null(ExcludeRowsDataSet)) {
    newEvent <- dataSet [-ExcludeRowsDataSet, as.character(censModelFormula) [2]]
    newTime <- dataSet [-ExcludeRowsDataSet, as.character(survModelFormula) [2]]
  }
  else {
    newEvent <- dataSet [, as.character(censModelFormula) [2]]
    newTime <- dataSet [, as.character(survModelFormula) [2]]
  }
  RET <- sapply(marker, spec)
  
  # Output
  tempDat <- data.frame(cutoff = sort(marker), fpr = 1-RET [order(marker)])
  rownames(tempDat) <- 1:dim(tempDat) [1]
  RET <- list(Output=tempDat, 
              Input=list(timepoint=timepoint, dataSet=dataSet, trainIndices=trainIndices, 
                         survModelFormula=survModelFormula, censModelFormula=censModelFormula, 
                         linkFunc=linkFunc, idColumn=idColumn, Short=FALSE, timeAsFactor=timeAsFactor))
  class(RET) <- "discSurvFprUno"
  return(RET)
}

print.discSurvFprUno <- function (x, ...) {
  print(round(x$Output, 4))
}

plot.discSurvFprUno <- function (x, ...) {
  plot(x=x$Output [, "cutoff"], y=x$Output [, "fpr"], xlab="Cutoff", ylab="Fpr", las=1, type="l", main=paste("Fpr(c, t=", x$Input$timepoint, ")", sep=""), ...)
}

#######################
# fprUnoShort
#######################

# Description
# Estimates the false positive rate given prior estimated marker values

fprUnoShort <- function (timepoint, marker, newTime, newEvent) {
  
  # Help function
#   spec <- function(k) {
#     
#     specNum <- sum((marker <= k) * (newTime > timepoint), na.rm = TRUE)
#     specDenom <- sum(newTime > timepoint, na.rm = TRUE)
#     
#     if (specDenom > 0)
#       return(specNum / specDenom) else
#         return(0)
#   }
  spec <- function(k){
    
    specNum <- sum((marker <= k) * ((newTime > timepoint) | (newTime==timepoint & newEvent==0) ), na.rm = TRUE)
    specDenom <- sum(((newTime > timepoint) | (newTime==timepoint & newEvent==0) ), na.rm = TRUE)
    
    if (specDenom > 0)
      return(specNum / specDenom) else
        return(0)
  }
  
  # Output
  RET <- sapply(marker, spec)
  tempDat <- data.frame(cutoff = sort(marker), fpr = 1 - RET [order(marker)])
  rownames(tempDat) <- 1:dim(tempDat) [1]
  RET <- list(Output=tempDat, Input=list(timepoint=timepoint, marker=marker, newTime=newTime, Short=TRUE))
  class(RET) <- "discSurvFprUno"
  return(RET)
}

######################
# aucUno
######################

# Description
# Computes the auc (area under the curve) measure given timepoint t
# Both input object must have identical input parameters, e. g. same timepoint, formula!

# Input
# tprObj: Object of class "discSurvTprUno"
# fprObj: Object of class "discSurvFprUno"

# Output
# auc value (numeric scalar) given timepoint of *tprObj* and *fprObj*

aucUno <- function (tprObj, fprObj) {
  
  # Input checks
  if(class(tprObj)!="discSurvTprUno") {stop("This object has not the appropriate class! Please specify an object of class *discSurvTprUno*.")}
  if(class(fprObj)!="discSurvFprUno") {stop("This object has not the appropriate class! Please specify an object of class *discSurvFprUno*.")}
  if(tprObj$Input$Short!=fprObj$Input$Short) {stop("Tpr and fpr were computed using different functions! Please ensure that both are estimated either by the cross validated version or the short version.")}
  if(!tprObj$Input$Short) {
    InputCheck <- identical(tprObj$Input, fprObj$Input)
    if(!InputCheck) {stop("Some input parameters of *tprObj* or *fprObj* are not identical! Please check if both objects were estimated using exact identical input values.")}
  }
  else {
    InputCheck1 <- identical(tprObj$Input$timepoint, fprObj$Input$timepoint)
    InputCheck2 <- identical(tprObj$Input$marker, fprObj$Input$marker)
    InputCheck3 <- identical(tprObj$Input$newTime, fprObj$Input$newTime)
    InputCheck <- all(InputCheck1, InputCheck2, InputCheck3)
    if(!InputCheck) {stop("Some input parameters of *tprObj* or *fprObj* are not identical! Please check if both objects were estimated using exact identical input values.")}
  }
  
  tpr <- c(1, tprObj$Output$tpr)
  fpr <- c(1, fprObj$Output$fpr)
  
  trapz <- function (x, y){ # from package caToosl
    idx = 2:length(x)
    return(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2)
  }
  
  Output <- - trapz(fpr, tpr)
  names(Output) <- paste("AUC(t=", tprObj$Input$timepoint, ")", sep="")
  auc <- list(Output = Output, Input=list(tprObj=tprObj, fprObj=fprObj))
  class(auc) <- "discSurvAucUno"
  return(auc)
}

print.discSurvAucUno <- function (x, ...) {
  print(round(x$Output, 4))
}

plot.discSurvAucUno <- function (x, ...) {
  tprVal <- x$Input$tprObj$Output [, "tpr"]
  fprVal <- x$Input$fprObj$Output [, "fpr"]
  plot(x=fprVal, y=tprVal, xlab="Fpr", ylab="Tpr", las=1, type="l", main=paste("ROC(c, t=", x$Input$tprObj$Input$timepoint, ")", sep=""), ...)
  lines(x=seq(0, 1, length.out=500), y=seq(0, 1, length.out=500), lty=2)
}

########################################
# concordanceIndex
########################################

# Description
# Calculates the Concordance index (independent measure of time)

# Input
# auc: Double vector with length equal to the number of time intervals
# SurvT: Double vector of estimated survival curve P(T>t) for each time interval

# Output
# Weighted integrated auc over time as measure of accuracy. The format is a double scalar value.

concorIndex <- function (aucObj) {

  # Input checks
  if(class(aucObj)!="discSurvAucUno") {stop("This object has not the appropriate class! Please specify an object of class *discSurvAucUno*.")}

  if(aucObj$Input$tprObj$Input$Short) {
    
    # Estimate AUC for all t
    marker <- aucObj$Input$tprObj$Input$marker
    newTime <- aucObj$Input$tprObj$Input$newTime
    newEvent <- aucObj$Input$tprObj$Input$newEvent
    trainTime <- aucObj$Input$tprObj$Input$trainTime
    trainEvent <- aucObj$Input$tprObj$Input$trainEvent
    MaxTime <- max(trainTime)
    AUCalltime <- vector("numeric", MaxTime)
    for(i in 1:MaxTime) {
      tempTPR <- tprUnoShort (timepoint=i, marker=marker, newTime=newTime, newEvent=newEvent, 
                              trainTime=trainTime, trainEvent=trainEvent)
      tempFPR <- fprUnoShort (timepoint=i, marker=marker, newTime=newTime, newEvent=newEvent)
      AUCalltime [i] <- as.numeric(aucUno (tprObj=tempTPR, fprObj=tempFPR)$Output)
      cat("Timepoint =", i, "done", "\n")
    }

    # Estimate nonparametric survival function S(T=t) and marginal probabilities P(T=t)
    tempLifeTab <- lifeTable (dataSet=data.frame(trainTime=trainTime, trainEvent=trainEvent), timeColumn="trainTime", censColumn="trainEvent")
    MargHaz <- tempLifeTab [[1]] [, "hazard"]
    MargSurv <- estSurv(MargHaz)
    MargProb <- estMargProb(MargHaz)
  }
  
  else {

    # Estimate AUC for all t
    MaxTime <- max(aucObj$Input$tprObj$Input$dataSet [, as.character(aucObj$Input$tprObj$Input$survModelFormula) [2] ])
    DataSet <- aucObj$Input$tprObj$Input$dataSet
    TrainIndices <- aucObj$Input$tprObj$Input$trainIndices
    SurvModelFormula <- aucObj$Input$tprObj$Input$survModelFormula
    CensModelFormula <- aucObj$Input$tprObj$Input$censModelFormula
    LinkFunc <- aucObj$Input$tprObj$Input$linkFunc
    IdColumn <- aucObj$Input$tprObj$Input$idColumn
    timeAsFactor <- aucObj$Input$tprObj$Input$timeAsFactor
    AUCalltime <- vector("numeric", MaxTime)
    for(i in 1:MaxTime) {
      tempTPR <- tprUno (timepoint=i, dataSet=DataSet, trainIndices=TrainIndices, survModelFormula=SurvModelFormula, censModelFormula=CensModelFormula, linkFunc=LinkFunc, idColumn=IdColumn, timeAsFactor=timeAsFactor)
      tempFPR <- fprUno (timepoint=i, dataSet=DataSet, trainIndices=TrainIndices,  survModelFormula=SurvModelFormula, censModelFormula=CensModelFormula, linkFunc=LinkFunc, idColumn=IdColumn, timeAsFactor=timeAsFactor)
      AUCalltime [i] <- aucUno (tprObj=tempTPR, fprObj=tempFPR)$Output
      cat("Timepoint =", i, "done", "\n")
    }
  
    # Estimate survival function and marginal probabilities without covariates
    if(!is.null(IdColumn)) {
      TrainLongFull <- dataLongTimeDep (dataSet=DataSet, timeColumn=as.character(SurvModelFormula) [2], 
                                        censColumn=as.character(CensModelFormula) [2], idColumn=IdColumn, timeAsFactor=timeAsFactor)
    }
    else {
      TrainLongFull <- dataLong (dataSet=DataSet, timeColumn=as.character(SurvModelFormula) [2], 
                                 censColumn=as.character(CensModelFormula) [2], timeAsFactor=timeAsFactor)
    }
    
    # Estimate marginal survival probability with glm
    MargFormula <- y ~ timeInt
    MargFit <- glm (formula=MargFormula, data=TrainLongFull, family=binomial(link=LinkFunc), control=glm.control(maxit=2500))
    if(timeAsFactor) {
      PredMargData <- data.frame(timeInt=factor(min(TrainLongFull [, as.character(SurvModelFormula) [2] ]):max(TrainLongFull [, as.character(SurvModelFormula) [2] ])))
    }
    else{
      PredMargData <- data.frame(timeInt=min(TrainLongFull [, as.character(SurvModelFormula) [2] ]):max(TrainLongFull [, as.character(SurvModelFormula) [2] ]))
    }
    MargHaz <- as.numeric(predict(MargFit, PredMargData, type="response"))
    # Survival function S(T=t) = P(T>t) = \prod_{j=1}^{t1} (1-\lambda (T=j))
    MargSurv <- estSurv(MargHaz)
    # Marginal Probability P(T=t) = \lambda (T=t) \prod_{j=1}^{t-1} (1-\lambda (T=j))
    MargProb <- estMargProb(MargHaz)
  }
  
  # Calcualte concordance index
  weights <- MargProb*MargSurv / sum(MargProb*MargSurv)
  # Last weight is zero and can therefore be omitted
  Concor <- sum(AUCalltime * weights [-length(weights)])
  names(Concor) <- "C*"
  names(AUCalltime) <- paste("AUC(t=", 1:MaxTime, "|x)", sep="")
  Output <- list(Output=Concor, Input=list(aucObj=aucObj, AUC=AUCalltime, MargProb=MargProb, MargSurv=MargSurv))
  class(Output) <- "discSurvConcorIndex"
  return(Output)
}

print.discSurvConcorIndex <- function (x, ...) {
  print(round(x$Output, 4))
}

summary.discSurvConcorIndex <- function (object, ...) {
  cat("Concordance: Should be higher than 0.5 (random assignment)", "\n")
  print(round(object$Output, 4))
  cat("\n", "AUC(t): Should be higher than 0.5 for all time points (random assignment)", "\n")
  print(round(object$Input$AUC, 4))
  cat("\n", "Marginal P(T=t) without covariates (used in weighting)", "\n")
  print(round(object$Input$MargProb, 4))
  cat("\n", "Marginal S(T=t) without covariates (used in weighting)", "\n")
  print(round(object$Input$MargSurv, 4))
}

######################
# predErrDisc

# # 1. Estimate survival function of model
# # 2. Estimate censoring survival function of a covariate free model
# # 3. Calculate observed survival function
# # 4. Calculate prediction error curve given timepoint t
# 
# Problems with this function:
# - Survival function for censoring process is only evaluated with training data!
# - 
# - Function is too complex! Should be simplified! 
# 
# predErrDisc <- function(timepoints, dataSet, trainIndices, survModelFormula, censModelFormula, linkFunc="logit", idColumn=NULL) {
#   
#   # Input Checks
#   if(!is.vector(timepoints) | !all(timepoints==floor(timepoints))) {stop("Argument *timepoints* is not in the correct format! Please specify as integer vector.")}
#   if(!is.data.frame(dataSet)) {stop("Argument *dataSet* is not in the correct format! Please specify as data.frame object.")}
#   if(!is.list(trainIndices)) {stop("Argument *trainIndices* is not in the correct format! Please specify a list.")}
#   InputCheck1 <- all(sapply(1:length(trainIndices), function (x) is.integer(trainIndices [[x]])))
#   if(!InputCheck1) {stop("Sublists of *trainIndices* are not all integer values! Please specify a list of integer Indices.")}
#   if(length(trainIndices)!=1) {
#     InputCheck2 <- all(sort(as.numeric(do.call(c, lapply(trainIndices, function (x) setdiff(1:dim(dataSet) [1], x)))))==(1:dim(dataSet) [1]))
#   }
#   else {
#     InputCheck2 <- all(trainIndices [[1]]==(1:dim(dataSet) [1]))
#   }
#   if(!InputCheck2) {stop("Argument *trainIndices* does not contain cross validation samples! Please ensure that the union of all test indices equals the indices of the complete data set.")}
#   if(!("formula" %in% class(censModelFormula))) {stop("*censModelFormula* is not of class formula! Please specify a valid formula, e. g. yCens ~ 1")}
#   if(!("formula" %in% class(survModelFormula))) {stop("*survModelFormula* is not of class formula! Please specify a valid formula, e. g. y ~ x")}
#   if(!(any(names(dataSet)==idColumn) | is.null(idColumn))) {stop("Argument *idColumn* is not available in *dataSet*! Please specify the correct column name of the identification number.")}
#   
#   Help functions
#     ObsSurvFunc <- function (x) {
#       resp <- TrainLongFullExcSplit [[x]]
#       if(resp [length(resp)]==1) {
#         temp <- c(rep(1, length(resp)-1), 0)
#         result <- tryCatch(temp [timepoints [k] ], error=function (e) NA)
#         return(result)
#       }
#       else {
#         temp <- rep(1, length(resp))
#         result <- tryCatch(temp [timepoints [k] ], error=function (e) NA)
#         return (result)
#       }
#     }
#   
#   WeightFunction <- function () {
#     PartialSum1 <- newEvent * (1 - Sobs) / GT
#     PartialSum2 <- Sobs / GTfixed
#     return(PartialSum1 + PartialSum2)
#   }
#   
#   predErr <- function () {
#     sum(weights [[k]] * (STfixed - Sobs)^2, na.rm=TRUE) / length(weights [[k]] [!is.na(weights [[k]])])
#   }
#   
#   # Loop across all training data sets
#   ExcludeRowsCensList <- vector("list", length(trainIndices))
#   ExcludeRowsDataSetList <- vector("list", length(trainIndices))
#   oneMinusLambdaList <- vector("list", length(trainIndices))
#   oneMinusLambdaSurvList <- vector("list", length(trainIndices))
#   
#   # Convert full sample to long format
#   if(!is.null(idColumn)) {
#     TrainLongFull <- dataLongTimeDep (dataSet=dataSet, timeColumn=as.character(survModelFormula) [2], censColumn=as.character(censModelFormula) [2], idColumn=idColumn)
#   }
#   else {
#     TrainLongFull <- dataLong (dataSet=dataSet, timeColumn=as.character(survModelFormula) [2], censColumn=as.character(censModelFormula) [2])
#   }
#   
#   for(i in 1:length(trainIndices)) {
#     
#     # 0. Extract training, test data and responses
#     TrainSet <- dataSet [trainIndices [[i]], ]
#     if(length(trainIndices)!=1) {
#       TestSet <- dataSet [-trainIndices [[i]], ]
#     }
#     else {
#       TestSet <- TrainSet
#     }
#     
#     # 1. Convert training data to long format
#     if(!is.null(idColumn)) {
#       TrainLong <- dataLongTimeDep (dataSet=TrainSet, timeColumn=as.character(survModelFormula) [2], censColumn=as.character(censModelFormula) [2], idColumn=idColumn)
#     }
#     else {
#       TrainLong <- dataLong (dataSet=TrainSet, timeColumn=as.character(survModelFormula) [2], censColumn=as.character(censModelFormula) [2])
#     }
#     
#     # 2. Convert response in training data to censoring variable
#     TrainLong <- dataCensoring (dataSetLong=TrainLong, respColumn="y", idColumn="obj")
#     
#     # 3. Convert test data to long format
#     if(!is.null(idColumn)) {
#       TestLong <- dataLongTimeDep (dataSet=TestSet, timeColumn=as.character(survModelFormula) [2], censColumn=as.character(censModelFormula) [2], idColumn=idColumn)
#     }
#     else {
#       TestLong <- dataLong (dataSet=TestSet, timeColumn=as.character(survModelFormula) [2], censColumn=as.character(censModelFormula) [2])
#     }
#     
#     # 4. Fit censoring model on training data in long format
#     CensnewFormula <- update(censModelFormula, yCens ~ timeInt + .)
#     CensFit <- glm (formula=CensnewFormula, data=TrainLong, family=binomial(link=linkFunc), control=glm.control(maxit=2500))
#     
#     # 5. Estimate censoring survival curves inputs on test data
#     # Exclude cases with new factor levels in test data in long format
#     Check <- "error" %in% class(tryCatch(predict(CensFit, TestLong, type="response"), error= function (e) e))
#     if(Check) {
#       
#       # Which columns are factors in test data?
#       IndexFactor <- which(sapply(1:dim(TestLong)[2], function (x) is.factor(TestLong [, x]))==TRUE)
#       
#       # What are the levels of these factors?
#       TestLevelsFactor <- sapply(IndexFactor, function (x) levels(TestLong [, x]))
#       
#       # First column does not count (censoring process)
#       # What are the levels of the corresponding factors in the training data?
#       TrainLevelsFactor <- sapply(IndexFactor, function (x) levels(TrainLong [, x+1]))
#       
#       # Which levels of the test data exists in the training data factors?
#       InLevelsFactor <- lapply(1:length(TestLevelsFactor), function (x) which((TestLevelsFactor [[x]] %in% TrainLevelsFactor [[x]])==FALSE))
#       ExcludeRows <- lapply (1:length(IndexFactor), function (j) which(TestLong [, IndexFactor [j]] %in% TestLevelsFactor [[j]] [InLevelsFactor [[j]]]))
#       ExcludeRows <- do.call(c, ExcludeRows)
#       
#       # Convert Indices of left out test data to complete data set in long format
#       ExcludeRowsConv <- vector("integer", length(ExcludeRows))
#       for(j in 1:length(ExcludeRows)) {
#         I1 <- sapply(1:dim(TrainLongFull) [1], function(x) TrainLongFull [x, -1] == TestLong [ExcludeRows [j], -1])
#         ExcludeRowsConv [j] <- which(sapply(1:dim(I1) [2], function (x) all(I1 [, x]))==TRUE)
#       }
#       ExcludeRowsCensList [[i]] <- ExcludeRowsConv
#       
#       # Exclude rows of TestLong
#       TestLong <- TestLong [-ExcludeRows, ]
#     }
#     oneMinusLambdaList [[i]] <- 1 - predict(CensFit, TestLong, type="response")
#     
#     # 7. Estimate survival function inputs on test data
#     SurvnewFormula <- update(survModelFormula, y ~ timeInt + .)
#     SurvFit <- glm (formula=SurvnewFormula, data=TrainLong, family=binomial(link=linkFunc), control=glm.control(maxit=2500))
#     oneMinusLambdaSurvList [[i]] <- 1 - predict(SurvFit, TestLong, type="response")
#   }
#   
#   # 8. Estimate prediction error curve
#   # Estimate survival function of censoring process (complete data set)
#   oneMinuslambda <- do.call(c, oneMinusLambdaList)
#   # Merge excluded rows
#   ExcludeRowsCens <- do.call(c, ExcludeRowsCensList)
#   ExcludeRowsDataSet <- do.call(c, ExcludeRowsDataSetList)
#   if(!is.null(ExcludeRowsCens)) {
#     TrainLongFullExc <- TrainLongFull [-ExcludeRowsCens, ]
#   }
#   else {
#     TrainLongFullExc <- TrainLongFull
#   }
#   G <- aggregate(oneMinuslambda ~ obj, FUN = cumprod, data = TrainLongFullExc, simplify = FALSE)
#   if(!is.null(ExcludeRowsDataSet)) {
#     newEvent <- dataSet [-ExcludeRowsDataSet, as.character(censModelFormula) [2]]
#     newTime <- dataSet [-ExcludeRowsDataSet, as.character(survModelFormula) [2]]
#   }
#   else {
#     newEvent <- dataSet [, as.character(censModelFormula) [2]]
#     newTime <- dataSet [, as.character(survModelFormula) [2]]
#   }
#   n <- length(newEvent)
#   GT <- sapply(1:n, function(u){
#     if (newTime[u] > 1)
#       return(G[[2]] [u] [[1]] [newTime[u] - 1]) else
#         return(1) } )
#   
#   # Iterate over all time points
#   Lengthtimepoints <- length(timepoints)
#   weights <- vector("list", Lengthtimepoints)
#   results <- vector("numeric", Lengthtimepoints)
#   for(k in 1:Lengthtimepoints) {
#   
#   GTfixed <- sapply(1:n, function(u){
#   if (timepoints [k] > 1)
#     return(G[[2]] [u] [[1]] [timepoints [k] ]) else
#       return(1) } )
# 
#   # Estimate survival function
#   oneMinuslambdaSurv <- do.call(c, oneMinusLambdaSurvList)
#   S <- aggregate(oneMinuslambdaSurv ~ obj, FUN = cumprod, data = TrainLongFullExc, simplify = FALSE)
#   STfixed <- sapply(1:n, function(u){
#   if (timepoints [k] > 1)
#     return(S[[2]] [u] [[1]] [timepoints [k] ]) else
#       return(1) } )
#   
#   # Observed survival function must be evaluated on test data!
#   # Calculate observed survival function (complete data)
#   TrainLongFullExcSplit <- split(TrainLongFullExc$y, TrainLongFullExc$obj)
#   Sobs <- lapply(1:length(TrainLongFullExcSplit), ObsSurvFunc)
#   Sobs <- do.call(c, Sobs)
# 
#   # Exclude cases for which observed survival function is not available
#   IndexExclusion <- is.na(Sobs) | is.nan(Sobs) | is.infinite(Sobs)
#   if(!all(IndexExclusion==FALSE)) {
#     Sobs <- Sobs [!IndexExclusion]
#     STfixed <- STfixed [!IndexExclusion]
#     GTfixed <- GTfixed [!IndexExclusion]
#     GT <- GT [!IndexExclusion]
#     newEvent <- newEvent [!IndexExclusion]
#   }
#   
#   # Calculate prediction error curve given timepoints t
#   weights [[k]] <- WeightFunction ()
#   results [k] <- predErr ()
#   }
#   
#   # Combine outputs
#   RET <- list(Output=list(predErr = results, weights = weights), 
#               Input=list(timepoints=timepoints, dataSet=dataSet, trainIndices=trainIndices, survModelFormula=survModelFormula, censModelFormula=censModelFormula, linkFunc=linkFunc, idColumn=idColumn, Short=FALSE))
#   class(RET) <- "discSurvPredErrDisc"
#   return(RET)
# }

print.discSurvPredErrDisc <- function (x, ...) {
  print(round(x$Output$predErr, 4))
}

summary.discSurvPredErrDisc <- function (object, ...) {
  cat("Prediction error curve: Should be lower than 0.25 (random assignment) for all timepoints", "\n")
  print(round(object$Output$predErr, 4))
  cat("Corresponding weights given by the censoring survival function", "\n")
  tempDat <- lapply(1:length(object$Output$weights), function (x) {round(object$Output$weights [[x]], 4)})
  names(tempDat) <- paste("Time interval = ", 1:length(object$Output$weights), sep="")
  print(tempDat)
}

#####################################
# predErrDiscShort
# Short Version without CV but usable for arbitrary models

predErrDiscShort <- function (timepoints, estSurvList, newTime, newEvent, trainTime, trainEvent) {
  
  # Help functions
  WeightFunction <- function () {
    PartialSum1 <- newEventTemp * (1 - Sobs) / GT
    PartialSum2 <- Sobs / GTfixed
    return(PartialSum1 + PartialSum2)
  }
  predErr <- function () {
    sum(weights * (estSurvInd - Sobs)^2, na.rm=TRUE) / length(weights [!is.na(weights)])
  }
  
  #####################
  # Execution
  
  # Expand training data in long format with censoring variable
  dataSetLong <- dataLong (dataSet=data.frame(trainTime=trainTime, trainEvent=trainEvent), timeColumn="trainTime", censColumn="trainEvent")
  dataSetLongCens <- dataCensoring (dataSetLong=dataSetLong, respColumn="y", idColumn="obj")
  dataSetLongCens <- na.omit(dataSetLongCens)
  
  # Estimate glm with no covariates of censoring process
  glmCovariateFree <- glm(yCens ~ timeInt, data=dataSetLongCens, family=binomial(), control=glm.control(maxit = 2500))
  # Predict survival function of censoring process
  factorPrep <- factor(1:max(as.numeric(as.character(dataSetLongCens$timeInt))))
  GT_est <- cumprod(1 - predict(glmCovariateFree, newdata=data.frame(timeInt=factorPrep), type="response"))
  
  # Loop over all time points
  predErrorValues <- vector("numeric", length(timepoints))
  StoreWeights <- vector("list", length(timepoints))
  for(k in 1:length(timepoints)) {
    
    # Instable!
    # Restrict newTime and newEvent
    #   Check <- newTime < timepoints [k] & newEvent == 0
    #   newTimeTemp <- newTime [!Check]
    #   newEventTemp <- newEvent [!Check]
    #   if(length(newTimeTemp)==0) {
    #     return(0)
    #   }
    
    newTimeTemp <- newTime
    newEventTemp <- newEvent
    
    # Estimate survival functions of censoring process
    GT <- c(1, GT_est)
    GT <- GT [newTimeTemp]
    GTfixed <- GT_est [timepoints [k] ]
    
    # Generate observed survival function
    Sobs <- ifelse(timepoints [k] < newTimeTemp, 1, 0)
    
    # Filter all predicted values: First order is timepoint and second layer is individuals
    estSurvInd <- sapply(1:length(estSurvList), function (x) estSurvList [[x]] [timepoints [k] ])
    # estSurvInd <- estSurvInd [!Check]
    
    # Estimate weights of each person in test data
    weights <- WeightFunction ()
    StoreWeights [[k]] <- weights
    
    # Estimate prediction error
    predErrorValues [k] <- predErr ()
  }
  names(predErrorValues) <- paste("T=", timepoints, sep="")
  CheckRemove <- is.infinite(predErrorValues) | is.nan(predErrorValues) | is.na (predErrorValues)
  predErrorValues <- predErrorValues [!CheckRemove]
  
  # Combine outputs
  RET <- list(Output=list(predErr = predErrorValues, weights = StoreWeights), 
              Input=list(timepoints=timepoints, estSurvList=estSurvList, newTime=newTime, newEvent=newEvent, trainTime=trainTime, trainEvent=trainEvent, Short=TRUE))
  class(RET) <- "discSurvPredErrDisc"
  return(RET)
}

#######################
# intPredErrDisc

# Description
# Estimates prediction error curves

intPredErrDisc <- function (predErrObj, tmax=NULL) {
  
  # Input check
  if(!(class(predErrObj)=="discSurvPredErrDisc")) {stop("Object *predErrObj* is not of class *discSurvPredErrDisc*! Please give an appropriate objecte type as input.")}
  
#  if(predErrObj$Input$Short==FALSE) {
#  
#   # Help function
#   predErrDiscTime <- function (t) {
#     predErrDisc (timepoints=t, 
#                  dataSet=predErrObj$Input$dataSet, trainIndices=predErrObj$Input$trainIndices, survModelFormula=predErrObj$Input$survModelFormula, censModelFormula=predErrObj$Input$censModelFormula, linkFunc=predErrObj$Input$linkFunc, idColumn=predErrObj$Input$idColumn)
#   }
# 
#   # Estimate marginal probabilities P(T=t)
#   MaxTime <- max(predErrObj$Input$dataSet [, as.character(predErrObj$Input$survModelFormula) [2] ])
#   DataSet <- predErrObj$Input$dataSet
#   TrainIndices <- predErrObj$Input$trainIndices
#   SurvModelFormula <- predErrObj$Input$survModelFormula
#   CensModelFormula <- predErrObj$Input$censModelFormula
#   LinkFunc <- predErrObj$Input$linkFunc
#   IdColumn <- predErrObj$Input$idColumn
# 
#   # Estimate survival function and marginal probabilities without covariates
#   if(!is.null(IdColumn)) {
#     TrainLongFull <- dataLongTimeDep (dataSet=DataSet, timeColumn=as.character(SurvModelFormula) [2], censColumn=as.character(CensModelFormula) [2], idColumn=IdColumn)
#   }
#   else {
#     TrainLongFull <- dataLong (dataSet=DataSet, timeColumn=as.character(SurvModelFormula) [2], censColumn=as.character(CensModelFormula) [2])
#   }
#   MargFormula <- y ~ timeInt
#   MargFit <- glm (formula=MargFormula, data=TrainLongFull, family=binomial(link=LinkFunc), control=glm.control(maxit=2500))
#   PredMargData <- data.frame(timeInt=factor(min(TrainLongFull [, as.character(SurvModelFormula) [2] ]):max(TrainLongFull [, as.character(SurvModelFormula) [2] ])))
#   MargHaz <- as.numeric(predict(MargFit, PredMargData, type="response"))
#   # Marginal Probability P(T=t) = \lambda (T=t) \prod_{j=1}^{t-1} (1-\lambda (T=j))
#   MargProbs <- estMargProb(MargHaz)
# 
#   # Estimate prediction error curve over all time points
#   PredsErrorCurve <- predErrDiscTime (1:(MaxTime+1))$Output$predErr
#   IncludedIndices <- !(is.na(PredsErrorCurve) | is.nan (PredsErrorCurve) | is.infinite(PredsErrorCurve))
#   PredsErrorCurve <- PredsErrorCurve [IncludedIndices]
#   MargProbs <- as.numeric(MargProbs [IncludedIndices])
#   
#   # Output
#   Result <- sum(PredsErrorCurve * MargProbs)
#   return(c(IntPredErr=Result))
#   }
#  else {
    # Help function
    predErrDiscTime <- function (t) {
      predErrDiscShort (timepoints= t, estSurvList=EstSurvList, newTime=NewTime, newEvent=NewEvent, trainTime=TrainTime, trainEvent=TrainEvent)
    }
    
    # Estimate marginal probabilities P(T=t)
    MaxTime <- max(predErrObj$Input$trainTime)
    if(!is.null(tmax)) {
      if(tmax <= MaxTime) {
        MaxTime <- tmax
      }
      else {
        warning("Argument *tmax* is higher than the latest observed interval in training data. 
                Only prediction errors up to the latest observed interval time are given.")
      }
    }
    EstSurvList <- predErrObj$Input$estSurvList
    NewTime <- predErrObj$Input$newTime
    NewEvent <- predErrObj$Input$newEvent
    TrainTime <- predErrObj$Input$trainTime
    TrainEvent <- predErrObj$Input$trainEvent
    
    # Estimate survival function and marginal probabilities without covariates
    TrainLongFull <- dataLong (dataSet=data.frame(TrainTime=TrainTime, TrainEvent=TrainEvent), timeColumn="TrainTime", censColumn="TrainEvent")
    MargFormula <- y ~ timeInt
    MargFit <- glm (formula=MargFormula, data=TrainLongFull, family=binomial(), control=glm.control(maxit=2500))
    PredMargData <- data.frame(timeInt=factor(min(TrainTime):max(TrainTime)))
    MargHaz <- as.numeric(predict(MargFit, PredMargData, type="response"))
    # Marginal Probability P(T=t) = \lambda (T=t) \prod_{j=1}^{t-1} (1-\lambda (T=j))
    MargProbs <- estMargProb(MargHaz)
    
    # Estimate prediction error curve over all time points
    PredsErrorCurve <- predErrDiscTime (1:(MaxTime+1))$Output$predErr
    IncludedIndices <- !(is.na(PredsErrorCurve) | is.nan (PredsErrorCurve) | is.infinite(PredsErrorCurve))
    PredsErrorCurve <- PredsErrorCurve [IncludedIndices]
    MargProbs <- as.numeric(MargProbs [IncludedIndices])
    
    # Enlarge PredsErrorCurve with 0 to match length of MargProbs
    if(length(PredsErrorCurve)!=length(MargProbs)) {
      PredsErrorCurve <- c(PredsErrorCurve, rep(0, length(MargProbs) - length(PredsErrorCurve)))
    }
    
    # Output
    Result <- sum(PredsErrorCurve [1:MaxTime] * MargProbs [1:MaxTime]) / sum(MargProbs [1:MaxTime])
    return(c(IntPredErr=Result))
#  }
}

######################
# martingaleResid

# Description
# Calculates the martingale residuals

martingaleResid <- function (dataSet, survModelFormula, censColumn, linkFunc="logit", idColumn=NULL) {
  
  # Input checks
  if(!is.data.frame(dataSet)) {stop("Argument *dataSet* is not in the correct format! Please specify as data.frame object.")}
  if(!("formula" %in% class(survModelFormula))) {stop("*survModelFormula* is not of class formula! Please specify a valid formula, e. g. y ~ x + z.")}
  if(!any(names(dataSet)==censColumn)) {stop("Argument *censColumn* is not available in *dataSet*! Please specify the correct column name of the event indicator.")}
  if(!(any(names(dataSet)==idColumn) | is.null(idColumn))) {stop("Argument *idColumn* is not available in *dataSet*! Please specify the correct column name of the identification numbers of persons.")}
  
  # Convert to long format
  if(!is.null(idColumn)) {
    dataSetLong <- dataLongTimeDep (dataSet=dataSet, timeColumn=as.character(survModelFormula) [2], censColumn=censColumn, idColumn=idColumn)
  }
  else {
    dataSetLong <- dataLong (dataSet=dataSet, timeColumn=as.character(survModelFormula) [2], censColumn=censColumn)
  }
  
  # Fit generalized, linear model
  NewFormula <- update(survModelFormula, y ~ timeInt + .)
  glmFit <- glm(formula=NewFormula, data=dataSetLong, family=binomial(link=linkFunc), control=glm.control(maxit=2500))
  hazards <- predict(glmFit, type="response")
  
  # Calculate residuals
  splitHazards <- split(hazards, dataSetLong$obj)
  splitY <- split(dataSetLong$y, dataSetLong$obj)
  martResid <- sapply(1:length(splitY), function (x) sum(splitY [[x]] - splitHazards [[x]]))
  Output <- list(Output=list(MartingaleResid=martResid, GlmFit=glmFit), 
                 Input=list(dataSet=dataSet, survModelFormula=survModelFormula, censColumn=censColumn, linkFunc=linkFunc, idColumn=idColumn))
  class(Output) <- "discSurvMartingaleResid"
  return(Output)
}

print.discSurvMartingaleResid <- function (x, ...) {
  print(round(x$Output$MartingaleResid, 4))
}

plot.discSurvMartingaleResid <- function (x, ...) {
  
  # Convert to long format
  if(!is.null(x$Input$idColumn)) {
    dataSetLong <- dataLongTimeDep (dataSet=x$Input$dataSet, timeColumn=as.character(x$Input$survModelFormula) [2], censColumn=x$Input$censColumn, idColumn=x$Input$idColumn)
  }
  else {
    dataSetLong <- dataLong (dataSet=x$Input$dataSet, timeColumn=as.character(x$Input$survModelFormula) [2], censColumn=x$Input$censColumn)
  }
  
  splitDataSetLong <- split(dataSetLong, dataSetLong$obj)
  tailSplitDataSetLong <- lapply(splitDataSetLong, function (x) tail(x, 1))
  tailSplitDataSetLong <- do.call(rbind, tailSplitDataSetLong)
  
  # Covariates of interest
  LengthSurvFormula <- length(attr(terms(x$Input$survModelFormula), "term.labels"))
  CovarSurvFormula <- attr(terms(x$Input$survModelFormula), "term.labels")
  
  # Create plots of covariates
  for(i in 1:LengthSurvFormula) {
    tempData <- data.frame(x=tailSplitDataSetLong [, CovarSurvFormula [i] ], y=x$Output$MartingaleResid)
    tempData <- tempData [order(tempData$x), ]
    plot(x=tempData$x, y=tempData$y, las=1, xlab=CovarSurvFormula [i], ylab="Martingale Residuals", ...)
    if(is.numeric(tempData$x)) {
      loessPred <- predict(loess(formula=y ~ x, data=tempData))
      lines(x=tempData$x, y=loessPred)
    }
    abline(h=0, lty=2)
  }
}

########################
# devResid

# Description
# Computes the root of the squared deviance residual

devResid <- function (dataSet, survModelFormula, censColumn, linkFunc="logit", idColumn=NULL) {
  
  # Input checks
  if(!is.data.frame(dataSet)) {stop("Argument *dataSet* is not in the correct format! Please specify as data.frame object.")}
  if(!("formula" %in% class(survModelFormula))) {stop("*survModelFormula* is not of class formula! Please specify a valid formula, e. g. y ~ x + z.")}
  if(!any(names(dataSet)==censColumn)) {stop("Argument *censColumn* is not available in *dataSet*! Please specify the correct column name of the event indicator.")}
  if(!(any(names(dataSet)==idColumn) | is.null(idColumn))) {stop("Argument *idColumn* is not available in *dataSet*! Please specify the correct column name of the identification numbers of persons.")}
  
  # Help function
  SquDevResid <- function (x) {-2*sum(splitY [[x]] * log(splitHazards [[x]]) + (1 - splitY [[x]]) * log(1 - splitHazards [[x]] ))}
  
  # Convert to long format
  if(!is.null(idColumn)) {
    dataSetLong <- dataLongTimeDep (dataSet=dataSet, timeColumn=as.character(survModelFormula) [2], censColumn=censColumn, idColumn=idColumn)
  }
  else {
    dataSetLong <- dataLong (dataSet=dataSet, timeColumn=as.character(survModelFormula) [2], censColumn=censColumn)
  }
  
  # Fit generalized, linear model
  NewFormula <- update(survModelFormula, y ~ timeInt + .)
  glmFit <- glm(formula=NewFormula, data=dataSetLong, family=binomial(link=linkFunc), control=glm.control(maxit=2500))
  hazards <- predict(glmFit, type="response")
  
  # Calculate residuals
  splitHazards <- split(hazards, dataSetLong$obj)
  splitY <- split(dataSetLong$y, dataSetLong$obj)
  Residuals <- sapply(1:length(splitY), SquDevResid)
  Output <- list(Output=list(DevResid=sqrt(Residuals), GlmFit=glmFit), 
                 Input=list(dataSet=dataSet, survModelFormula=survModelFormula, censColumn=censColumn, linkFunc=linkFunc, idColumn=idColumn))
  class(Output) <- "discSurvDevResid"
  return(Output)
}

print.discSurvDevResid <- function (x, ...) {
  print(round(x$Output$DevResid, 4))
}

########################
# adjDevResid

# Description
# Calculates the adjusted deviance residuals. Should be normal distributed, in the case of a well fitting model

adjDevResid <- function (dataSet, survModelFormula, censColumn, linkFunc="logit", idColumn=NULL) {
  # Input checks
  if(!is.data.frame(dataSet)) {stop("Argument *dataSet* is not in the correct format! Please specify as data.frame object.")}
  if(!("formula" %in% class(survModelFormula))) {stop("*survModelFormula* is not of class formula! Please specify a valid formula, e. g. y ~ x + z.")}
  if(!any(names(dataSet)==censColumn)) {stop("Argument *censColumn* is not available in *dataSet*! Please specify the correct column name of the event indicator.")}
  if(!(any(names(dataSet)==idColumn) | is.null(idColumn))) {stop("Argument *idColumn* is not available in *dataSet*! Please specify the correct column name of the identification numbers of persons.")}
  
  # Help function
  AdjDevResid <- function (x) { 
    LogTerm1 <- ifelse(splitY [[x]]==1, -log(splitHazards [[x]]), 0)
    LogTerm2 <- ifelse(splitY [[x]]==0, -log (1 - splitHazards [[x]]), 0)
    FirstPartialSum <- sum(sign(splitY [[x]] - splitHazards [[x]]) * (sqrt(splitY [[x]] * LogTerm1 + (1 - splitY [[x]]) * LogTerm2)))
    SecondPartialSum <- sum( (1 - 2*splitHazards [[x]]) / sqrt (splitHazards [[x]] * (1 - splitHazards [[x]]) * 36) )
    return(FirstPartialSum + SecondPartialSum)
  }
  
  # Convert to long format
  if(!is.null(idColumn)) {
    dataSetLong <- dataLongTimeDep (dataSet=dataSet, timeColumn=as.character(survModelFormula) [2], censColumn=censColumn, idColumn=idColumn)
  }
  else {
    dataSetLong <- dataLong (dataSet=dataSet, timeColumn=as.character(survModelFormula) [2], censColumn=censColumn)
  }
  
  # Fit generalized, linear model
  NewFormula <- update(survModelFormula, y ~ timeInt + .)
  glmFit <- glm(formula=NewFormula, data=dataSetLong, family=binomial(link=linkFunc))
  hazards <- predict(glmFit, type="response")
  
  # Calculate residuals
  splitHazards <- split(hazards, dataSetLong$obj)
  splitY <- split(dataSetLong$y, dataSetLong$obj)
  Residuals <- sapply(1:length(splitY), AdjDevResid)
  Output <- list(Output=list(AdjDevResid=Residuals, GlmFit=glmFit), 
                 Input=list(dataSet=dataSet, survModelFormula=survModelFormula, censColumn=censColumn, linkFunc=linkFunc, idColumn=idColumn))
  class(Output) <- "discSurvAdjDevResid"
  return(Output)
}

print.discSurvAdjDevResid <- function (x, ...) {
  print(round(x$Output$AdjDevResid, 4))
}

plot.discSurvAdjDevResid <- function (x, ...) {
  qqnorm (y=x$Output$AdjDevResid, las=1, ...)
  qqline(y=x$Output$AdjDevResid, ...)
}

########################
# adjDevResidShort

# Description
# Calculates the adjusted deviance residuals for arbitrary hazard prediction models. Should be normal distributed, in the case of a well fitting model

adjDevResidShort <- function (dataSet, hazards) {
  # Input checks
  if(!is.data.frame(dataSet)) {stop("Argument *dataSet* is not in the correct format! Please specify as data.frame object.")}
  if(!all(hazards>=0 & hazards<=1)) {stop("Argument *hazards* must contain probabilities in the closed interval from zero to one. Please verify that *hazards* are estimated hazard rates")}
  if(!(dim(dataSet)[1]==length(hazards))) {stop("The length of argument *hazards* must match the number of observations")}
  
  # Help function
  AdjDevResid <- function (x) { 
    LogTerm1 <- ifelse(splitY [[x]]==1, -log(splitHazards [[x]]), 0)
    LogTerm2 <- ifelse(splitY [[x]]==0, -log (1 - splitHazards [[x]]), 0)
    FirstPartialSum <- sum(sign(splitY [[x]] - splitHazards [[x]]) * (sqrt(splitY [[x]] * LogTerm1 + (1 - splitY [[x]]) * LogTerm2)))
    SecondPartialSum <- sum( (1 - 2*splitHazards [[x]]) / sqrt (splitHazards [[x]] * (1 - splitHazards [[x]]) * 36) )
    return(FirstPartialSum + SecondPartialSum)
  }

  # Calculate residuals
  splitHazards <- split(hazards, dataSet$obj)
  splitY <- split(dataSet$y, dataSet$obj)
  Residuals <- sapply(1:length(splitY), AdjDevResid)
  Output <- list(Output=list(AdjDevResid=Residuals), 
                 Input=list(dataSet=dataSet, hazards=hazards))
  class(Output) <- "discSurvAdjDevResid"
  return(Output)
}

########################
# devResidShort

# Description
# Computes the root of the squared deviance residual

devResidShort <- function (dataSet, hazards) {
  
  # Input checks
  if(!is.data.frame(dataSet)) {stop("Argument *dataSet* is not in the correct format! Please specify as data.frame object.")}
  if(!all(hazards>=0 & hazards<=1)) {stop("Argument *hazards* must contain probabilities in the closed interval from zero to one. Please verify that *hazards* are estimated hazard rates")}
  if(!(dim(dataSet)[1]==length(hazards))) {stop("The length of argument *hazards* must match the number of observations")}
  
  # Help function
  SquDevResid <- function (x) {-2*sum(splitY [[x]] * log(splitHazards [[x]]) + (1 - splitY [[x]]) * log(1 - splitHazards [[x]] ))}

  # Calculate residuals
  splitHazards <- split(hazards, dataSet$obj)
  splitY <- split(dataSet$y, dataSet$obj)
  Residuals <- sapply(1:length(splitY), SquDevResid)
  Output <- list(Output=list(DevResid=sqrt(Residuals)), 
                 Input=list(dataSet=dataSet, hazards=hazards))
  class(Output) <- "discSurvDevResid"
  return(Output)
}
