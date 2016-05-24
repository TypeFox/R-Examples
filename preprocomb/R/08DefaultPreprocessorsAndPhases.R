#' @include 07AnalysisComponent.R
NULL

# DEFAULT PREPROCESSORS ==========================

# NEAR ZERO VARIANCE

nezevar <- function(dataobject){
  temp <- caret::nearZeroVar(dataobject@x)
  if (length(temp) !=0) {
    x_new <- dataobject@x[,-temp]
    dataobject <- initializedataclassobject(data.frame(x_new, dataobject@y))
}
  return(dataobject)
}

setpreprocessor("nearzerovar", "nezevar(dataobject)")

## IMPUTATION

naomit <- function(dataobject){
  temp <- apply(dataobject@x, 1, function(x) any(is.na(x)))
  if (any(temp==TRUE)){
  x_new <- dataobject@x[-temp==FALSE,]
  y_new <- dataobject@y[-temp==FALSE]
  dataobject <- initializedataclassobject(data.frame(x_new, y_new))
}
  return(dataobject)
}

setpreprocessor("naomit", "naomit(dataobject)")

# Mean imputation

meanimpute_aux <- function(x){
  x[is.na(x)] <- mean(x, na.rm=TRUE) #convert the item with NA to median value from the column
  x
}

meanimpute <- function(dataobject){
  temp <- apply(dataobject@x, 1, function(x) any(is.na(x)))
  if (any(temp==TRUE)){
  x_new <- data.frame(apply(dataobject@x, 2, function(x) meanimpute_aux(x)))
  dataobject <- initializedataclassobject(data.frame(x_new, dataobject@y))
}
  return(dataobject)
}

setpreprocessor("meanimpute", "meanimpute(dataobject)")

# mean class imputation

meanclass <- function(dataobject){
  temp <- apply(dataobject@x, 1, function(x) any(is.na(x)))
  if (any(temp==TRUE)){
    x_new <- zoo::na.aggregate(dataobject@x, dataobject@y, FUN=mean)
    dataobject <- initializedataclassobject(data.frame(x_new, dataobject@y))
  }
  return(dataobject)
}

setpreprocessor("meanclassimpute", "meanclass(dataobject)")

# knnimputation

knnimpute <- function(dataobject){
  temp <- apply(dataobject@x, 1, function(x) any(is.na(x)))
  if (any(temp==TRUE)){
  x_new <- DMwR::knnImputation(dataobject@x, k=5)
  dataobject <- initializedataclassobject(data.frame(x_new, dataobject@y))
  }
  return(dataobject)
}

setpreprocessor("knnimpute", "knnimpute(dataobject)")

# random forest imputation

rfimputefunc <- function(dataobject){
  temp <- apply(dataobject@x, 1, function(x) any(is.na(x)))
  if (any(temp==TRUE)){
  res <- randomForest::rfImpute(dataobject@y ~ ., dataobject@x)
  x_new <- res[,2:ncol(res)]
  y_new <- res[,1]
  dataobject <- initializedataclassobject(data.frame(x_new, y_new))
  }
  return(dataobject)
}

setpreprocessor("randomforestimpute", "rfimputefunc(dataobject)")

## SCALING

# basic

basicscaling <- function(dataobject){
  dataobject <- initializedataclassobject(data.frame(x=scale(dataobject@x, center=FALSE), dataobject@y))
}

setpreprocessor("basicscale", "basicscaling(dataobject)")

# center

centerscaling <- function(dataobject){
  dataobject <- initializedataclassobject(data.frame(x=scale(dataobject@x, center=TRUE), dataobject@y))
}

setpreprocessor("centerscale", "centerscaling(dataobject)")

# min-max scaling

range01 <- function(y){
  newrange <- (y-min(y))/(max(y)-min(y))
}

minmaxscaling <- function(dataobject){

  x_new <- data.frame(apply(dataobject@x, 2, range01))
  dataobject <- initializedataclassobject(data.frame(x=x_new, dataobject@y))
}

setpreprocessor("minmaxscale", "minmaxscaling(dataobject)")

decimalscaling <- function(dataobject){
  m <- apply(dataobject@x,2,max)
  v <- round(log10(m),0)
  c <- 10^v
  x_new <- sweep(dataobject@x, 2, c, `/`)
dataobject <- initializedataclassobject(data.frame(x=x_new, y=dataobject@y))
}

setpreprocessor("decimalscale", "decimalscaling(dataobject)")

softmaxscaling <- function(dataobject){
  dataobject <- initializedataclassobject(data.frame(x=data.frame(apply(dataobject@x, 2, DMwR::SoftMax)), dataobject@y))
}

setpreprocessor("softmaxscale", "softmaxscaling(dataobject)")


# OUTLIER REMOVAL

orhcut <- function(dataobject){
  x_original <- dataobject@x
  y_original <- dataobject@y
  orh_score <- suppressMessages(DMwR::outliers.ranking(x_original))
  orh_rank <- orh_score$prob.outliers[orh_score$rank.outliers]
  orh_cut <- quantile(orh_rank, .95)
  orh_obs <- as.integer(names(which(orh_rank >= orh_cut)))
  x_preprocessed <- x_original[-orh_obs,]
  y_preprocessed <- y_original[-orh_obs]
  result <- initializedataclassobject(data.frame(x_preprocessed, y=y_preprocessed))
}

setpreprocessor("orhoutlier", "orhcut(dataobject)")

## CLASS IMBALANCE CORRECTIONS

# oversampling

oversample <- function(dataobject){

  data <- data.frame(dataobject@x, y=dataobject@y)

  if (nlevels(data$y) > 2) {stop("Oversampling can only be applied to binary class.")}

  freq <- table(data$y)
  temp <- order(freq)
  nsample <- freq[temp[2]]-freq[temp[1]]
  indexes <- which(data$y==names(freq)[temp[1]])
  tempsample <- sample(indexes, nsample, replace=TRUE)
  newdata <- initializedataclassobject(data.frame(rbind(data, data[tempsample,])))
  return(newdata)
}

setpreprocessor("oversample", "oversample(dataobject)")

# undersampling

undersample <- function(dataobject){

  data <- data.frame(dataobject@x, y=dataobject@y)

  if (nlevels(data$y) > 2) {stop("Undersampling can only be applied to binary class.")}

  freq <- table(data$y)
  temp <- order(freq)
  nsample <- freq[temp[1]]
  indexes <- which(data$y==names(freq)[temp[2]])
  indexes2 <- which(data$y==names(freq)[temp[1]])
  tempsample <- sample(indexes, nsample)
  finalsample <- c(indexes2, tempsample)
  newdata <- initializedataclassobject(data[finalsample,])
  return(newdata)
}

setpreprocessor("undersample", "undersample(dataobject)")

# smote sampling

#smotesample <- function(dataobject){
#  temp <- data.frame(dataobject@x, y=dataobject@y)
#
#  if (nlevels(temp$y) > 2) {stop("SMOTE can only be applied to binary class.")}
#
#  temp1 <- as.numeric(table(temp$y))
#  temp2 <- order(temp1)
#  temp3 <- temp1[temp2[2]]/temp1[temp2[1]]
#  temp4 <- temp1[temp2[2]]/(temp1[temp2[2]]-temp1[temp2[1]])
#  newData <- DMwR::SMOTE(y ~ ., temp, perc.over = 100*temp3, perc.under=100*temp4)
#  dataobject <- initializedataclassobject(newData)
#}

setpreprocessor("smotesample", "smotesample(dataobject)")

## FEATURE SELECTION

# random forest importance

rfimportance <- function(dataobject, qt){
  rf.imp <- randomForest::randomForest(dataobject@y ~ ., data=dataobject@x, ntree=100, keep.forest=FALSE, importance=TRUE)
  temp <- data.frame(randomForest::importance(rf.imp))
  temp1 <- temp$MeanDecreaseAccuracy > quantile(temp$MeanDecreaseAccuracy, qt)
  dataobject@x  <- dataobject@x[,temp1]
  return(dataobject)
}

setpreprocessor("rfselect75", "rfimportance(dataobject, .25)")
setpreprocessor("rfselect50", "rfimportance(dataobject, .50)")

# smoothing with lowess

smoothlowess <- function(dataobject){
  y <- dataobject@x
  result <-data.frame(round(apply(y, 2, function(y) lowess(y[order(y)], f=1/2)[[2]][match(y, y[order(y)])]),2))
  dataobject <- initializedataclassobject(data.frame(result, y=dataobject@y))
  }

setpreprocessor("lowesssmooth", "smoothlowess(dataobject)")

smoothcoarse <- function(dataobject){
   x_new <- data.frame(apply(dataobject@x, 2, function(y) round(y, digits=-log10(abs(y))+1)))
   dataobject <- initializedataclassobject(data.frame(x=x_new, y=dataobject@y))
}

setpreprocessor("coarsesmooth", "smoothcoarse(dataobject)")

setpreprocessor("noaction", "identity(dataobject)")

# DEFAULT PHASES ==========================


#### PHASES

imputation <- setphase("imputation", c("naomit", "meanimpute", "meanclassimpute", "knnimpute", "randomforestimpute"), TRUE)
variance <- setphase("variance", c("noaction", "nearzerovar"), FALSE)
smoothing <- setphase("smoothing", c("noaction", "coarsesmooth", "lowesssmooth"), FALSE)
scaling <- setphase("scaling", c("noaction", "basicscale", "centerscale", "minmaxscale", "decimalscale", "softmaxscale"), FALSE)
outliers <- setphase("outliers", c("noaction", "orhoutlier"), FALSE)
sampling <- setphase("imbalance", c("noaction", "oversample", "undersample"), FALSE) # add:  "smotesample"
selection <- setphase("selection", c("noaction", "rfselect50", "rfselect75"), FALSE)

#' preprodefault
#'
#' Default phases
#' @details preprodefault object can be used as default phases for setgrid
#' @export
preprodefault <- c("imputation", "variance", "smoothing", "scaling", "outliers", "sampling", "selection")
