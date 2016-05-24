#' @include 01BaseClass.R
NULL


# DEFAULT SUBCLASSES DEFINITION AND ASSOCIATED METHODS ====================================

constructfeature("MissingValueShare", "apply(data, 1, function(x) sum(is.na(x))/ncol(data))", impute=TRUE)

aux0 <- function(data){
  nacount <- apply(data, 1, function(x) sum(is.na(x)))
  temp1 <- data[sapply(data,is.factor)][,1]
  temp2 <- data.frame(nacount, temp1)
  temp3 <- tapply(temp2[,1], temp2[,2], mean)
  levels(temp1) <- as.numeric(temp3)
  temp4 <- as.numeric(as.character(temp1))
  result <- nacount/ temp4
}

constructfeature("MissingValueToClass", "aux0(data)", impute=TRUE)

constructfeature("LOFScore", "DMwR::lofactor(data, k=5)", mode="numeric")

constructfeature("DistanceToNearest", "apply(as.matrix(dist(data)), 1, function(x) min(x[x>0]))", mode="numeric")

constructfeature("LenghtOfIQR", "apply(data, 1, function(x) IQR( (x-min(x))/(max(x)-min(x))) )", "numeric")

aux1 <- function(data) 
{
  temp <- data[sapply(data,is.numeric)]
  temp1 <- data[sapply(data,is.factor)]
  nearest <- apply(as.matrix(dist(temp)), 1, function(x) which.min(x[x>0]))
  sequence <- seq(1, nrow(temp), 1)
  result <- as.integer(temp1[nearest,1]==temp1[sequence,1])
}

constructfeature("NearestPointPurity", "aux1(data)")

aux2 <- function(data) 
{
  
  temp <- data[sapply(data,is.numeric)]
  temp1 <- data[sapply(data,is.factor)]
  area <- round(nrow(temp)/10,0)
  nearest <- apply(as.matrix(dist(temp)), 1, function(x) order(x))
  nearest1 <- data.frame(nearest[, c(1:area)])
  result <- as.numeric(apply(nearest1, 1, function(x) max(table(temp1[x,1]))/area))
}

constructfeature("NeighborhoodDiversity", "aux2(data)")

aux3 <- function(data){

  temp1 <- data[sapply(data,is.factor)][,1]
  model.rf <- randomForest::randomForest(temp1 ~ ., data=data, importance=TRUE, proximity=TRUE)
  predict.rf <- predict(model.rf, type="prob")
  result <- apply(predict.rf, 1, function(x) max(x)/1)
}

constructfeature("ClassificationCertainty", "aux3(data)")

aux4 <- function(data){
  tryCatch({
  temp <- data[sapply(data,is.numeric)]
  temp1 <- data[sapply(data,is.factor)][,1]
  temp2 <- as.character(temp1)
  distance <- caret::classDist(temp, temp1, pca=TRUE)
  result <- data.frame(predict(distance, temp), temp2, stringsAsFactors=FALSE)
  temp3 <- gsub("dist.", "", colnames(result))
  temp4 <- cbind(a=seq(1,nrow(temp),1), b=match(temp2, temp3))
  result1 <- abs(as.numeric(result[temp4]))
  }, error= function(e) return(rep(nrow(data),1)) )
}

constructfeature("MahalanobisToClassCenter", "aux4(data)")



aux6 <- function(data){
  result <- numeric(nrow(data))
  area <- round(nrow(data)/10,1)
  temp <- data[sapply(data, is.numeric)]
  temp1 <- data[sapply(data, is.factor)][,1]
  distance <- as.data.frame(as.matrix(dist(temp, diag=TRUE, upper=TRUE)))
  sequence <- data.frame(apply(distance, 1, function(x) order(x)[1:area])) 

  aux_temp <- function(x)
  {
    counter <- 0
    for (i in 1:(area-1))
    {
      if (temp1[x[i]]!=temp1[x[i+1]]) {counter <- counter+1}
    }
    return(counter)
  }
  result <- apply(sequence, 2, function(x) aux_temp(x))

}

constructfeature("ScatterCounter", "aux6(data)")


#library(doSNOW)
#registerDoSNOW(makeCluster(2, type = "SOCK"))

#aux5 <- function(data){
#  result <- numeric(nrow(data))
#  temp <- data[sapply(data, is.numeric)]
#  base <- clusterttend::hopkins(temp, n=(nrow(data)-1))$H
#  result <- foreach(i = 1:nrow(data), .combine = "c", .export=c("hopkins")) %dopar%{
#    comp <- (hopkins(temp[-i,], n=nrow(data)-2)$H)/base
#    #result[i] <- comp/base
#  }
#  return(result)
#}

#constructfeature("ClusteringTendency", "aux5(data)")