CreateTrainTestSet <- function(data, fractionTrain=0.5, type="regression"){
  selecttrain <- sample(1:length(data$y), size=floor(fractionTrain*length(data$y)))
  if(type=="survival"){
    train <- list(x=data$x[,selecttrain], y=data$y[selecttrain],
      censoring.status=data$censoring.status[selecttrain])
    test <- list(x=data$x[,-selecttrain], y=data$y[-selecttrain],
      censoring.status=data$censoring.status[-selecttrain])
  } else if(type=="regression"){
    train <- list(x=data$x[,selecttrain], y=data$y[selecttrain])
    test <- list(x=data$x[,-selecttrain], y=data$y[-selecttrain])
  } else if (type=="two class"){
    whichers <- which(data$y==1)
    testindices <- sample(whichers, floor(length(whichers)*(1-fractionTrain)))
    whichers <- which(data$y==2)
    testindices <- c(testindices, sample(whichers, floor(length(whichers)*(1-fractionTrain))))
    train <- list(x=data$x[,-testindices], y=data$y[-testindices])
    test <- list(x=data$x[,testindices], y=data$y[testindices])
  } else if (type=="multiclass"){
    testindices <- NULL
    for(i in unique(data$y)){
      whichers <- which(data$y==i)
      testindices <- c(testindices, sample(whichers, floor(length(whichers))*(1-fractionTrain)))
    }
    train <- list(x=data$x[,-testindices], y=data$y[-testindices])
    test <- list(x=data$x[,testindices], y=data$y[testindices])
  } else stop("Not an acceptable type.")
  return(list(train=train, test=test))
}
