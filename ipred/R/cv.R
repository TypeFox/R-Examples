#$Id: cv.R,v 1.21 2004/02/11 09:13:51 peters Exp $

cv <- function(y, ...) {
  if(is.null(class(y)))
    class(y) <- data.class(y)
  UseMethod("cv", y)
}

cv.default <- function(y, ...) {
  stop(paste("Do not know how to handle objects of class", class(y)))
}

cv.integer <- function(y, ...) {
  cv.numeric(y, ...)
}

cv.factor <- function(y, formula, data, model, predict, k=10, random=TRUE, 
                      strat=FALSE, predictions=NULL, getmodels=NULL, list.tindx = NULL, ...) {

  # k-fold cross-validation of misclassification error

  if (!is.data.frame(data)) stop("data is not of class data.frame")
 
  N <- length(y)
  classes <- levels(y)

  if (is.null(k)) k <- 10
  if (is.null(random)) random <- TRUE
  if (is.null(strat)) strat <- FALSE
  if (is.null(predictions)) predictions <- FALSE
  if (is.null(getmodels)) getmodels <- FALSE
  USEPM <- FALSE

  if(!is.null(list.tindx)) k <- length(list.tindx)
  if(!is.null(list.tindx)) {
    random <- FALSE
  }
  
  # to reproduce results, either use `set.seed' or a fixed partition of 
  # the samples
  if (random) 
    myindx <- sample(1:N, N)
  else 
    myindx <- 1:N

  y <- y[myindx]
  data <- data[myindx,]
  
    # determine an appropriate splitting for the sample size into
    # k roughly equally sized parts
  
  mysplit <- ssubset(y, k, strat=strat)
  allpred <- vector(mode="character", length=N)
  fu <- function(x) levels(x)[as.integer(x)]
  nindx <- 1:N

  if (getmodels)
    models <- vector(k, mode="list")

  for(i in 1:k) {

    if(!is.null(list.tindx)) {
      tindx <- list.tindx[[i]]
    } else {
      tindx <- mysplit[[i]]
    }
    
    folddata <- subset(data, !(nindx %in% tindx))
    mymodel <- model(formula, data=folddata, ...)
    if (getmodels) models[[i]] <- mymodel

    # check of mymodel is a function which should be used instead of
    # predict
    if (is.function(mymodel)) {
      if(!is.null(predict) & i == 1) 
        warning("model returns a function and predict is specified, using models output")
      predict <- mymodel
      USEPM <- TRUE
    }

    # we assume predict to return factor levels
    if (USEPM)
      pred <- predict(newdata=data)
    else 
      pred <- predict(mymodel, newdata = data)
    if (!is.factor(pred)) stop("predict does not return factor values")
    pred <- factor(pred, levels=classes)
    
    # <FIXME>
    # there is no c() for factors which preserves the levels, isn't it?
    # use characters
    allpred[tindx] <- fu(pred[tindx])
    # </FIXME>
  }
  allpred <- factor(allpred, levels=classes)
  allpred <- allpred[order(myindx)]
  err <- mean(allpred != y[order(myindx)], na.rm = TRUE)
  if (predictions)
    RET <- list(error = err, k = k, predictions=allpred)
  else 
    RET <- list(error = err, k = k)
  if (getmodels)
    RET <- c(RET, models=list(models))
  class(RET) <- "cvclass"
  RET
}

cv.numeric <- function(y, formula, data, model, predict, k=10, random=TRUE,
                       predictions=NULL, strat=NULL, getmodels=NULL, list.tindx = NULL, ...) {

  # k-fold cross-validation of mean squared error 

  if (!is.data.frame(data)) stop("data is not of class data.frame")
  if(!is.null(list.tindx)) k <- length(list.tindx)
  N <- length(y)

  if (is.null(k)) k <- 10
  if (is.null(random)) random <- TRUE
  if (is.null(predictions)) predictions <- FALSE
  if (is.null(getmodels)) getmodels <- FALSE   
  USEPM <- FALSE
  # determine an appropriate splitting for the sample size into
  # k roughly equally sized parts

#  if(is.null(list.tindx)) {
    a <- kfoldcv(k, N)
    # to reproduce results, either use `set.seed' or a fixed partition of
    # the samples
    if (random)
      myindx <- sample(1:N, N)
    else
      myindx <- 1:N
    nindx <- 1:N
#  }
  
  if (getmodels)
    models <- vector(k, mode="list")

  allpred <- rep(0, N)
  for(i in 1:k) {
    if(!is.null(list.tindx)) {
      tindx <- list.tindx[[i]]
    } else {
      if (i > 1)
        tindx <- myindx[(sum(a[1:(i-1)])+1):sum(a[1:i])]
      else
        tindx <- myindx[1:a[1]]
    }
    folddata <- subset(data, !(nindx %in% tindx))
    mymodel <- model(formula, data=folddata, ...)
    if (getmodels) models[[i]] <- mymodel

    # check of mymodel is a function which should be used instead of
    # predict
    if (is.function(mymodel)) {   
      if(!is.null(predict) & i == 1) 
        warning("model returns a function and predict is specified, using models output")
      predict <- mymodel  
      USEPM <- TRUE  
    }  

    outfolddata <- subset(data, nindx %in% tindx)
    if (USEPM)
      pred <- predict(newdata=outfolddata)
    else
      pred <- predict(mymodel, newdata = outfolddata)
    if (!is.numeric(pred)) stop("predict does not return numerical values")
    allpred[sort(tindx)] <- pred
  }
  err <- sqrt(mean((allpred - y)^2, na.rm = TRUE))
  if (predictions)
    RET <- list(error = err, k = k, predictions=allpred)
  else
    RET <- list(error = err, k = k)
  if (getmodels) 
    RET <- c(RET, models=list(models))
  class(RET) <- "cvreg" 
  RET  
}

cv.Surv <- function(y, formula, data=NULL, model, predict, k=10, random=TRUE,
                    predictions=FALSE, strat=FALSE, getmodels=NULL, list.tindx = NULL, ...) {

  # k-fold cross-validation of Brier's score

  if (is.null(predictions)) predictions <- FALSE
  if(is.null(random)) random <- TRUE
  if (is.null(predictions)) predictions <- FALSE
  if (is.null(strat)) strat <- FALSE
  if (is.null(getmodels)) getmodels <- FALSE   
  USEPM <- FALSE

  N <- length(y[,1])
  nindx <- 1:N
  if(is.null(random)) random <- TRUE
  if(is.null(k)) k <- 10
  if (is.null(data)) data <- rep(1, N)

  if(!is.null(list.tindx)) k <- length(list.tindx)
  if(is.null(k)) stop("k for k-fold cross-validation is missing")
 
  # determine an appropriate splitting for the sample size into
  # k roughly equally sized parts

 # if(is.null(list.tindx)) {
    a <- kfoldcv(k, N)
    # to reproduce results, either use `set.seed' or a fixed partition of
    # the samples
    if (random)
      myindx <- sample(1:N, N)
    else
      myindx <- 1:N
 # }

  if (getmodels)
    models <- vector(k, mode="list")

  cverr <- c()
  for(i in 1:k) {
    if(!is.null(list.tindx)) {
      tindx <- list.tindx[[i]]
    } else {
      if (i > 1)
        tindx <- myindx[(sum(a[1:(i-1)])+1):sum(a[1:i])]
      else
        tindx <- myindx[1:a[1]]
    }

    folddata <- subset(data, !(nindx %in% tindx))
    mymodel <- model(formula, data=folddata, ...)
    if (getmodels) models[[i]] <- mymodel

    # check if mymodel is a function which should be used instead of
    # predict
    if (is.function(mymodel)) {   
      if(!is.null(predict) & i == 1) 
        warning("model returns a function and predict is specified, using models output")
      predict <- mymodel  
      USEPM <- TRUE  
    }  

    outfolddata <- subset(data, (nindx %in% tindx))
    if (USEPM)
      pred <- predict(newdata=outfolddata)
    else
      pred <- predict(mymodel, newdata = outfolddata)
    if (is.list(pred)) {
      if (!inherits(pred[[1]], "survfit") && !inherits(pred, "survfit"))
        stop("predict does not return a list of survfit objects")
    } else {
      stop("predict does not return a list of survfit objects")
    }

    err <- sbrier(y[sort(tindx)], pred)
    cverr <- c(cverr,rep(err, length(tindx)))
  }
  RET <- list(error = mean(cverr), k=k)
  if (getmodels) 
    RET <- c(RET, models=list(models))
  class(RET) <- "cvsurv" 
  RET  
}
       
