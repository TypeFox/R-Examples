# $Id: bootest.R,v 1.18 2004/02/09 08:08:21 peters Exp $

bootest <- function(y, ...) {
  if(is.null(class(y)))
    class(y) <- data.class(y)
  UseMethod("bootest", y)
}

bootest.default <- function(y, ...) {
  stop(paste("Do not know how to handle objects of class", class(y)))
}

bootest.integer <- function(y, ...) {
  bootest.numeric(y, ...)
}

bootest.factor <- function(y, formula, data, model, predict, 
                           nboot=25, bc632plus = FALSE, list.tindx = NULL, predictions = FALSE, both.boot = FALSE, ...) {

  # bootstrap estimator of misclassification error
  
  N <- length(y)
  nindx <- 1:N
  if(!is.null(list.tindx)) nboot <- length(list.tindx)
  
  bootindx <- matrix(NA, ncol=nboot, nrow=N)
  if(predictions) {
    BOOTINDX <- data.frame(matrix(NA, ncol=nboot, nrow=N))
  }
  
  classes <- levels(y)
  USEPM <- FALSE
  
  if(!is.data.frame(data)) stop("data is not a data.frame")
  if(nboot <=2) stop("to small number of bootstrap replications")
  if(is.null(nboot)) stop("number of bootstrap replications is missing")
  if(!is.null(list.tindx) & length(list.tindx) != nboot) stop(paste("List of selected observations per bootstrap sample has to be", nboot))
  
  for(i in 1:nboot) {
    if(!is.null(list.tindx)) {
      tindx <- list.tindx[[i]]
      if(length(tindx) > N) warning("number of selected observations is larger than the sample size")
    } else {
      tindx <- sample(nindx, N, replace = TRUE)
    }
    
    mymodel <- model(formula, data = data[tindx,], ...)

    # check if mymodel is a function which should be used instead of   
    # predict
    if (is.function(mymodel)) {
      if(!is.null(predict) & i == 1) 
        warning("model returns a function and predict is specified, using models output")
      predict <- mymodel
      USEPM <- TRUE
    }

    if (USEPM) 
      pred <- predict(newdata=data)
    else 
      pred <- predict(mymodel, newdata = data)
    if (!is.factor(pred)) stop("predict does not return factor values")
    pred <- factor(pred, levels=classes)[-tindx]
    if (length(pred) != length(y[-tindx]))
        stop("different length of data and prediction")
    if(predictions) {
      BOOTINDX[,i] <- factor(BOOTINDX[,i],levels = classes)
      BOOTINDX[-tindx, i] <- pred
    }
    bootindx[-tindx, i] <- (pred != y[-tindx])
  }

  fun <- function(x)
       ifelse(all(is.na(x)), NA, mean(as.integer(x), na.rm = TRUE))

  one <- mean(apply(bootindx, 1, fun), na.rm = TRUE)

  if (bc632plus) {
    full.model <- model(formula, data = data, ...)
    # check if full.model is a function which should be used instead of
    # predict
    if (is.function(full.model)) {
      predict <- full.model
      USEPM <- TRUE
    }

    if (USEPM)
      full.pred <- predict(newdata=data)
    else

    full.pred <- predict(full.model, newdata = data)
    resubst <- mean(full.pred != y, na.rm = TRUE)

    err632 <- 0.368*resubst + 0.632*one

    y <- y[!is.na(y) & !is.na(full.pred)]
    full.pred <- full.pred[!is.na(y) & !is.na(full.pred)]
    gamma <- sum(outer(y, full.pred, function(x, y) ifelse(x==y, 0, 1) ))/
                 (length(y)^2)
    r <- (one - resubst)/(gamma - resubst)
    r <- ifelse(one > resubst & gamma > resubst, r, 0)
    errprime <- min(one, gamma)
    #    weight <- .632/(1-.368*r)
    #    err <- (1-weight)*resubst + weight*one
    err <- err632 + (errprime - resubst)*(0.368*0.632*r)/(1-0.368*r)
    if(predictions) 
      RET <- list(error = err, nboot = nboot, bc632plus = TRUE, predictions = BOOTINDX)
    else
      RET <- list(error = err, nboot=nboot, bc632plus = TRUE)
    if(both.boot){
      bc632plus <- FALSE
      RETbc <- RET
    }
  }

  if(!bc632plus) {
    err <- one
    expb <- rep(0, nboot)
    for(i in 1:nboot)
      expb[i] <- mean(apply(bootindx[,-i], 1, fun), na.rm = TRUE)
    
    sdint <- sqrt( ((nboot - 1)/nboot)*sum((expb - mean(expb))^2) )
    if(predictions) 
      RET <- list(error = err, sd = sdint, bc632plus = FALSE, nboot = nboot, predictions = BOOTINDX)
    else
      RET <- list(error = err, sd=sdint, bc632plus=FALSE, nboot=nboot)
    if(both.boot){
      RET <- list("boot" = RET, "632plus" = RETbc)
    }
  }
  class(RET) <- "bootestclass"
  RET
}

bootest.numeric <- function(y, formula, data, model, predict, 
                           nboot=25, bc632plus=FALSE, list.tindx = NULL, predictions = FALSE, ...) {
  
  # bootstrap estimator of root of mean squared error 

  if (bc632plus) stop("cannot compute 632+ estimator of mean squared error")
  if(!is.null(list.tindx)) nboot <- length(list.tindx)
  if (nboot <=2) stop("to small number of bootstrap replications")

  ##FIX: nrow = 
  N <- length(y)
  nindx <- 1:N
  
  bootindx <- matrix(NA, ncol=nboot, nrow=N)
  if(predictions) BOOTINDX <- matrix(NA, ncol=nboot, nrow=N)
  USEPM <- FALSE
  
  if (!is.data.frame(data)) stop("data is not a data.frame")
  if(is.null(nboot)) stop("number of bootstrap replications is missing")
  if(!is.null(list.tindx) & length(list.tindx) != nboot) stop(paste("List of selected observations per bootstrap sample has to be", nboot))
  

  for(i in 1:nboot) {
    if(!is.null(list.tindx)) {
      tindx <- list.tindx[[i]]
      if(length(tindx) > N) warning("number of selected observations is larger than the sample size")
    } else {
      tindx <- sample(nindx, N, replace = TRUE)
    }
#    tindx <- ifelse(!is.null(list.tindx), list.tindx[[i]], sample(nindx, N, replace = TRUE))
    mymodel <- model(formula, data = data[tindx,], ...)
    outbootdata <- subset(data, !(nindx %in% tindx))
    # check if mymodel is a function which should be used instead of
    # predict
    if (is.function(mymodel)) {
      if(!is.null(predict) & i == 1) 
        warning("model returns a function and predict is specified, using models output")
      predict <- mymodel
      USEPM <- TRUE
    }

    if (USEPM)
      pred <- predict(newdata=outbootdata)
    else
      pred <- predict(mymodel, newdata = outbootdata)
    if (!is.numeric(pred)) stop("predict does not return numerical values")
    if (length(pred) != length(y[-tindx]))
        stop("different length of data and prediction")
    if(predictions) BOOTINDX[-tindx, i] <- pred 
    bootindx[-tindx, i] <- (pred - y[-tindx])^2
  }

  fun <- function(x)
        ifelse(all(is.na(x)), NA, mean(x, na.rm = TRUE))

  err <- sqrt(mean(apply(bootindx, 1, fun), na.rm = TRUE))
  if(predictions) 
    RET <- list(error = err, nboot = nboot, predictions = BOOTINDX)
  else
    RET <- list(error = err, nboot=nboot)
  class(RET) <- "bootestreg"
  RET
}

bootest.Surv <- function(y, formula, data=NULL, model, predict, 
                           nboot=25, bc632plus=FALSE, list.tindx = NULL, predictions = FALSE, ...) {
  
  # bootstrap estimator of Brier's score

  if (bc632plus) stop("cannot compute 632+ estimator of Brier's score")
 
  N <- dim(y)[1]
  if(!is.null(list.tindx)) nboot <- length(list.tindx)
  nindx <- 1:N

  bootindx <- matrix(NA, ncol=nboot, nrow=N)
  if(predictions) BOOTINDX <- matrix(NA, ncol=nboot, nrow=N)
  USEPM <- FALSE

  if(is.null(nboot)) stop("number of bootstrap replications is missing")
  if (nboot <=2) stop("to small number of bootstrap replications")
  if (is.null(data)) data <- as.data.frame(rep(1, N))
  if (!is.data.frame(data)) stop("data is not a data.frame")
  if(!is.null(list.tindx)) nboot <- length(list.tindx)
     
  
  for(i in 1:nboot) {
    if(!is.null(list.tindx)) {
      tindx <- list.tindx[[i]]
      if(tindx > N) warning("number of selected observations is larger than the sample size")
    } else {
      tindx <- sample(nindx, N, replace = TRUE)
    }
    #tindx <- ifelse(!is.null(list.tindx), list.tindx[[i]], sample(nindx, N, replace = TRUE))
    #tindx <- sample(nindx, N, replace = TRUE)
    mymodel <- model(formula, data=data[tindx,], ...)
    outbootdata <- subset(data, !(nindx %in% tindx))
    # check if mymodel is a function which should be used instead of
    # predict
    if (is.function(mymodel)) {
      if(!is.null(predict) & i == 1) 
        warning("model returns a function and predict is specified, using models output")
      predict <- mymodel
      USEPM <- TRUE
    }

    if (USEPM)
      pred <- predict(newdata=outbootdata)
    else
      pred <- predict(mymodel, newdata = outbootdata)

    if (is.list(pred)) {
      if (!inherits(pred[[1]], "survfit") && !inherits(pred, "survfit"))
        stop("predict does not return a list of survfit objects")
    } else {
      stop("predict does not return a list of survfit objects")
    }
    if(predictions) BOOTINDX[-tindx, i] <- sbrier(y[-tindx], pred) ###???
    bootindx[-tindx, i] <- sbrier(y[-tindx], pred)
  }

  fun <- function(x)
        ifelse(all(is.na(x)), NA, mean(x, na.rm = TRUE))

  err <- mean(apply(bootindx, 1, fun), na.rm = TRUE)
  if(predictions) 
    RET <- list(error = err, nboot = nboot, predictions = BOOTINDX)
  else
    RET <- list(error = err, nboot=nboot)
  class(RET) <- "bootestsurv"
  RET
}

