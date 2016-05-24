cvsmac <- function(x,y, kfold = 5, lambda = NULL, nlambda = 100, lambda.min = NULL, seed = 0, weight = NULL,...){
  if(length(y) != nrow(x)){
    stop("The dimension of covariates must be equal.")
  }
  if(is.na(sum(x))|is.nan(sum(x))){
    stop('The matrix/data frame should not contain NA/NaN.')
  }
  
  set.seed(seed)

  np <- ncol(x)
  nobs <- nrow(x)
  if (length(y)!=nobs) {stop("The dimension of covariates should match the length of the label.")}
  if (is.na(nlambda) | is.nan(nlambda)) stop("nlambda should not be NA/NaN.")
  if (!is.numeric(nlambda)) {stop("nlambda should be numeric.")}
  if (nlambda!=round(nlambda)) {stop("nlambda should be an integer.")}
  if (nlambda<0) {stop("nlambda should be at least 1.")}
  if (length(as.vector(nlambda))>1) {stop("nlambda should be a scalar.")}
  
  if (kfold <= 0){
    stop('The value for kfold should be positive.')
  }
  if ((kfold%%1) != 0){
    stop('The value for kfold should be an integer.')
  }
  
  for(i in levels(as.factor(y))){
    if(length(y[which(y == i)]) == 1){
      stop(paste0('Class ', i, ' has only one observation.'))
    }
    if(length(y[which(y == i)]) < kfold){
	warns=paste("Class ",i," has less observation than kfold.",sep="")
      stop(warns)
      }
    if(length(y[which(y == i)])%/% kfold < 5){
      warns=paste("Class ",i," has too few observations.",sep="")
      warning(warns)
    }
  }
  
  # Finding lambda's in case lambda is not specified
  if(is.null(lambda.min)){lambda.min = ifelse(nobs < np, 0.05, 1e-03)}
  
  # Normalising Weight
  if(is.null(weight)){
	weight <- rep(x = 1, times = length(y))
	} else {
  weight <- (weight/sqrt(mean(weight^2)))
	warning("The weight vector is re-scaled.")
	}

  # Getting cvlambda
  if(is.null(lambda)){
    a <- smac(y=y,x=x, nlambda = 100, weight = weight, ...)
    cvlambda <- a$lambda
  } else {
    cvlambda <- lambda
    a <- smac(y=y,x=x, lambda=cvlambda, weight = weight, ...)
  }
  
  # Partitioning the data set for cross validation
  x_nlst <- paste('x',1:kfold,sep = '')
  y_nlst <- paste('y',1:kfold,sep = '')
  w_nlst <- paste('w',1:kfold,sep = '')
  # Creating list containing NULL
  for(i in 1:kfold){
    x_lst <- sapply(x_nlst, function(x){NULL})
    y_lst <- sapply(y_nlst, function(x){NULL})
    w_lst <- sapply(w_nlst, function(x){NULL})
  }
  
  # Checking if y is numeric. If y is not numeric, then
  # convert y to numeric values
  if(!is.numeric(y)){
    # Make a copy of the data
    y.factor <- as.factor(y)
    y.temp <- NULL
    k = 1
    for(i in levels(y.factor)){
      y.temp[which(y == i)] <- k
      k <- k+1
    }
    
    x.copy <- x; y.copy <- as.factor(y.temp)
  } else if (is.numeric(y)){
    # Make a copy of the data
    x.copy <- x; y.copy <- as.factor(y)
  }
  
  # Name Holder
  names <- paste('v',1:length(levels(y.copy)),sep = '')
  vec <- NULL
  for(nj in names){
    j <- levels(y.copy)[as.numeric(substr(nj,start = 2, stop = 3))]
    vecj <- assign(x = nj, value = which(y.copy == j))
    vec <- c(vec,sample(x = vecj, size = length(vecj)))
  }
  

  # for loop to partition data
  for(i in 1:length(vec)){
    j <- i%%kfold + 1
      y_lst[[j]] <- c(y_lst[[j]], y[vec[i]])
      x_lst[[j]] <- rbind(x_lst[[j]],x[vec[i],])
      w_lst[[j]] <- c(w_lst[[j]],weight[vec[i]])
      }
  
  # Calculating the errors
  err <- rep(x = 0, times = length(cvlambda)) # vector containing each error for each lambda
  
  for(i in 1:kfold){
    x.train <- NULL; y.train <- NULL; w.train <- NULL
    hold <- which(1:kfold != i)
    for(j in hold){
      x.train <- rbind(x.train,x_lst[[j]])
      y.train <- c(y.train,y_lst[[j]])
      w.train <- c(w.train,w_lst[[j]])
    }
    x.test <- x_lst[[i]]
    y.test <- y_lst[[i]]
    w.test <- w_lst[[i]]

    model <- smac(x = x.train, y = y.train, lambda = cvlambda, weight = w.train,...)
    pred <- predict(model, new.x = x.test, lambda = cvlambda)
    
    err <- err + sapply(pred$pred.y,function(x){
      (sum(w.test*as.numeric(y.test!=x)))
    })
  }
  
  # Finding minimum lambda and best models
  min_err <- min(err) 
  best_lam <- cvlambda[which(err == min_err)]
  best <- list(lambda = cvlambda, beta0 = a$beta0, beta = a$beta, error = err, best.lambda = best_lam, best.beta0 = a$beta0[which(cvlambda %in% best_lam)], best.beta = a$beta[which(cvlambda %in% best_lam)], model = a, min.error = min_err)
  
  this.call = match.call()  
  best$call <- this.call  
  class(best) <- 'cvsmac'

  return(best)
}
