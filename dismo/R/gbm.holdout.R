#
# j leathwick, j elith - October 2006
#
# version 2.9 - developed in R 2.3.1
#
# calculates a gradient boosting (gbm)object in which model complexity is 
# determined using a training set with predictions made to a withheld set
# an initial set of trees is fitted, and then trees are progressively added
# testing performance # along the way, using gbm::gbm.perf until the optimal
# number of trees is identified
#
# as any structured ordering of the data should be avoided, a copy of the data set 
# BY DEFAULT is randomly reordered each time the function is run
#
# takes as input a dataset and args selecting x and y variables, and degree of interaction depth
#
# requires gbm
#

gbm.holdout <-
function (data,                        # the input data frame
   gbm.x,                              # indices of predictor variables
   gbm.y,                              # index of response variable
   learning.rate = 0.001,              # typically varied between 0.1 and 0.001
   tree.complexity = 1,                # sometimes called interaction depth
   family = "bernoulli",               # "bernoulli","poisson", etc. as for gbm
   n.trees = 200,                      # initial number of trees
   add.trees = n.trees,                # number of trees to add at each increment
   max.trees = 20000,                  # maximum number of trees to fit
   verbose = TRUE,                     # controls degree of screen reporting
   train.fraction = 0.8,               # proportion of data to use for training
   permute = TRUE,                     # reorder data to start with
   prev.stratify = TRUE,               # stratify selection for p/a data
   var.monotone = rep(0, length(gbm.x)),# allows constraining of response to monotone 
   site.weights = rep(1, nrow(data)),  # set equal to 1 by default
   refit = TRUE,                       # refit the model with the full data but id'd no of trees
   keep.data = TRUE)                   # keep copy of the data
{
#
    if (! requireNamespace('gbm') ) { stop ('you need to install the gbm package to run this function') }
	requireNamespace('splines')
# setup input data and assign to position one

  cv.folds <- 0 

  if (permute) { 
    print("",quote=FALSE)
    print("WARNING - data is being randomly reordered to avoid confounding effects",quote=FALSE)
    print("of inherent structure as submitted - use permute = FALSE to turn off this option",quote=FALSE)
    n.rows <- nrow(data)

    if (prev.stratify == TRUE & family == "bernoulli") {

      presence.mask <- data[,gbm.y] == 1
      absence.mask <- data[,gbm.y] == 0
      n.pres <- sum(presence.mask)
      n.abs <- sum(absence.mask)

      selector <- seq(1,n.rows)

      temp <- sample(selector[presence.mask],size = n.pres * train.fraction)
      selector[temp] <- 0

      temp <- sample(selector[absence.mask],size = n.abs * train.fraction)
      selector[temp] <- 0 

      sort.vector <- sort(selector,index.return = TRUE)[[2]]
    }

    else {
      sort.vector <- sample(seq(1,n.rows),n.rows,replace=FALSE)
    }

    sort.data <- data[sort.vector,]
    x.data <- sort.data[, gbm.x, drop=FALSE]                 #form the temporary datasets
    y.data <- sort.data[, gbm.y]
  } else {
    x.data <- data[, gbm.x, drop=FALSE]                 #form the temporary datasets
    y.data <- data[, gbm.y]
  }

  names(x.data) <- names(data)[gbm.x]
  sp.name <- names(data)[gbm.y]
    
#  assign("x.data", x.data, pos = 1)               #and assign them for later use
#  assign("y.data", y.data, pos = 1)

# fit the gbm model 

  print(paste("fitting initial gbm model of ",n.trees," trees for ",sp.name,sep=""),quote=FALSE)
  print(" and expanding using withheld data for evaluation",quote=FALSE) 

  gbm.call <- paste("gbm::gbm(y.data ~ .,n.trees = n.trees, data=x.data, verbose = F, interaction.depth = tree.complexity, 
    weights = site.weights, shrinkage = learning.rate, cv.folds = 0, distribution = as.character(family),
    train.fraction = train.fraction, var.monotone = var.monotone, keep.data = keep.data)", sep="")

  gbm.object <- eval(parse(text = gbm.call))

# identify the best number of trees using method appropriate to model

  best.trees <- gbm::gbm.perf(gbm.object, method = 'test', plot.it = FALSE) 

  n.fitted <- n.trees

  if (verbose) print("expanding model to find optimal no of trees...",quote=FALSE)

  while(gbm.object$n.trees - best.trees < n.trees & n.fitted < max.trees){

    gbm.object <- gbm::gbm.more(gbm.object, add.trees)
    best.trees <- gbm::gbm.perf(gbm.object, method = 'test', plot.it = FALSE)
    n.fitted <- n.fitted + add.trees

    if (n.fitted %% 100 == 0){ #report times along the way
      if (verbose) print(paste("fitted trees = ", n.fitted, sep = ""), quote = FALSE)
    }
  }

  if (verbose) print(paste("fitting stopped at ",best.trees," trees",sep=""),quote=FALSE)

  if (refit) {  # we are refitting the model with fixed tree size
    print(paste("refitting the model to the full dataset using ",best.trees," trees",sep=""),quote=FALSE)

    x.data <- data[, gbm.x, drop=FALSE]                 #form the temporary datasets
    y.data <- data[, gbm.y]
     
    gbm.call <- eval(paste("gbm::gbm(y.data ~ .,n.trees = best.trees, data=x.data, verbose = F, interaction.depth = tree.complexity, 
      weights = site.weights, shrinkage = learning.rate, cv.folds = 0, distribution = as.character(family),
      var.monotone = var.monotone, keep.data = keep.data)", sep=""))

    gbm.object <- eval(parse(text = gbm.call))

  }

#extract fitted values and summary table
  
  fitted.values <- gbm::predict.gbm(gbm.object,x.data,n.trees = best.trees,type="response")
  gbm.summary <- summary(gbm.object,n.trees = best.trees, plotit = FALSE)

  y_i <- y.data
  u_i <- fitted.values

  if (family == "poisson") {
    deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
    resid.deviance <- 2 * sum(deviance.contribs)
    residuals <- sqrt(abs(deviance.contribs * 2))
    residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)
    u_i <- mean(y.data)
    total.deviance <- 2 * sum(ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i))
  }

  if (family == "bernoulli") {
    deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
    resid.deviance <- -2 * sum(deviance.contribs)
    residuals <- sqrt(abs(deviance.contribs * 2))
    residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)
    u_i <- mean(y.data)
    total.deviance <- -2 * sum((y_i * log(u_i)) + ((1-y_i) * log(1 - u_i)))
  }

  if (verbose) {
    print(paste("total deviance = ",round(total.deviance,2),sep=""),quote=F)
    print(paste("residual deviance = ",round(resid.deviance,2),sep=""),quote=F)
  }

# now assemble data to be returned

  gbm.detail <- list(dataframe = data, gbm.x = gbm.x, predictor.names = names(x.data), 
    gbm.y = gbm.y, response.name = sp.name, tree.complexity = tree.complexity, n.trees = best.trees, 
    learning.rate = learning.rate, best.trees = best.trees, cv.folds = cv.folds, 
    family = family, train.fraction = train.fraction, var.monotone = var.monotone )

  gbm.object$fitted <- fitted.values
  gbm.object$residuals <- residuals
  gbm.object$contributions <- gbm.summary
  gbm.object$deviances <- list(null.deviance = total.deviance, resid.deviance = resid.deviance) 
  gbm.object$weights <- weights
  gbm.object$gbm.call <- gbm.detail

#  rm(x.data,y.data, pos=1)           #finally, clean up the temporary dataframes

  return(gbm.object)
}


