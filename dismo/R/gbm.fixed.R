#
# copyright  j leathwick, j elith - 6th May 2007
#
# version 2.9 - developed in R 2.0
#
# calculates a gradient boosting (gbm)object with a fixed number of trees
# with the number of trees identified using gbm.step or some other procedure
#
# mostly used as a utility function, e.g., when being called by gbm.simplify
#
# takes as input a dataset and args selecting x and y variables, learning rate and tree complexity
#
# updated 13/6/05 to accommodate weighting of sites when calculating total and residual deviance
#
# updated 10/8/05 to correct how site.weights are returned
#
# requires gbm
#
#

gbm.fixed <- 
function (data,                        # the input dataframe
  gbm.x,                               # indices of the predictors in the input dataframe
  gbm.y,                               # index of the response in the input dataframe
  tree.complexity = 1,                 # the tree depth - sometimes referred to as interaction depth
  site.weights = rep(1, nrow(data)),   # by default set equal
  verbose = TRUE,                      # to control reporting
  learning.rate = 0.001,               # controls speed of the gradient descent
  n.trees = 2000,                      # default number of trees
  bag.fraction = 0.5,                  # varies random sample size for each new tree
  family = "bernoulli",                # can be any of bernoulli, poisson, gaussian, laplace - note quotes
  keep.data = FALSE,                   # keep original data
  var.monotone = rep(0, length(gbm.x)) # constrain to positive (1) or negative monontone (-1)
  ) 
{
    train.fraction = 1

    if (! requireNamespace('gbm') ) { stop ('you need to install the gbm package to run this function') }
	
	
# setup input data and assign to position one

#  dataframe.name <- deparse(substitute(data))   # get the dataframe name

  x.data <- data[, gbm.x, drop=FALSE]                 #form the temporary datasets
  y.data <- data[, gbm.y]
  sp.name <- names(data)[gbm.y]
    

  # assign("x.data", x.data, pos = 1)             #and assign them for later use
  # assign("y.data", y.data, pos = 1)

# fit the gbm model 

  z1 <- unclass(Sys.time())

  gbm.call <- paste("gbm::gbm(y.data ~ .,n.trees = n.trees, data=x.data, verbose = F, interaction.depth = tree.complexity, 
    weights = site.weights, shrinkage = learning.rate, distribution = as.character(family),
    var.monotone = var.monotone, bag.fraction = bag.fraction, keep.data = keep.data)", sep="")

  if (verbose) {
    print(paste("fitting gbm model with a fixed number of ", n.trees, " trees for ", sp.name,sep=""),quote=FALSE) }

  gbm.object <- eval(parse(text = gbm.call))

  best.trees <- n.trees

#extract fitted values and summary table
  
  fitted.values <- gbm::predict.gbm(gbm.object,x.data,n.trees = n.trees,type="response")
  gbm.summary <- summary(gbm.object,n.trees = n.trees, plotit = FALSE)

  y_i <- y.data
  u_i <- fitted.values

  if (family == "poisson") {
    deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
    resid.deviance <- 2 * sum(deviance.contribs * site.weights)
    residuals <- sqrt(abs(deviance.contribs * 2))
    residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)

    u_i <- sum(y.data * site.weights) / sum(site.weights)
    deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
    total.deviance <- 2 * sum(deviance.contribs * site.weights)
  }

  if (family == "bernoulli") {
    deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
    resid.deviance <- -2 * sum(deviance.contribs * site.weights)
    residuals <- sqrt(abs(deviance.contribs * 2))
    residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)

    u_i <- sum(y.data * site.weights) / sum(site.weights)
    deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
    total.deviance <- -2 * sum(deviance.contribs * site.weights)
  }

  if (family == "laplace") {
    resid.deviance <- sum(abs(y_i - u_i))
    residuals <- y_i - u_i
    u_i <- mean(y.data)
    total.deviance <- sum(abs(y_i - u_i))
  }

  if (family == "gaussian") {
    resid.deviance <- sum((y_i - u_i) * (y_i - u_i))
    residuals <- y_i - u_i
    u_i <- mean(y.data)
    total.deviance <- sum((y_i - u_i) * (y_i - u_i))
  }

  if (verbose) {
    print(paste("total deviance = ",round(total.deviance,2),sep=""),quote=F)
    print(paste("residual deviance = ",round(resid.deviance,2),sep=""),quote=F)}

# now assemble data to be returned

  z2 <- unclass(Sys.time())
  elapsed.time.minutes <- round((z2 - z1)/ 60,2)  #calculate the total elapsed time

  gbm.detail <- list(data = data, gbm.x = gbm.x, predictor.names = names(x.data), 
    gbm.y = gbm.y, reponse.name = names(y.data), family = family, tree.complexity = tree.complexity, 
    learning.rate = learning.rate, bag.fraction = bag.fraction, cv.folds = 0, n.trees = n.trees, best.trees = best.trees, 
    train.fraction = train.fraction, var.monotone = var.monotone, date = date(), elapsed.time.minutes = elapsed.time.minutes)

  gbm.object$gbm.call <- gbm.detail
  gbm.object$fitted <- fitted.values
  gbm.object$residuals <- residuals
  gbm.object$contributions <- gbm.summary
  gbm.object$self.statistics <- list(null.deviance = total.deviance, resid.deviance = resid.deviance)
  gbm.object$weights <- site.weights

#  rm(x.data,y.data, pos=1)           #finally, clean up the temporary dataframes

  return(gbm.object)
}


