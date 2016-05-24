#
# j. leathwick/j. elith - 19th September 2005
#
# version 2.9
#
# function to assess optimal no of boosting trees using k-fold cross validation 
#
# implements the cross-validation procedure described on page 215 of 
# Hastie T, Tibshirani R, Friedman JH (2001) The Elements of Statistical Learning: 
# Data Mining, Inference, and Prediction Springer-Verlag, New York.
#
# divides the data into 10 subsets, with stratification by prevalence if required for pa data
# then fits a gbm model of increasing complexity along the sequence from n.trees to n.trees + (n.steps * step.size)
# calculating the residual deviance at each step along the way
# after each fold processed, calculates the average holdout residual deviance and its standard error
# then identifies the optimal number of trees as that at which the holdout deviance is minimised 
# and fits a model with this number of trees, returning it as a gbm model along with additional information 
# from the cv selection process
#
# updated 13/6/05 to accommodate weighting of sites
#
# updated 19/8/05 to increment all folds simultaneously, allowing the stopping rule
# for the maxinum number of trees to be fitted to be imposed by the data, 
# rather than being fixed in advance
#
# updated 29/8/05 to return cv test statistics, and deviance as mean
# time for analysis also returned via unclass(Sys.time())
#
# updated 5/9/05 to use external function calc.deviance 
# and to return cv test stats via predictions formed from fold models 
# with n.trees = target.trees
#
# updated 15/5/06 to calculate variance of fitted and predicted values across folds
# these can be expected to approximate the variance of fitted values
# as would be estimated for example by bootstrapping
# as these will underestimate the true variance 
# they are corrected by multiplying by (n-1)2/n 
# where n is the number of folds
#
# updated 25/3/07 tp allow varying of bag fraction
#
# requires gbm library from Cran
# requires roc and calibration scripts of J Elith
# requires calc.deviance script of J Elith/J Leathwick
#
#


gbm.step <- function (
  data,                                     # the input dataframe
  gbm.x,                                    # the predictors
  gbm.y,                                    # and response
  offset = NULL,                            # allows an offset to be specified
  fold.vector = NULL,                       # allows a fold vector to be read in for CV with offsets,
  tree.complexity = 1,                      # sets the complexity of individual trees
  learning.rate = 0.01,                     # sets the weight applied to inidivudal trees
  bag.fraction = 0.75,                      # sets the proportion of observations used in selecting variables
  site.weights = rep(1, nrow(data)),        # allows varying weighting for sites
  var.monotone = rep(0, length(gbm.x)),     # restricts responses to individual predictors to monotone 
  n.folds = 10,                             # number of folds
  prev.stratify = TRUE,                     # prevalence stratify the folds - only for p/a data
  family = "bernoulli",                     # family - bernoulli (=binomial), poisson, laplace or gaussian
  n.trees = 50,                             # number of initial trees to fit
  step.size = n.trees,                      # numbers of trees to add at each cycle
  max.trees = 10000,                        # max number of trees to fit before stopping
  tolerance.method = "auto",                # method to use in deciding to stop - "fixed" or "auto"
  tolerance = 0.001,                        # tolerance value to use - if method == fixed is absolute, 
                                            # if auto is multiplier * total mean deviance
  plot.main = TRUE,                         # plot hold-out deviance curve
  plot.folds = FALSE,                       # plot the individual folds as well
  verbose = TRUE,                           # control amount of screen reporting
  silent = FALSE,                           # to allow running with no output for simplifying model)
  keep.fold.models = FALSE,                 # keep the fold models from cross valiation
  keep.fold.vector = FALSE,                 # allows the vector defining fold membership to be kept
  keep.fold.fit = FALSE,                    # allows the predicted values for observations from CV to be kept
  ...)                                      # allows for any additional plotting parameters
{

	if (! requireNamespace('gbm') ) { stop ('you need to install the gbm package to run this function') }
	requireNamespace('splines')

	if (silent) verbose <- FALSE

# initiate timing call

	z1 <- Sys.time()

# setup input data and assign to position one

#  dataframe.name <- deparse(substitute(data))   # get the dataframe name
#  data <- eval(data)
	x.data <- data[, gbm.x, drop=FALSE]                 #form the temporary datasets
	#names(x.data) <- names(data)[gbm.x]
	y.data <- data[, gbm.y]
	sp.name <- names(data)[gbm.y]
	if (family == "bernoulli") {
		prevalence <- mean(y.data)
	}

#  assign("x.data", x.data, env = globalenv())               #and assign them for later use
#  assign("y.data", y.data, env = globalenv())

	#offset.name <- deparse(substitute(offset))   # get the dataframe name
	#offset <- eval(offset)

	n.cases <- nrow(data)
	n.preds <- length(gbm.x)

	if (!silent) {
		cat("\n","\n","GBM STEP - version 2.9","\n","\n")
		cat("Performing cross-validation optimisation of a boosted regression tree model \n")
		cat("for", sp.name, "and using a family of",family,"\n")
		cat("Using",n.cases,"observations and",n.preds,"predictors \n")
	}

# set up the selector variable either with or without prevalence stratification

	if (is.null(fold.vector)) {

		if (prev.stratify & family == "bernoulli") {
			presence.mask <- data[,gbm.y] == 1
			absence.mask <- data[,gbm.y] == 0
			n.pres <- sum(presence.mask)
			n.abs <- sum(absence.mask)

# create a vector of randomised numbers and feed into presences
			selector <- rep(0, n.cases)
			temp <- rep(seq(1, n.folds, by = 1), length = n.pres)
			temp <- temp[order(runif(n.pres, 1, 100))]
			selector[presence.mask] <- temp

# and then do the same for absences
			temp <- rep(seq(1, n.folds, by = 1), length = n.abs)
			temp <- temp[order(runif(n.abs, 1, 100))]
			selector[absence.mask] <- temp
		} else {  #otherwise make them random with respect to presence/absence
			selector <- rep(seq(1, n.folds, by = 1), length = n.cases)
			selector <- selector[order(runif(n.cases, 1, 100))]
		}
    } else {
		if (length(fold.vector) != n.cases) {
			stop("supplied fold vector is of wrong length")
		}
		cat("loading user-supplied fold vector \n")
		selector <- fold.vector
    }

# set up the storage space for results

	pred.values <- rep(0, n.cases)

	cv.loss.matrix <- matrix(0, nrow = n.folds, ncol = 1)
	training.loss.matrix <- matrix(0, nrow = n.folds, ncol = 1)
	trees.fitted <- n.trees

	model.list <- list(paste("model", 1:n.folds, sep=""))     # dummy list for the tree models

# set up the initial call to gbm

	if (is.null(offset)) {
		gbm.call <- paste("gbm::gbm(y.subset ~ .,data=x.subset, n.trees = n.trees, interaction.depth = tree.complexity, shrinkage = learning.rate, bag.fraction = bag.fraction, weights = weight.subset, distribution = as.character(family), var.monotone = var.monotone, verbose = FALSE)", sep="")
    } else {
		gbm.call <- paste("gbm::gbm(y.subset ~ . + offset(offset.subset), data=x.subset, n.trees = n.trees, interaction.depth = tree.complexity, shrinkage = learning.rate, bag.fraction = bag.fraction, weights = weight.subset, distribution = as.character(family), var.monotone = var.monotone, verbose = FALSE)", sep="")
    } 

	n.fitted <- n.trees

# calculate the total deviance
               
	y_i <- y.data

	u_i <- sum(y.data * site.weights) / sum(site.weights)
	u_i <- rep(u_i,length(y_i))

	total.deviance <- calc.deviance(y_i, u_i, weights = site.weights, family = family, calc.mean = FALSE)

	mean.total.deviance <- total.deviance/n.cases

	tolerance.test <- tolerance

	if (tolerance.method == "auto") {
		tolerance.test <- mean.total.deviance * tolerance
	}

# now step through the folds setting up the initial call

	if (!silent){ 
		cat("creating",n.folds,"initial models of",n.trees,"trees","\n")
		if (prev.stratify & family == "bernoulli") {
			cat("\n","folds are stratified by prevalence","\n")
		} else {
			cat("\n","folds are unstratified","\n")
		}
		cat ("total mean deviance = ",round(mean.total.deviance,4),"\n")
		cat("tolerance is fixed at ",round(tolerance.test,4),"\n")
		if (tolerance.method != "fixed" & tolerance.method != "auto") {
			stop("invalid argument for tolerance method - should be auto or fixed")
		}
	}

	if (verbose) {
		cat("ntrees resid. dev.","\n")
	}
	
	for (i in 1:n.folds) {

		model.mask <- selector != i  #used to fit model on majority of data
		pred.mask <- selector == i   #used to identify the with-held subset

		y.subset <- y.data[model.mask]
		x.subset <- x.data[model.mask, ,drop=FALSE]
		weight.subset <- site.weights[model.mask]

		if (!is.null(offset)) {
			offset.subset <- offset[model.mask] 
		} else {
			offset.subset <- NULL
		}

		model.list[[i]] <- eval(parse(text = gbm.call))

		fitted.values <- model.list[[i]]$fit  #predict.gbm(model.list[[i]], x.subset, type = "response", n.trees = n.trees)
		if (!is.null(offset)) {
			fitted.values <- fitted.values + offset[model.mask]
		}
		if (family == "bernoulli") {
			fitted.values <- exp(fitted.values)/(1 + exp(fitted.values))
		} else if (family == "poisson") {
			fitted.values <- exp(fitted.values)
		}
		
		pred.values[pred.mask] <- gbm::predict.gbm(model.list[[i]], x.data[pred.mask, ,drop=FALSE], n.trees = n.trees)
		
		if (!is.null(offset)) {
			pred.values[pred.mask] <- pred.values[pred.mask] + offset[pred.mask]
		}
		if (family == "bernoulli") {
			pred.values[pred.mask] <- exp(pred.values[pred.mask])/(1 + exp(pred.values[pred.mask]))
		} else if (family == "poisson") {
			pred.values[pred.mask] <- exp(pred.values[pred.mask])
		}

	# calc training deviance

		y_i <- y.subset
		u_i <- fitted.values
		weight.fitted <- site.weights[model.mask]
		training.loss.matrix[i,1] <- calc.deviance(y_i, u_i, weight.fitted, family = family)

# calc holdout deviance

		y_i <- y.data[pred.mask]
		u_i <- pred.values[pred.mask]
		weight.preds <- site.weights[pred.mask]
		cv.loss.matrix[i,1] <- calc.deviance(y_i, u_i, weight.preds, family = family)

	} # end of first loop

# now process until the change in mean deviance is =< tolerance or max.trees is exceeded

	delta.deviance <- 1

	cv.loss.values <- apply(cv.loss.matrix,2,mean)
	if (verbose) {
		cat(n.fitted,"  ",round(cv.loss.values,4),"\n")
	}
	if (!silent) {
		cat("now adding trees...","\n")
	}
	
	j <- 1

	while (delta.deviance > tolerance.test & n.fitted < max.trees) { 
		# beginning of inner loop

		# add a new column to the results matrice..

		training.loss.matrix <- cbind(training.loss.matrix,rep(0,n.folds))
		cv.loss.matrix <- cbind(cv.loss.matrix,rep(0,n.folds))

		n.fitted <- n.fitted + step.size
		trees.fitted <- c(trees.fitted,n.fitted)
  
		j <- j + 1
 
		for (i in 1:n.folds) {

			model.mask <- selector != i  #used to fit model on majority of data
			pred.mask <- selector == i   #used to identify the with-held subset

			y.subset <- y.data[model.mask]
			x.subset <- x.data[model.mask, ,drop=FALSE]
			weight.subset <- site.weights[model.mask]
			if (!is.null(offset)) {
				offset.subset <- offset[model.mask] 
			}
			model.list[[i]] <- gbm::gbm.more(model.list[[i]], weights = weight.subset, step.size)

			fitted.values <- model.list[[i]]$fit # predict.gbm(model.list[[i]],x.subset, type = "response", n.trees = n.fitted) 
			if (!is.null(offset)) {
				fitted.values <- fitted.values + offset[model.mask]
			}
			if (family == "bernoulli") {
				fitted.values <- exp(fitted.values)/(1 + exp(fitted.values))
			} else if (family == "poisson") {
				fitted.values <- exp(fitted.values)
			}
			pred.values[pred.mask] <- gbm::predict.gbm(model.list[[i]], x.data[pred.mask, ,drop=FALSE], n.trees = n.fitted)
			
			if (!is.null(offset)) {
				pred.values[pred.mask] <- pred.values[pred.mask] + offset[pred.mask]
			}
			
			if (family == "bernoulli") {
				pred.values[pred.mask] <- exp(pred.values[pred.mask])/(1 + exp(pred.values[pred.mask]))
			} else if (family == "poisson") {
				pred.values[pred.mask] <- exp(pred.values[pred.mask])
			}
# calculate training deviance

			y_i <- y.subset
			u_i <- fitted.values
			weight.fitted <- site.weights[model.mask]
			training.loss.matrix[i,j] <- calc.deviance(y_i, u_i, weight.fitted, family = family)

# calc holdout deviance

			u_i <- pred.values[pred.mask]
			y_i <- y.data[pred.mask]
			weight.preds <- site.weights[pred.mask]
			cv.loss.matrix[i,j] <- calc.deviance(y_i, u_i, weight.preds, family = family)

		}  # end of inner loop

		cv.loss.values <- apply(cv.loss.matrix,2,mean)

		if (j < 5) {
			if (cv.loss.values[j] > cv.loss.values[j-1]) {
				if (!silent) {
					cat("restart model with a smaller learning rate or smaller step size...")
				}
				return()
			}
		}

		if (j >= 20) {   #calculate stopping rule value
			test1 <- mean(cv.loss.values[(j-9):j])
			test2 <- mean(cv.loss.values[(j-19):(j-9)])
			delta.deviance <- test2 - test1
		}

		if (verbose) {
			cat(n.fitted," ",round(cv.loss.values[j],4),"\n") 
			flush.console()
		}
	} # end of while loop

# now begin process of calculating optimal number of trees

	training.loss.values <- apply(training.loss.matrix,2,mean)

	cv.loss.ses <- rep(0,length(cv.loss.values))
	cv.loss.ses <- sqrt(apply(cv.loss.matrix,2,var)) / sqrt(n.folds)

# find the target holdout deviance

	y.bar <- min(cv.loss.values) 


# identify the optimal number of trees 

	target.trees <- trees.fitted[match(TRUE, cv.loss.values == y.bar)]

# plot out the resulting curve of holdout deviance 
	if (plot.main) {

		y.min <- min(cv.loss.values - cv.loss.ses)  #je added multiplier 10/8/05
		y.max <- max(cv.loss.values + cv.loss.ses)  #je added multiplier 10/8/05 }

		if (plot.folds) {
			y.min <- min(cv.loss.matrix)
			y.max <- max(cv.loss.matrix) 
		}

		plot(trees.fitted, cv.loss.values, type = 'l', ylab = "holdout deviance", xlab = "no. of trees", ylim = c(y.min,y.max), ...)
		abline(h = y.bar, col = 2)

		lines(trees.fitted, cv.loss.values + cv.loss.ses, lty=2)  
		lines(trees.fitted, cv.loss.values - cv.loss.ses, lty=2)  

		if (plot.folds) {
			for (i in 1:n.folds) {
				lines(trees.fitted, cv.loss.matrix[i,],lty = 3)
			}
		}
		abline(v = target.trees, col=3)
		title(paste(sp.name,", d - ",tree.complexity,", lr - ",learning.rate, sep=""))
	}

# estimate the cv deviance and test statistics
# includes estimates of the standard error of the fitted values added 2nd may 2005

	cv.deviance.stats <- rep(0, n.folds)
	cv.roc.stats <- rep(0, n.folds)
	cv.cor.stats <- rep(0, n.folds)
	cv.calibration.stats <- matrix(0, ncol=5, nrow = n.folds)
	if (family == "bernoulli") {
		threshold.stats <- rep(0, n.folds)
	}
	fitted.matrix <- matrix(NA, nrow = n.cases, ncol = n.folds)  # used to calculate se's
	fold.fit <- rep(0, n.cases)

	for (i in 1:n.folds) {

		pred.mask <- selector == i   #used to identify the with-held subset
		model.mask <- selector != i  #used to fit model on majority of data

		fits <- gbm::predict.gbm(model.list[[i]], x.data[model.mask, ,drop=FALSE], n.trees = target.trees)
		if (!is.null(offset)) {
			fits <- fits + offset[model.mask]
		}
		if (family == "bernoulli") {
			fits <- exp(fits)/(1 + exp(fits))
		} else if (family == "poisson") {
			fits <- exp(fits)
		}
		fitted.matrix[model.mask,i] <- fits

		fits <- gbm::predict.gbm(model.list[[i]], x.data[pred.mask, ,drop=FALSE], n.trees = target.trees)
		if (!is.null(offset)) fits <- fits + offset[pred.mask]
		fold.fit[pred.mask] <- fits  # store the linear predictor values
		if (family == "bernoulli") {
			fits <- exp(fits)/(1 + exp(fits))
		} else if (family == "poisson") {
			fits <- exp(fits)
		}
		fitted.matrix[pred.mask,i] <- fits

		y_i <- y.data[pred.mask] 
		u_i <- fitted.matrix[pred.mask,i]  #pred.values[pred.mask]
		weight.preds <- site.weights[pred.mask]

		cv.deviance.stats[i] <- calc.deviance(y_i, u_i, weight.preds, family = family)

		cv.cor.stats[i] <- cor(y_i,u_i)

		if (family == "bernoulli") {
			cv.roc.stats[i] <- .roc(y_i,u_i)
			cv.calibration.stats[i,] <- .calibration(y_i,u_i,"binomial")
			threshold.stats[i] <- approx(ppoints(u_i), sort(u_i,decreasing = T), prevalence)$y
		}

		if (family == "poisson") {
			cv.calibration.stats[i,] <- .calibration(y_i,u_i,"poisson")
		}
	}

	fitted.vars <- apply(fitted.matrix,1, var, na.rm = TRUE)

# now calculate the mean and se's for the folds

	cv.dev <- mean(cv.deviance.stats, na.rm = TRUE)
	cv.dev.se <- sqrt(var(cv.deviance.stats)) / sqrt(n.folds)

	cv.cor <- mean(cv.cor.stats, na.rm = TRUE)
	cv.cor.se <- sqrt(var(cv.cor.stats, use = "complete.obs")) / sqrt(n.folds)

	cv.roc <- 0.0
	cv.roc.se <- 0.0 

	if (family == "bernoulli") {
		cv.roc <- mean(cv.roc.stats,na.rm=TRUE)
		cv.roc.se <- sqrt(var(cv.roc.stats, use = "complete.obs")) / sqrt(n.folds)  
		cv.threshold <- mean(threshold.stats, na.rm = T)
		cv.threshold.se <- sqrt(var(threshold.stats, use = "complete.obs")) / sqrt(n.folds)  
	}
 
	cv.calibration <- 0.0
	cv.calibration.se <- 0.0

	if (family == "poisson" | family == "bernoulli") {
		cv.calibration <- apply(cv.calibration.stats,2,mean)
		cv.calibration.se <- apply(cv.calibration.stats,2,var)
		cv.calibration.se <- sqrt(cv.calibration.se) / sqrt(n.folds) 
	}

# fit the final model

	if (is.null(offset)) {
		gbm.call <- paste("gbm::gbm(y.data ~ .,data=x.data, n.trees = target.trees, interaction.depth = tree.complexity, shrinkage = learning.rate, bag.fraction = bag.fraction, weights = site.weights, distribution = as.character(family), var.monotone = var.monotone, verbose = FALSE)", sep="")
	} else {
		gbm.call <- paste("gbm::gbm(y.data ~ . + offset(offset),data=x.data, n.trees = target.trees, interaction.depth = tree.complexity, shrinkage = learning.rate, bag.fraction = bag.fraction, weights = site.weights, distribution = as.character(family), var.monotone = var.monotone,  verbose = FALSE)", sep="")
    } 

	if (!silent) {
		message("fitting final gbm model with a fixed number of ", target.trees, " trees for ", sp.name) 
	}
	gbm.object <- eval(parse(text = gbm.call))

	best.trees <- target.trees

#extract fitted values and summary table
  
	gbm.summary <- summary(gbm.object,n.trees = target.trees, plotit = FALSE)

	fits <- gbm::predict.gbm(gbm.object, x.data, n.trees = target.trees)
	if (!is.null(offset)) fits <- fits + offset
	if (family == "bernoulli") {
		fits <- exp(fits)/(1 + exp(fits))
	} else if (family == "poisson") {
		fits <- exp(fits)
	}
	fitted.values <- fits

	y_i <- y.data
	u_i <- fitted.values
	resid.deviance <- calc.deviance(y_i, u_i, weights = site.weights, family = family, calc.mean = FALSE)

	self.cor <- cor(y_i,u_i)
	self.calibration <- 0.0
	self.roc <- 0.0

	if (family == "bernoulli") {  # do this manually as we need the residuals
		deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
		residuals <- sqrt(abs(deviance.contribs * 2))
		residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)
		self.roc <- .roc(y_i,u_i)
		self.calibration <- .calibration(y_i,u_i,"binomial")
	}

	if (family == "poisson") {   # do this manually as we need the residuals
		deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
		residuals <- sqrt(abs(deviance.contribs * 2))
		residuals <- ifelse((y_i - u_i) < 0, 0 - residuals, residuals)
		self.calibration <- .calibration(y_i,u_i,"poisson")
	}

	if (family == "gaussian" | family == "laplace") {
		residuals <- y_i - u_i
	}

	mean.resid.deviance <- resid.deviance/n.cases

	z2 <- Sys.time()
	elapsed.time.minutes <- round(as.numeric(z2 - z1)/ 60, 2)  #calculate the total elapsed time

	if (verbose) {
		cat("\n")
		cat("mean total deviance =", round(mean.total.deviance,3),"\n")
		cat("mean residual deviance =", round(mean.resid.deviance,3),"\n","\n")
		cat("estimated cv deviance =", round(cv.dev,3),"; se =", round(cv.dev.se,3),"\n","\n")
		cat("training data correlation =",round(self.cor,3),"\n")
		cat("cv correlation = ",round(cv.cor,3),"; se =",round(cv.cor.se,3),"\n","\n")
		if (family == "bernoulli") {
			cat("training data AUC score =",round(self.roc,3),"\n")
			cat("cv AUC score =",round(cv.roc,3),"; se =",round(cv.roc.se,3),"\n","\n")
		}
		cat("elapsed time - ",round(elapsed.time.minutes,2),"minutes","\n")
	}

	if (n.fitted == max.trees & !silent) {
		cat("\n","########### warning ##########","\n","\n")
		cat("maximum tree limit reached - results may not be optimal","\n")
		cat("  - refit with faster learning rate or increase maximum number of trees","\n") 
	}
    
# now assemble data to be returned

	gbm.detail <- list(dataframe = data, gbm.x = gbm.x, predictor.names = names(x.data), 
    gbm.y = gbm.y, response.name = sp.name, offset = offset, family = family, tree.complexity = tree.complexity, 
    learning.rate = learning.rate, bag.fraction = bag.fraction, cv.folds = n.folds, 
    prev.stratification = prev.stratify, max.fitted = n.fitted, n.trees = target.trees, 
    best.trees = target.trees, train.fraction = 1.0, tolerance.method = tolerance.method, 
    tolerance = tolerance, var.monotone = var.monotone, date = date(), 
    elapsed.time.minutes = elapsed.time.minutes)

	training.stats <- list(null = total.deviance, mean.null = mean.total.deviance, 
    resid = resid.deviance, mean.resid = mean.resid.deviance, correlation = self.cor, 
    discrimination = self.roc, calibration = self.calibration)

	cv.stats <- list(deviance.mean = cv.dev, deviance.se = cv.dev.se, 
    correlation.mean = cv.cor, correlation.se = cv.cor.se,
    discrimination.mean = cv.roc, discrimination.se = cv.roc.se,
    calibration.mean = cv.calibration, calibration.se = cv.calibration.se)

	if (family == "bernoulli") {
		cv.stats$cv.threshold <- cv.threshold
		cv.stats$cv.threshold.se <- cv.threshold.se
	}

#  rm(x.data,y.data, envir = globalenv())           #finally, clean up the temporary dataframes

# and assemble results for return

	gbm.object$gbm.call <- gbm.detail
	gbm.object$fitted <- fitted.values
	gbm.object$fitted.vars <- fitted.vars
	gbm.object$residuals <- residuals
	gbm.object$contributions <- gbm.summary
	gbm.object$self.statistics <- training.stats
	gbm.object$cv.statistics <- cv.stats
	gbm.object$weights <- site.weights
	gbm.object$trees.fitted <- trees.fitted
	gbm.object$training.loss.values <- training.loss.values
	gbm.object$cv.values <- cv.loss.values
	gbm.object$cv.loss.ses <- cv.loss.ses
	gbm.object$cv.loss.matrix <- cv.loss.matrix
	gbm.object$cv.roc.matrix <- cv.roc.stats

	if (keep.fold.models) {
		gbm.object$fold.models <- model.list
	} else {
		gbm.object$fold.models <- NULL
	}

	if (keep.fold.vector) {
		gbm.object$fold.vector <- selector
	} else {
		gbm.object$fold.vector <- NULL
	}
	if (keep.fold.fit) {
		gbm.object$fold.fit <- fold.fit
	} else {
		gbm.object$fold.fit <- NULL
	}
	
	return(gbm.object)  
}



