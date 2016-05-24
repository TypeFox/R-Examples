
gbm.simplify <- function(
  gbm.object,                 # a gbm object describing sample intensity 
  n.folds = 10,               # number of times to repeat the analysis
  n.drops = "auto",           # can be automatic or an integer specifying the number of drops to check
  alpha = 1,                  # controls stopping when n.drops = "auto"
  prev.stratify = TRUE,       # use prevalence stratification in selecting evaluation data
  eval.data = NULL,           # an independent evaluation data set - leave here for now
  plot = TRUE )               # plot results
{
# function to simplify a brt model fitted using gbm.step
# 
# version 2.9 - J. Leathwick/J. Elith - June 2007
#
# starts with an inital cross-validated model as produced by gbm.step
# and then assesses the potential to remove predictors using k-fold cv
# does this for each fold, removing the lowest contributing predictor, 
# and repeating this process for a set number of steps 
# after the removal of each predictor, the change in predictive deviance
# is computed relative to that obtained when using all predictors
# it returns a list containing the mean change in deviance and its se
# as a function of the number of variables removed
# having completed the cross validation, it then identifies the sequence
# of variable to remove when using the full data set, testing this
# up to the number of steps used in the cross-validation phase of the analysis
# with results reported to the screen - it then returns 
# a table containing the order in which variables are to be removed
# and a list of vectors, each of which specifies the predictor col numbers
# in the original dataframe  - the latter can be used as an argument to gbm.step
# e.g., gbm.step(data = data, gbm.x = simplify.object$pred.list[[4]]...
# would implement a new analysis with the original predictor set, minus its
# four lowest contributing predictors
#

	if (! requireNamespace('gbm') ) { stop ('you need to install the gbm package to run this function') }
	requireNamespace('splines')
# first get the original analysis details..

	gbm.call <- gbm.object$gbm.call
	data <- gbm.call$dataframe
#	eval(parse(text=gbm.call$dataframe))
	n.cases <- nrow(data)
	gbm.x <- gbm.call$gbm.x
	gbm.y <- gbm.call$gbm.y
	family <- gbm.call$family
	lr <- gbm.call$learning.rate
	tc <- gbm.call$tree.complexity
	start.preds <- length(gbm.x)
	max.drops <- start.preds - 2
	response.name <- gbm.call$response.name
	predictor.names <- gbm.call$predictor.names
	n.trees <- gbm.call$best.trees
	pred.list <- list(initial = gbm.x)
	weights <- gbm.object$weights

	if (n.drops == "auto") {
		auto.stop <- TRUE
	} else {
		auto.stop <- FALSE
	}
# take a copy of the original data and starting predictors

	orig.data <- data
	orig.gbm.x <- gbm.x

#  if (!is.null(eval.data)) independent.test <- TRUE
#    else independent.test <- FALSE

# extract original performance statistics...

	original.deviance <- round(gbm.object$cv.statistics$deviance.mean,4)
	original.deviance.se <- round(gbm.object$cv.statistics$deviance.se,4)

	message("gbm.simplify - version 2.9 \nsimplifying gbm.step model for ", response.name, " with ", start.preds, " predictors and ", n.cases, " observations \noriginal deviance = ", original.deviance, "(",original.deviance.se,")")

# check that n.drops is less than n.preds - 2 and update if required

	if (auto.stop) {
		message("variable removal will proceed until average change exceeds the original se")
		n.drops <- 1 
	} else {
		if (n.drops > start.preds - 2) {
			message("value of n.drops (",n.drops,") is greater than permitted\nresetting value to ",start.preds - 2)
			n.drops <- start.preds - 2
		} else {
			message("a fixed number of ", n.drops, " drops will be tested")
		}
	}

# set up storage for results

	dev.results <- matrix(0, nrow = n.drops, ncol = n.folds)
	dimnames(dev.results) <- list(paste("drop.",1:n.drops,sep=""), paste("rep.",1:n.folds,sep=""))

	drop.count <- matrix(NA, nrow = start.preds, ncol = n.folds)
	dimnames(drop.count) <- list(predictor.names,paste("rep.",1:n.folds,sep=""))

	original.deviances <- rep(0,n.folds)

	model.list <- list(paste("model",c(1:n.folds),sep=""))     # dummy list for the tree models

# create gbm.fixed function call

	gbm.call.string <- paste("try(gbm.fixed(data=train.data,gbm.x=gbm.new.x,gbm.y=gbm.y,",sep="")
	gbm.call.string <- paste(gbm.call.string,"family=family,learning.rate=lr,tree.complexity=tc,",sep="")
	gbm.call.string <- paste(gbm.call.string,"n.trees = ",n.trees,", site.weights = weights.subset,verbose=FALSE))",sep="")
  
# now set up the fold structure

	if (prev.stratify & family == "bernoulli") {
		presence.mask <- data[,gbm.y] == 1
		absence.mask <- data[,gbm.y] == 0
		n.pres <- sum(presence.mask)
		n.abs <- sum(absence.mask)

# create a vector of randomised numbers and feed into presences
		selector <- rep(0,n.cases)
		temp <- rep(seq(1, n.folds, by = 1), length = n.pres)
		temp <- temp[order(runif(n.pres, 1, 100))]
		selector[presence.mask] <- temp

# and then do the same for absences
		temp <- rep(seq(1, n.folds, by = 1), length = n.abs)
		temp <- temp[order(runif(n.abs, 1, 100))]
		selector[absence.mask] <- temp
    }  else {  
	#otherwise make them random with respect to presence/absence
		selector <- rep(seq(1, n.folds, by = 1), length = n.cases)
		selector <- selector[order(runif(n.cases, 1, 100))]
    }

# now start by creating the intial models for each fold

	message("creating initial models...")

	gbm.new.x <- orig.gbm.x

	for (i in 1:n.folds) {        
  
# create the training and prediction folds

		train.data <- orig.data[selector!=i,]
		weights.subset <- weights[selector != i]
		eval.data <- orig.data[selector==i,]

		model.list[[i]] <- eval(parse(text=gbm.call.string))  # create a fixed size object

# now make predictions to the withheld fold

		u_i <- eval.data[,gbm.y]
		y_i <- gbm::predict.gbm(model.list[[i]], eval.data, n.trees, "response")

		original.deviances[i] <- round(calc.deviance(u_i,y_i, family = family, calc.mean = TRUE),4)
 
	} # end of creating initial models

	n.steps <- 1
 
	message("dropping predictor:", appendLF = FALSE)
	
	while (n.steps <= n.drops & n.steps <= max.drops) {
	
		message(" ", n.steps, appendLF = FALSE)

		for (i in 1:n.folds) {

# get the right data

			train.data <- orig.data[selector!=i,]
			eval.data <- orig.data[selector==i,]
			weights.subset <- weights[selector != i]

# get the current model details

			gbm.x <- model.list[[i]]$gbm.call$gbm.x
			n.preds <- length(gbm.x)
			these.pred.names <- model.list[[i]]$gbm.call$predictor.names
			contributions <- model.list[[i]]$contributions

# get the index number in pred.names of the last variable in the contribution table

			last.variable <- match(as.character(contributions[n.preds,1]),these.pred.names)
			gbm.new.x <- gbm.x[-last.variable]

# and keep a record of what has been dropped

			last.variable <- match(as.character(contributions[n.preds,1]),predictor.names)
			drop.count[last.variable,i] <- n.steps

			model.list[[i]] <- eval(parse(text=gbm.call.string))  # create a fixed size object

			u_i <- eval.data[,gbm.y]
			y_i <- gbm::predict.gbm(model.list[[i]],eval.data,n.trees,"response")

			deviance <- round(calc.deviance(u_i,y_i, family = family, calc.mean = TRUE),4)
    
# calculate difference between intial and new model by subtracting new from old because we want to minimise deviance

			dev.results[n.steps,i] <- round(deviance - original.deviances[i] ,4)
		}

		if (auto.stop){ # check to see if delta mean is less than original deviance error estimate

			delta.mean <- mean(dev.results[n.steps,])

			if (delta.mean < (alpha * original.deviance.se)) {
				n.drops <- n.drops + 1
				dev.results <- rbind(dev.results, rep(0,n.folds))
			}
		}   
		n.steps <- n.steps + 1
	}
	message("")

# now label the deviance matrix

	dimnames(dev.results) <- list(paste("drop.",1:n.drops,sep=""), paste("rep.", 1:n.folds, sep=""))

# calculate mean changes in deviance and their se

	mean.delta <- apply(dev.results,1,mean)
	se.delta <- sqrt(apply(dev.results,1,var))/sqrt(n.folds)

###########################

	if (plot) {
		y.max <- 1.5 * max(mean.delta + se.delta)
		y.min <- 1.5 * min(mean.delta - se.delta)
		plot(seq(0,n.drops),c(0,mean.delta),xlab="variables removed", ylab = "change in predictive deviance",type='l',ylim=c(y.min,y.max))
		lines(seq(0,n.drops),c(0,mean.delta) + c(0,se.delta),lty = 2)
		lines(seq(0,n.drops),c(0,mean.delta) - c(0,se.delta),lty = 2)
		abline(h = 0 , lty = 2, col = 3)
		min.y <- min(c(0,mean.delta))
		min.pos <- match(min.y,c(0,mean.delta)) - 1 # subtract one because now zero base
		abline(v = min.pos, lty = 3, col = 2)
		abline(h = original.deviance.se, lty = 2, col = 2)
		title(paste("RFE deviance - ",response.name," - folds = ",n.folds,sep=""))
	}

# and do a final backwards drop sequence from the original model

	message("processing final dropping of variables with full data")

	gbm.call.string <- paste("try(gbm.fixed(data=orig.data,gbm.x=gbm.new.x,gbm.y=gbm.y,",sep="")
	gbm.call.string <- paste(gbm.call.string,"family=family,learning.rate=lr,tree.complexity=tc,",sep="")
	gbm.call.string <- paste(gbm.call.string,"n.trees = ",n.trees,", site.weights = weights,verbose=FALSE))",sep="")

	n.steps <- n.steps - 1 #decrement by one to reverse last increment in prev loop

	final.model <- gbm.object  # restore the original model and data
	train.data <- orig.data

# and set up storage

	final.drops <- matrix(NA, nrow = start.preds, ncol = 1)
	dimnames(final.drops) <- list(predictor.names,"step")

	for (i in 1:n.steps) {

# get the current model details

		gbm.x <- final.model$gbm.call$gbm.x
		n.preds <- length(gbm.x)
		these.pred.names <- final.model$gbm.call$predictor.names
		contributions <- final.model$contributions

		message(i, "-", as.character(contributions[n.preds,1]))

# get the index number in pred.names of the last variable in the contribution table

		last.variable <- match(as.character(contributions[n.preds,1]),these.pred.names)
		gbm.new.x <- gbm.x[-last.variable]

# and keep a record of what has been dropped

		last.variable <- match(as.character(contributions[n.preds,1]),predictor.names)
		final.drops[last.variable] <- i

		final.model <- eval(parse(text=gbm.call.string))  # create a fixed size object

	}

#and then the corresponding numbers

	removal.list <- dimnames(final.drops)[[1]]
	removal.list <- removal.list[order(final.drops)]
	removal.list <- removal.list[1:n.drops]

	removal.numbers <- rep(0,n.steps)

# construct predictor lists to faciliate final model fitting

	for (i in 1:n.steps) {
		removal.numbers[i] <- match(removal.list[i],predictor.names)
		pred.list[[i]] <- orig.gbm.x[0-removal.numbers[1:i]]
		names(pred.list)[i] <- paste("preds.",i,sep="")
	}

	deviance.summary <- data.frame(mean = round(mean.delta,4), se = round(se.delta,4))

	final.drops <- data.frame("preds" = dimnames(final.drops)[[1]][order(final.drops)], "order" = final.drops[order(final.drops)])

	return(list(deviance.summary = deviance.summary, deviance.matrix = dev.results, drop.count = drop.count, final.drops = final.drops, pred.list = pred.list, gbm.call = gbm.call))

}
