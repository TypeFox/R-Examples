## Copyright (c) 2004-2012, Ph. Grosjean <phgrosjean@sciviews.org>
##
## This file is part of ZooImage
##
## ZooImage is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## ZooImage is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with ZooImage.  If not, see <http://www.gnu.org/licenses/>.

## Create basically a mlearning object, but with predicted and cvpredicted added
## to it, and the '+other+' level added at the end of all levels
ZIClass <- function (formula, data, method = getOption("ZI.mlearning",
"mlRforest"), calc.vars = getOption("ZI.calcVars", calcVars), drop.vars = NULL,
drop.vars.def = dropVars(), cv.k = 10, cv.strat = TRUE, ...,
subset, na.action = na.omit)
{	
	## Rework calc.vars to freeze the list of variables to drop
	## Default (minimal) calc.vars function if none is provided
	if (!length(calc.vars)) {
		## Create a simple calc.vars function that just drop variables
		calc.vars <- function(x, drop.vars = NULL, drop.vars.def = dropVars()) {
			## Eliminate variables that are not predictors... and use Id as rownames
			Id <- x$Id
			if (length(Id)) rownames(x) <- Id
	
			## Variables to drop
			dropAll <- unique(as.character(c(drop.vars, drop.vars.def)))
			for (dropVar in dropAll) x[[dropVar]] <- NULL
	
			## Return the recalculated data frame
			x
		}
	}
	
	## Freeze data for drop.vars and drop.vars.def arguments of calc.vars
	if (!length(drop.vars)) drop.vars <- character(0) else
		drop.vars <- as.character(drop.vars)
	formals(calc.vars)$drop.vars <- drop.vars
	if (!length(drop.vars.def)) drop.vars.def <- character(0) else
		drop.vars.def <- as.character(drop.vars.def)
	formals(calc.vars)$drop.vars.def <- drop.vars.def
	
	## Check calc.vars and use it on data
	if (!is.function(calc.vars)) {
		stop("'calc.vars' must be a function or NULL")
	} else data <- calc.vars(data)
	
	## Train the machine learning algorithm
	ZI.class <- mlearning(formula, data = data, method = method,
		model.args = list(formula  = formula, data = substitute(data),
		subset = substitute(subset)), call = match.call(), ...,
		subset = subset, na.action = substitute(na.action))
		
	## Add ZIClass as class of the object
	class(ZI.class) <- c("ZIClass", class(ZI.class))
	
	## Save our customized calc.vars function in the object
	attr(ZI.class, "calc.vars") <- calc.vars

	## Get useful attributes from ZITrain
	attr(ZI.class, "traindir") <- attr(data, "traindir")
	attr(ZI.class, "path") <- attr(data, "path")

	## Calculate predictions with full training set
    attr(ZI.class, "predict") <- predict(ZI.class, data, calc = FALSE)

	## Possibly make a k-fold cross-validation and check results
	if (length(cv.k)) {
		attr(ZI.class, "cvpredict") <- cvpredict(ZI.class, type = "class",
			cv.k = cv.k, cv.strat = cv.strat)
		attr(ZI.class, "k") <- cv.k
		attr(ZI.class, "strat") <- cv.strat
	}
	
	ZI.class
}

print.ZIClass <- function (x, ...)
{
	algorithm <- attr(x, "algorithm")
	classes <- attr(x, "response")
	lclasses <- levels(classes)
    predicted <- attr(x, "predict")
	if (is.list(predicted)) predicted <- predicted$class
	k <- attr(x, "k")
	strat <- attr(x, "strat")
	cat("A 'ZIClass' object predicting for", length(lclasses), "classes:\n")
	print(lclasses)
	Confu <- table(classes, predicted)
	SelfConsist <- 100 * (sum(diag(Confu)) / sum(Confu))

	## Change the number of digits to display
	oldDigits <- options(digits = 4)
	on.exit(options(oldDigits))
	cat("\nAlgorithm used:", algorithm, "\n")
	cat("Self-consistency: ", SelfConsist, "%\n", sep = "")
	if (!is.null(k)) {
		if (isTRUE(strat)) msg <- ", stratified" else msg <- ""
    	cat("K-fold cross validation error estimation (k = ", k, msg, "):\n",
			sep = "")
		cvpredicted <- attr(x, "cvpredict")
		if (is.list(cvpredicted)) cvpredicted <- cvpredicted$class
		prior <- table(classes)
		ok <- diag(table(classes, cvpredicted))
		err <- 100 * (1 - (sum(ok) / sum(prior)))
		cat("Error rate: ", err, "%\n", sep = "")
		cat("\nError per class:\n")
		`Error (%)` <- sort(1 - (ok / prior)) * 100
		print(as.data.frame(`Error (%)`))
	}
	invisible(x)
}

summary.ZIClass <- function(object, sort.by = "Fscore", decreasing = TRUE,
na.rm = FALSE, ...)
{
	## Get the confusion object out of a ZIClass object and calc stats from there
	summary(confusion(object, response(object)), sort.by = sort.by, decreasing = decreasing,
		na.rm = na.rm, ...)
}

predict.ZIClass <- function (object, newdata, calc = TRUE, class.only = TRUE,
type = "class", ...)
{
	## Make sure we have correct objects
	if (!inherits(object, "ZIClass"))
		stop("'object' must be a 'ZIClass' object")
	
	if (!missing(newdata)) {
		if (!inherits(newdata, c("ZIDat", "data.frame")))
			stop("'newdata' must be a 'ZIDat' or 'data.frame' object")
		data <- as.data.frame(newdata)
		if (isTRUE(as.logical(calc)))
			data <- attr(object, "calc.vars")(data)
	}
    
	class(object) <- class(object)[-1]
		
	class.only <- isTRUE(as.logical(class.only))
	
	type <- as.character(type)[1]
	
	if (class.only && type != "class") {
		warning("with class.only == TRUE, type can only be 'class' and is force to it")
		type <- "class"
	}
	
	## Perform the prediction
	if (missing(newdata)) {
		res <- predict(object, ...)
	} else res <- predict(object, newdata = data, ...)
	
	## Return either the prediction, or the ZIDat object with Predicted
	## column append/replaced
	if (class.only) res else {
		newdata$Predicted <- res
		newdata
	}
}

confusion.ZIClass <- function (x, y = response(x),
labels = c("Actual", "Predicted"), useNA = "ifany", prior, use.cv = TRUE, ...) {
	## Check labels
	labels <- as.character(labels)
	if (length(labels) != 2)
		stop("You must provide exactly 2 character strings for 'labels'")
	
	## Extract class2: cvpredict or predict from the object
	if (isTRUE(as.logical(use.cv))) {
		class2 <- attr(x, "cvpredict")
		if (is.list(class2)) class2 <- class2$class
		if (is.null(class2))
			stop("No or wrong cross-validated predictions in this ZIClass object")
	} else { # Use predict
		class2 <- attr(x, "predict")
		if (is.list(class2)) class2 <- class2$class
	}
	
	## Check that both variables are of same length and same levels
	if (length(y) != length(class2))
		stop("lengths of 'x' and 'y' are not the same")
	
	## Full list of levels is in (cv)predict in class2...
	## Response in y may have dropped levels! 
	lev1 <- levels(y)
	lev2 <- levels(class2)
	if (!all(lev1  %in% lev2))
		stop("levels of 'x' and 'y' do not match")
	
	## Rework levels in y to make sure they match perfectly thos in class2
	y <- factor(as.character(y), levels = lev2)
	
	## Construct the confusion object
	if (missing(prior)) {
		getNamespace("mlearning")$.confusion(data.frame(class1 = y,
			class2 = class2), labels = labels, useNA = useNA, ...)
	} else {
		getNamespace("mlearning")$.confusion(data.frame(class1 = y,
			class2 = class2), labels = labels, useNA = useNA, prior = prior, ...)
	}
}
