response <- function (object, ...)
	UseMethod("response")
	
response.default <- function (object, ...)
	attr(object, "response")
	
train <- function (object, ...)
	UseMethod("train")
	
train.default <- function (object, ...)
	attr(object, "train")

## TODO: test performances of optimized code for Class ~ .
mlearning <- function (formula, data, method, model.args, call = match.call(),
..., subset, na.action = na.fail)
{
	## Our own construction of response vector and terms matrix
	if (missing(model.args))
		model.args <- list(formula  = formula, data = substitute(data),
			subset = substitute(subset))
	
	## Get data and initial number of cases
	data <- eval.parent(model.args$data)
	nobs <- NROW(data)
	
	## Special case for formula like response ~ . which speeds up calc and
	## uses less memory than model.frame()
	isSimpleFormula <- function (formula) {
		vars <- all.vars(formula)
		(length(vars) == 2 && vars[2] == ".") || # Supervised (response ~ .)
		(length(vars) == 1 && vars[1] == ".")	 # Unsupervised (~ .)
	}	
	optim <- isSimpleFormula(model.args$formula)
	if (optim) {
		## data do not need to be changed... except for subset or na.action
		if (model.args$subset != "")
			data <- data[eval.parent(model.args$subset), ] 
		if (missing(na.action) || as.character(na.action) == "") {
			## Use same rules as model.frame():
			## (1) any na.action attribute of data
			na.action <- attr(data, "na.action")
			## (2) option na.action, or (3) na.fail
			if (is.null(na.action))
				na.action <- getOption("na.action", na.fail)
		}
		## Apply provided na.action
		data <- match.fun(na.action)(data)
		if (is.function(na.action)) na.action <- substitute(na.action)
		na.action <- as.character(na.action)
		model.terms <- terms(formula, data = data[1, ])
		attr(data, "terms") <- model.terms
	} else { # Use model.frame()
		if (missing(na.action) || as.character(na.action) == "") {
			data <- do.call("model.frame", model.args)
			na.action <- as.character(attr(data, "na.action"))
			if (!length(na.action)) {
				na.action <- "na.pass" # If not provided, either pass, or no NAs!
			} else na.action <- paste("na", class(na.action), sep = ".")
		} else {
			model.args$na.action <- na.action
			data <- do.call("model.frame", model.args)
			if (is.function(na.action)) na.action <- substitute(na.action)
			na.action <- as.character(na.action)
		}
		model.terms <- attr(data, "terms")
	}
	## Final number of observations
	nobs[2] <- NROW(data)
	names(nobs) <- c("initial", "final")
	
	## Construct the matrix of numeric predictors and the response
	term.labels <- attr(model.terms, "term.labels")
	response.pos <- attr(model.terms, "response")
	if (!response.pos) {
		response.label <- NULL
		train <- data
		response <- NULL
		lev <- NULL
		type <- "unsupervised"
	} else { # Supervised classification or regression
		response.label <- deparse(attr(model.terms, "variables")
			[[response.pos + 1]])
		response <- data[, response.label]
		if (is.factor(response)) {
			lev <- levels(response)
			response <- droplevels(response)
			type <- "classification"
		} else {
			if (!is.numeric(response))
				stop("response variable must be factor or numeric")
			lev <- NULL
			type <- "regression"
		}
		train <- data[, term.labels]
	}
	
	## Calculate weights
	w <- model.weights(data)
	if (length(w) == 0L) w <- rep(1, nrow(train))
	
	## Pass special arguments to the default method
	args <- list()
	args$formula <- formula
	args$levels <- lev
	args$n <- nobs
	args$weights <- w
	args$optim <- optim
	args$type <- type
	args$na.action <- substitute(na.action)
	args$mlearning.call <- call
	args$method <- method
	
	## Construct the mlearning object
	match.fun(method)(train = train, response = response, .args. = args, ...)
}

print.mlearning <- function (x, ...)
{
	cat("A mlearning object of class ", class(x)[1], " (",
		attr(x, "algorithm"), "):\n", sep = "")
	type <- attr(x, "type")
	switch(type,
		regression = cat("[regression variant]\n"),
		unsupervised = cat("[unsupervised classification variant]\n"))
	cat("Call: ", deparse(attr(x, "mlearning.call")),
		"\n", sep = "")
	
	if (type == "classification") {
		## Number of cases used
		n <- attr(x, "n")
		if (n["final"] < n["initial"]) {
			msg <- paste("Trained using", n["final"], "out of",
				n["initial"], "cases:")
		} else msg <- paste("Trained using", n["final"], "cases:")
	
		## Categories
		classes <- attr(x, "response")
		levUsed <- levels(classes)
		levIni <- levels(x)
		if (length(levUsed) < length(levIni)) {
			cat("Levels with no cases in the training set that were eliminated:\n")
			print(levIni[!levIni %in% levUsed])
		}

		## Number of cases per used categories
		print(table(classes, dnn = msg))
	}
	return(invisible(x))
}

summary.mlearning <- function (object, ...)
{
	train <- attr(object, "train")
	response <- attr(object, "response")
	mlearning.class <- class(object)[1]
	class(object) <- class(object)[-(1:2)]
	## Summary is sometimes implemented as print() for some machine
	## learning algorithms... this is is 'summary' attribute
	sumfun <- attr(object, "summary")
	if (length(sumfun)) {
		res <- match.fun(sumfun)(object, ...)
	} else res <- object
	class(res) <- c("summary.mlearning", class(res))
	attr(res, "mlearning.class") <- mlearning.class
	attr(res, "algorithm") <- attr(object, "algorithm")
	attr(res, "type") <- attr(object, "type")
	attr(res, "mlearning.call") <- attr(object, "mlearning.call")
	res
}
		
print.summary.mlearning <- function (x, ...)
{
	cat("A mlearning object of class ", attr(x, "mlearning.class"), " (",
		attr(x, "algorithm"), "):\n", sep = "")
	type <- attr(x, "type")
	switch(type,
		regression = cat("[regression variant]\n"),
		unsupervised = cat("[unsupervised classification variant]\n"))
	cat("Initial call: ", deparse(attr(x, "mlearning.call")),
		"\n", sep = "")
	X <- x
	class(X) <- class(x)[-1]
	print(X)
	invisible(x)
}

plot.mlearning <- function (x, y, ...)
{
	train <- attr(x, "train")
	response <- attr(x, "response")
	class(x) <- class(x)[-(1:2)]
	plot(x, ...)
}

.membership <- function (x, levels, scale = TRUE)
{
	## Make sure x is a matrix of numerics
	x <- as.matrix(x)
	if (!is.numeric(x))
		stop("'x' must be numeric")
	
	## Make sure all columns are named with names in levels
	nms <- colnames(x)
	if (!length(nms))
		stop("missing column names in 'x'")
	if (any(!nms %in% levels))
		stop("One or more column in 'x' not in 'levels'")
	
	## Add columns of zeros for inexistant levels
	toAdd <- levels[!levels %in% nms]
	if (length(toAdd)) {
		xAdd <- matrix(0, nrow = NROW(x), ncol = length(toAdd))
		colnames(xAdd) <- toAdd
		x <- cbind(x, xAdd)
	}
	
	## Make sure columns are in the same order as levels
	x <- x[, levels]
	
	## Possibly scale to one, row-wise
	if (isTRUE(as.logical(scale)))
		x <- x / apply(x, 1, sum)
	
	x
}

.expandFactor <- function (f, n, ndrop)
{
	if (!length(ndrop) || class(ndrop) != "exclude") return(f)
	res <- factor(rep(NA, n), levels = levels(f))
	res[-ndrop] <- f
	res
}
	
.expandMatrix <- function (m, n, ndrop)
{
	if (!length(ndrop) || class(ndrop) != "exclude") return(m)
	res <- matrix(NA, nrow = n, ncol = ncol(m))
	res[-ndrop, ] <- m
	res
}

predict.mlearning <- function(object, newdata,
type = c("class", "membership", "both"), method = c("direct", "cv"),
na.action = na.exclude, ...)
{
	## Not usable for unsupervised type
	if (attr(object, "type") == "unsupervised")
		stop("no predict() method for unsupervised version")
	
	## If method == "cv", delegate to cvpredict()
	if (as.character(method)[1] == "cv") {
		if (!missing(newdata))
			stop("cannot handle new data with method = 'cv'")
		return(cvpredict(object = object, type = type, ...))
	}
	
	## Recalculate newdata according to formula...
	if (missing(newdata)) { # Use train
		newdata <- attr(object, "train")
	} else if (attr(object, "optim")) { # Use optimized approach
		## Just keep vars similar as in train
		vars <- names(attr(object, "train"))
		if (!all(vars %in% names(newdata)))
			stop("one or more missing variables in newdata")
		newdata <- newdata[, vars]
	} else { # Use model.frame
		newdata <- model.frame(formula = attr(object, "formula"),
			data = newdata, na.action = na.pass)[, names(attr(object, "train"))]
	}
	## Do we need only numerical predictors
	if (attr(object, "numeric.only"))
		newdata <- sapply(as.data.frame(newdata), as.numeric)
	
	## Determine how many data and perform na.action
	n <- NROW(newdata)
	newdata <- match.fun(na.action)(newdata)
	ndrop <- attr(newdata, "na.action")
		
	## Delegate to the original predict() method
	class(object) <- class(object)[-(1:2)]
	if (attr(object, "type") == "regression")
		return(predict(object, newdata = newdata, ...))
	
	## Otherwise, this is a supervised classification
	type <- as.character(type)[1]
	## Special case for both
	if (type == "both") type <- c("class", "membership")
	## Check that type is supported and look for corresponding type name
	## in original predict() method
	pred.type <- attr(object, "pred.type")
	if (!all(type %in% names(pred.type)))
		stop("unsupported predict type")
	
	if (length(type) == 2) {
		## Special case where we predict both class and membership
		classes <- predict(object, newdata = newdata,
			type = pred.type["class"], ...)
		members <- predict(object, newdata = newdata,
			type = pred.type["membership"], ...)
		## Create a list with both res
		levels <- levels(object)
		return(list(class = .expandFactor(factor(as.character(classes),
			levels = levels), n, ndrop),
			membership = .expandMatrix(.membership(members, levels = levels),
			n, ndrop)))
	} else {
		res <- predict(object, newdata = newdata, type = pred.type[type], ...)
	}
	
	## Rework result according to initial levels (before drop of empty ones)
	res <- switch(type,
		class = .expandFactor(factor(as.character(res), levels = levels(object)),
			n, ndrop),
		membership = .expandMatrix(.membership(res, levels = levels(object)),
			n, ndrop),
		switch(class(res)[1],
			factor = .expandFactor(res, n, ndrop),
			matrix = .expandMatrix(res, n, ndrop),
			res))
		
	res
}

cvpredict <- function (object, ...)
	UseMethod("cvpredict")

cvpredict.mlearning <- function(object, type = c("class", "membership", "both"),
cv.k = 10, cv.strat = TRUE, ...)
{
	type <- switch(attr(object, "type"),
		regression = "class", # Another way to ignore 'type' for regressions
		classification = as.character(type)[1],
		stop("works only for classification or regression mlearning objects"))
	
	if (type == "class") {
		predictions <- TRUE
		getmodels <- FALSE
	} else if (type == "membership") {
		predictions <- FALSE
		getmodels <- TRUE
	} else if (type == "both") {
		predictions <- TRUE
		getmodels <- TRUE
	} else stop("type must be 'class', 'membership' or 'both'")
	
	## Create data, using numbers are rownames
	data <- data.frame(.response. = response(object), train(object))
	rn <- rownames(data)
	rownames(data) <- 1:NROW(data)
	
	## The predict() method with ... arguments added to the call
	constructPredict <- function (...) {
		fun <- function (object, newdata) return()
		body(fun) <- as.call(c(list(substitute(predict),
			object = substitute(object), newdata = substitute(newdata)), list(...)))
		fun
	}
	Predict <- constructPredict(...)
	
	## Perform cross-validation for prediction
	args <- attr(object, "args")
	if (!is.list(args)) args <- list()
	args$formula <- substitute(.response. ~ .)
	args$data <- substitute(data)
	args$model <- substitute(mlearning)
	args$method <- attr(object, "method")
	args$predict <- substitute(Predict)
	args$estimator <- "cv"
	args$est.para <- control.errorest(predictions = predictions,
		getmodels = getmodels, k = cv.k, strat = cv.strat)
	est <- do.call(errorest, args)
	
	## Only class
	if (type == "class") {
		res <- est$predictions
	} else {
		## Need to calculate membership
		predCV <- function (x, object, ...) {
			Train <- train(object)
			rownames(Train) <- 1:NROW(Train)
			suppressWarnings(predict(x, newdata =
				Train[-as.numeric(rownames(train(x))), ], ...))
		}
	
		## Apply predict on all model and collect results together
		membership <- lapply(est$models, predCV, object = object,
			type = "membership", na.action = na.exclude, ...)
	
		## Concatenate results
		membership <- do.call(rbind, membership)
	
		## Sort in correct order and replace initial rownames
		ord <- as.numeric(rownames(membership))
		## Sometimes, errorest() duplicates one or two items in two models
		## (rounding errors?) => eliminate them here
		notDup <- !duplicated(ord)
		membership <- membership[notDup, ]
		ord <- ord[notDup]
		
		# Restore order of the items
		rownames(membership) <- rn[ord]
		pos <- order(ord)
		membership <- membership[pos, ]
	
		if (type == "membership") {
			res <- membership
		} else {  # Need both class and membership
			## Because we don't know who is who in est$predictions in case of
			## duplicated items in est$models, we prefer to recalculate classes
			classes <- unlist(lapply(est$models, predCV, object = object,
				type = "class", na.action = na.exclude, ...))
			classes <- classes[notDup]
			classes <- classes[pos]
			
			## Check that both classes levels are the same!
			if (any(levels(classes) != levels(est$predictions)))
				warning("cross-validated classes do not match")

			res <- list(class = classes, membership = membership)
		}
	}
	
	## Add est object as "method" attribute, without predictions or models
	est$name <- "cross-validation"
	est$predictions <- NULL
	est$models <- NULL
	est$call <- match.call()
	est$strat <- cv.strat
	attr(res, "method") <- est
	
	res
}

## Note: ldahist() in MASS (when only one LD) seems to be broken!
mlLda <- function (...)
	UseMethod("mlLda")

mlLda.formula <- function (formula, data, ..., subset, na.action)
	mlearning(formula, data = data, method = "mlLda", model.args =
		list(formula  = formula, data = substitute(data),
		subset = substitute(subset)), call = match.call(), ...,
		subset = subset, na.action = substitute(na.action))

mlLda.default <- function (train, response, ...)
{
	if (!is.factor(response))
		stop("only factor response (classification) accepted for mlLda")

	dots <- list(...)
	.args. <- dots$.args.
	dots$.args. <- NULL
	if (!length(.args.)) .args. <- list(levels = levels(response),
		n = c(intial = NROW(train), final = NROW(train)),
		type = "classification", na.action = "na.pass",
		mlearning.call = match.call(), method = "mlLda")
	
	## Check if there are factor predictors
	if (any(sapply(train, is.factor)))
		warning("force conversion from factor to numeric; may be not optimal or suitable")
	
	## Return a mlearning object
	structure(MASS:::lda.default(x = sapply(train, as.numeric),
		grouping = response, ...), formula = .args.$formula, train = train,
		response = response, levels = .args.$levels, n = .args.$n, args = dots,
		optim = .args.$optim, numeric.only = TRUE, type = .args.$type,
		pred.type = c(class = "class", membership = "posterior", projection = "x"),
		summary = NULL, na.action = .args.$na.action,
		mlearning.call = .args.$mlearning.call, method = .args.$method,
		algorithm = "linear discriminant analysis",
		class = c("mlLda", "mlearning", "lda"))
}

predict.mlLda <- function(object, newdata,
type = c("class", "membership", "both", "projection"), prior = object$prior,
dimension, method = c("plug-in", "predictive", "debiased", "cv"), ...)
{
	if (!inherits(object, "mlLda"))
		stop("'object' must be a 'mlLda' object")
	
	## If method == "cv", delegate to cvpredict()
	method <- as.character(method)[1]
	if (method == "cv") {
		if (!missing(newdata))
			stop("cannot handle new data with method = 'cv'")
		if (missing(dimension)) {
			return(cvpredict(object = object, type = type, prior = prior, ...))	
		} else {
			return(cvpredict(object = object, type = type, prior = prior,
				dimension = dimension, ...))
		}
	}
	
	## Recalculate newdata according to formula...
	if (missing(newdata)) { # Use train
		newdata <- attr(object, "train")
	} else if (attr(object, "optim")) { # Use optimized approach
		## Just keep vars similar as in train
		vars <- names(attr(object, "train"))
		if (!all(vars %in% names(newdata)))
			stop("One or more missing variables in newdata")
		newdata <- newdata[, vars]
	} else { # Use model.frame
		newdata <- model.frame(formula = attr(object, "formula"),
			data = newdata, na.action = na.pass)[, names(attr(object, "train"))]
	}
	## Only numerical predictors
	newdata <- sapply(as.data.frame(newdata), as.numeric)
	
	## dimension
	if (missing(dimension)) {
        dimension <- length(object$svd)
	} else {
		dimension <- min(dimension, length(object$svd))
	}
	
	## Delegate to the MASS predict.lda method
	class(object) <- class(object)[-(1:2)]
	## I need to suppress warnings, because NAs produce ennoying warnings!
	res <- suppressWarnings(predict(object, newdata = newdata, prior = prior,
		dimen = dimension, method = method, ...))
	
	## Rework results according to what we want
	switch(as.character(type)[1],
		class = factor(as.character(res$class), levels = levels(object)),
		membership = .membership(res$posterior, levels = levels(object)),
		both = list(class = factor(as.character(res$class),
			levels = levels(object)), membership = .membership(res$posterior,
			levels = levels(object))),
		projection = res$x,
		stop("unrecognized 'type' (must be 'class', 'membership', 'both' or 'projection')"))
}

mlQda <- function (...)
	UseMethod("mlQda")

mlQda.formula <- function (formula, data, ..., subset, na.action)
	mlearning(formula, data = data, method = "mlQda", model.args =
		list(formula  = formula, data = substitute(data),
		subset = substitute(subset)), call = match.call(), ...,
		subset = subset, na.action = substitute(na.action))

mlQda.default <- function (train, response, ...)
{
	if (!is.factor(response))
		stop("only factor response (classification) accepted for mlQda")

	dots <- list(...)
	.args. <- dots$.args.
	dots$.args. <- NULL
	if (!length(.args.)) .args. <- list(levels = levels(response),
		n = c(intial = NROW(train), final = NROW(train)),
		type = "classification", na.action = "na.pass",
		mlearning.call = match.call(), method = "mlQda")
	
	## Check if there are factor predictors
	if (any(sapply(train, is.factor)))
		warning("force conversion from factor to numeric; may be not optimal or suitable")
	
	## Return a mlearning object
	structure(MASS:::qda.default(x = sapply(train, as.numeric),
		grouping = response, ...), formula = .args.$formula, train = train,
		response = response, levels = .args.$levels, n = .args.$n, args = dots,
		optim = .args.$optim, numeric.only = TRUE, type = .args.$type,
		pred.type = c(class = "class", membership = "posterior"),
		summary = NULL, na.action = .args.$na.action,
		mlearning.call = .args.$mlearning.call, method = .args.$method,
		algorithm = "quadratic discriminant analysis",
		class = c("mlQda", "mlearning", "qda"))
}

predict.mlQda <- function(object, newdata, type = c("class", "membership", "both"),
prior = object$prior, method = c("plug-in", "predictive", "debiased", "looCV",
"cv"), ...)
{
	if (!inherits(object, "mlQda"))
		stop("'object' must be a 'mlQda' object")
	
	## If method == "cv", delegate to cvpredict()
	method <- as.character(method)[1]
	if (method == "cv") {
		if (!missing(newdata))
			stop("cannot handle new data with method = 'cv'")
		return(cvpredict(object = object, type = type, prior = prior, ...))
	}
	
	## Recalculate newdata according to formula...
	if (missing(newdata)) { # Use train
		newdata <- attr(object, "train")
	} else if (attr(object, "optim")) { # Use optimized approach
		## Just keep vars similar as in train
		vars <- names(attr(object, "train"))
		if (!all(vars %in% names(newdata)))
			stop("One or more missing variables in newdata")
		newdata <- newdata[, vars]
	} else { # Use model.frame
		newdata <- model.frame(formula = attr(object, "formula"),
			data = newdata, na.action = na.pass)[, names(attr(object, "train"))]
	}
	## Only numerical predictors
	newdata <- sapply(as.data.frame(newdata), as.numeric)
		
	## Delegate to the MASS predict.qda method
	class(object) <- class(object)[-(1:2)]
	## I need to suppress warnings, because NAs produce ennoying warnings!
	res <- suppressWarnings(predict(object, newdata = newdata, prior = prior,
		method = method, ...))
	
	## Rework results according to what we want
	switch(as.character(type)[1],
		class = factor(as.character(res$class), levels = levels(object)),
		membership = .membership(res$posterior, levels = levels(object)),
		both = list(class = factor(as.character(res$class),
			levels = levels(object)), membership = .membership(res$posterior,
			levels = levels(object))),
		stop("unrecognized 'type' (must be 'class', 'membership' or 'both')"))
}

mlRforest <- function (...)
	UseMethod("mlRforest")

mlRforest.formula <- function (formula, data, ntree = 500, mtry,
replace = TRUE, classwt = NULL, ..., subset, na.action)
{
	if (missing(mtry)) {
		mlearning(formula, data = data, method = "mlRforest", model.args =
			list(formula  = formula, data = substitute(data),
			subset = substitute(subset)), call = match.call(), ntree = ntree,
			replace = replace, classwt = classwt, ...,
			subset = subset, na.action = substitute(na.action))	
	} else {
		mlearning(formula, data = data, method = "mlRforest", model.args =
			list(formula  = formula, data = substitute(data),
			subset = substitute(subset)), call = match.call(), ntree = ntree,
			mtry = mtry, replace = replace, classwt = classwt, ...,
			subset = subset, na.action = substitute(na.action))	
	}
}

mlRforest.default <- function (train, response, ntree = 500, mtry,
replace = TRUE, classwt = NULL, ...)
{
	dots <- list(...)
	.args. <- dots$.args.
	dots$.args. <- NULL
	if (!length(.args.)) {
		if (!length(response)) {
			type <- "unsupervised"
		} else if (is.factor(response)) {
			type <- "classification"
		} else type <- "regression"
		.args. <- list(levels = levels(response),
		n = c(intial = NROW(train), final = NROW(train)),
		type = type, na.action = "na.pass",
		mlearning.call = match.call(), method = "mlRforest")
	}
	dots$ntree <- ntree
	dots$replace <- replace
	dots$classwt <- classwt
	
	## Return a mlearning object
	if (missing(mtry) || !length(mtry)) {
		res <- randomForest:::randomForest.default(x = train,
		y = response, ntree = ntree, replace = replace,
		classwt = classwt, ...)
	} else {
		dots$mtry <- mtry
		res <- randomForest:::randomForest.default(x = train,
		y = response, ntree = ntree, mtry = mtry, replace = replace,
		classwt = classwt, ...)
	}
	 
	structure(res, formula = .args.$formula, train = train,
		response = response, levels = .args.$levels, n = .args.$n, args = dots,
		optim = .args.$optim, numeric.only = FALSE, type = .args.$type,
		pred.type = c(class = "response", membership = "prob", vote ="vote"),
		summary = NULL, na.action = .args.$na.action,
		mlearning.call = .args.$mlearning.call, method = .args.$method,
		algorithm = "random forest",
		class = c("mlRforest", "mlearning", "randomForest"))
}

predict.mlRforest <- function(object, newdata,
type = c("class", "membership", "both", "vote"), method = c("direct", "oob", "cv"),
...)
{
	type <- as.character(type)[1]
	
	## If method == "cv", delegate to cvpredict()
	method <- as.character(method)[1]
	if (method == "cv") {
		if (!missing(newdata))
			stop("cannot handle new data with method = 'cv'")
		return(cvpredict(object = object, type = type, ...))
	} else if (method == "oob") { # Get out-of-bag prediction!
		if (!missing(newdata))
			stop("you cannot provide newdata with method = 'oob'")
		
		toProps <- function (x, ntree) {
			if (sum(x[1, ] > 1)) {
				res <- t(apply(x, 1, "/", ntree))
			} else res <- x
			class(res) <- "matrix"
			res
		}

		toVotes <- function (x, ntree) {
			if (sum(x[1, ] < ntree - 1)) {
				res <- round(t(apply(x, 1, "*", ntree)))
			} else res <- x
			class(res) <- "matrix"
			res
		}
		
		res <- switch(type,
			class = factor(as.character(object$predicted),
				levels = levels(object)),
			membership = .membership(toProps(object$votes, object$ntree),
				levels = levels(object)),
			both = list(class = factor(as.character(object$predicted),
				levels = levels(object)),
				membership = .membership(toProps(object$votes, object$ntree),
				levels = levels(object))),
			vote = .membership(toVotes(object$votes, object$ntree),
						levels = levels(object)),
			stop("unknown type, must be 'class', 'membership', 'both' or 'vote'"))
		
		attr(res, "method") <- list(name = "out-of-bag")
		res
		
	} else predict.mlearning(object = object, newdata = newdata,
		type = type, norm.votes = FALSE, ...)
}

mlNnet <- function (...)
	UseMethod("mlNnet")

mlNnet.formula <- function (formula, data, size = NULL, rang = NULL, decay = 0,
maxit = 1000, ..., subset, na.action)
	mlearning(formula, data = data, method = "mlNnet", model.args =
		list(formula  = formula, data = substitute(data),
		subset = substitute(subset)), call = match.call(), size = size,
		rang = rang, decay = decay, maxit = maxit, ...,
		subset = subset, na.action = substitute(na.action))

mlNnet.default <- function (train, response, size = NULL, rang = NULL, decay = 0,
maxit = 1000, ...)
{
	if (!length(response))
		stop("unsupervised classification not usable for mlNnet")

	nnetArgs <- dots <- list(...)
	.args. <- nnetArgs$.args.
	dots$.args. <- NULL
	dots$size <- size
	dots$rang <- rang
	dots$decay <- decay
	dots$maxit <- maxit
	nnetArgs$.args. <- NULL
	if (!length(.args.)) .args. <- list(levels = levels(response),
		n = c(intial = NROW(train), final = NROW(train)),
		type = if (is.factor(response)) "classification" else "regression",
		na.action = "na.pass", mlearning.call = match.call(), method = "mlNnet")
	
	## Construct arguments list for nnet() call
	nnetArgs$x <- sapply(train, as.numeric)
	
	## Weights
	if (!length(nnetArgs$weights)) nnetArgs$weights <- .args.$weights
	
	## size
	if (!length(size))
		size <- length(levels(response)) - 1 # Is this a reasonable default?
	nnetArgs$size <- size
	
	## rang
	if (!length(rang)) {
		## default is 0.7 in original nnet code,
		## but the doc proposes something else
		rang <- round(1 / max(abs(nnetArgs$x)), 2)
		if (rang < 0.01) rang <- 0.01
		if (rang > 0.7) rang <- 0.7
	}
	nnetArgs$rang <- rang
	
	## decay and maxit
	nnetArgs$decay <- decay
	nnetArgs$maxit <- maxit
			
	## TODO: should I need to implement this???
	#x <- model.matrix(Terms, m, contrasts)
	#cons <- attr(x, "contrast")
	#xint <- match("(Intercept)", colnames(x), nomatch = 0L)
	#if (xint > 0L) 
	#    x <- x[, -xint, drop = FALSE]
			
	## Classification or regression?
	if (is.factor(response)) {
		if (length(levels(response)) == 2L) {
			nnetArgs$y <- as.vector(unclass(response)) - 1
			nnetArgs$entropy <- TRUE
			res <- do.call(nnet.default, nnetArgs)
			res$lev <- .args.$levels
		} else {
			nnetArgs$y <- class.ind(response)
			nnetArgs$softmax <- TRUE
			res <- do.call(nnet.default, nnetArgs)
			res$lev <- .args.$levels
		}
	} else { # Regression
		nnetArgs$y <- response
		res <- do.call(nnet.default, nnetArgs)	
	}
	
	## Return a mlearning object
	structure(res, formula = .args.$formula, train = train,
		response = response, levels = .args.$levels, n = .args.$n, args = dots,
		optim = .args.$optim, numeric.only = TRUE, type = .args.$type,
		pred.type = c(class = "class", membership = "raw"),
		summary = "summary", na.action = .args.$na.action,
		mlearning.call = .args.$mlearning.call, method = .args.$method,
		algorithm = "single-hidden-layer neural network",
		class = c("mlNnet", "mlearning", "nnet"))
}

mlLvq <- function (...)
	UseMethod("mlLvq")

mlLvq.formula <- function (formula, data, k.nn = 5, size, prior,
algorithm = "olvq1", ..., subset, na.action)
{
	if (missing(size)) {
		if (missing(prior)) {
			mlearning(formula, data = data, method = "mlLvq", model.args =
				list(formula  = formula, data = substitute(data),
				subset = substitute(subset)), call = match.call(), k.nn = k.nn,
				algorithm = algorithm, ...,
				subset = subset, na.action = substitute(na.action))
		} else {
			mlearning(formula, data = data, method = "lvq", model.args =
				list(formula  = formula, data = substitute(data),
				subset = substitute(subset)), call = match.call(), k.nn = k.nn,
				prior = prior, algorithm = algorithm, ...,
				subset = subset, na.action = substitute(na.action))
		}
	} else {
		if (missing(prior)) {
			mlearning(formula, data = data, method = "lvq", model.args =
				list(formula  = formula, data = substitute(data),
				subset = substitute(subset)), call = match.call(), k.nn = k.nn,
				size = size, algorithm = algorithm, ...,
				subset = subset, na.action = substitute(na.action))
		} else {
			mlearning(formula, data = data, method = "lvq", model.args =
				list(formula  = formula, data = substitute(data),
				subset = substitute(subset)), call = match.call(), k.nn = k.nn,
				size = size, prior = prior, algorithm = algorithm, ...,
				subset = subset, na.action = substitute(na.action))
		}
	}
}

mlLvq.default <- function (train, response, k.nn = 5, size, prior,
algorithm = "olvq1", ...)
{
	if (!is.factor(response))
		stop("only factor response (classification) accepted for mlLvq")

	dots <- list(...)
	.args. <- dots$.args.
	dots$.args. <- NULL
	dots$k.nn <- k.nn
	dots$algorithm <- algorithm
	if (!length(.args.)) .args. <- list(levels = levels(response),
		n = c(intial = NROW(train), final = NROW(train)),
		type = "classification", na.action = "na.pass",
		mlearning.call = match.call(), method = "mlQda")
	
	## matrix of numeric values
	if (any(sapply(train, is.factor))) {
		warning("force conversion from factor to numeric; may be not optimal or suitable")
		train <- sapply(train, as.numeric)
	}
			
	## Default values for size and prior, if not provided
	n <- nrow(train)
	if (missing(prior) || !length(prior)) {
		prior <- tapply(rep(1, length(response)), response, sum) / n
	} else dots$prior <- prior
	if (missing(size) || !length(size)) {
		np <- length(prior)  
		size <- min(round(0.4 * np * (np - 1 + ncol(train) / 2), 0), n)
	} else dots$size <- size
			
	## Initialize codebook
	init <- lvqinit(train, response, k = k.nn, size = size, prior = prior)
	
	## Calculate final codebook
	if (algorithm == "olvq1") times <- 40 else times <- 100
	niter <- dots$niter
	if (!length(niter)) niter <- times * nrow(init$x) # Default value
	alpha <- dots$alpha
	if (!length(alpha)) alpha <- if (algorithm == "olvq1") 0.3 else 0.03
	win <- dots$win
	if (!length(win)) win <- 0.3
	epsilon <- dots$epsilon
	if (!length(epsilon)) epsilon <- 0.1
	codebk <- switch(algorithm,
		olvq1 = olvq1(train, response, init, niter = niter,
			alpha = alpha),
		lvq1 = lvq1(train, response, init, niter = niter,
			alpha = alpha),
		lvq2 = lvq2(train, response, init, niter = niter,
			alpha = alpha, win = win),
		lvq3 = lvq3(train, response, init, niter = niter,
			alpha = alpha, win = win, epsilon = epsilon),
		stop("algorithm must be 'lvq1', 'lvq2', 'lvq3' or 'olvq1'"))
				
	## Return a mlearning object
	structure(codebk, formula = .args.$formula, train = train,
		response = response, levels = .args.$levels, n = .args.$n, args = dots,
		optim = .args.$optim, numeric.only = TRUE, type = .args.$type,
		pred.type = c(class = "class"), summary = "summary.lvq",
		na.action = .args.$na.action,
		mlearning.call = .args.$mlearning.call, method = .args.$method,
		algorithm = "learning vector quantization",
		class = c("mlLvq", "mlearning", class(codebk)))
}

summary.lvq <- function (object, ...)
	structure(cbind(Class = object$cl, as.data.frame(object$x)),
		class = c("summary.lvq", "data.frame"))

print.summary.lvq <- function (x, ...)
{
	cat("Codebook:\n")
	print(as.data.frame(x))
	invisible(x)
}

predict.mlLvq <- function (object, newdata, type = "class",
method = c("direct", "cv"), na.action = na.exclude, ...)
{
	if (!inherits(object, "mlLvq"))
		stop("'object' must be a 'mlLvq' object")
	if (type != "class") stop("Only 'class' currently supported for type")
    
	## If method == "cv", delegate to cvpredict()
	if (as.character(method)[1] == "cv") {
		if (!missing(newdata))
			stop("cannot handle new data with method = 'cv'")
		return(cvpredict(object = object, type = type, ...))
	}
	
	## Recalculate newdata according to formula...
	if (missing(newdata)) { # Use train
		newdata <- attr(object, "train")
	} else if (attr(object, "optim")) { # Use optimized approach
		## Just keep vars similar as in train
		vars <- colnames(attr(object, "train"))
		if (!all(vars %in% names(newdata)))
			stop("one or more missing variables in newdata")
		newdata <- newdata[, vars]
	} else { # Use model.frame
		newdata <- model.frame(formula = attr(object, "formula"),
			data = newdata, na.action = na.pass)[, names(attr(object, "train"))]
	}
	newdata <- sapply(as.data.frame(newdata), as.numeric)
	
	## Determine how many data and perform na.action
	n <- NROW(newdata)
	newdata <- match.fun(na.action)(newdata)
	ndrop <- attr(newdata, "na.action")
	
	.expandFactor(lvqtest(object, newdata), n, ndrop)
}

## svm from e1071 package
mlSvm <- function (...)
	UseMethod("mlSvm")

mlSvm.formula <- function(formula, data, scale = TRUE, type = NULL,
kernel = "radial", classwt = NULL, ..., subset, na.action)
	mlearning(formula, data = data, method = "mlSvm", model.args =
		list(formula  = formula, data = substitute(data),
		subset = substitute(subset)), call = match.call(),
		..., subset = subset, na.action = substitute(na.action))

mlSvm.default <- function (train, response, scale = TRUE, type = NULL,
kernel = "radial", classwt = NULL, ...)
{
	dots <- list(...)
	.args. <- dots$.args.
	dots$.args. <- NULL
	if (!length(.args.)) {
		if (is.factor(response)) {
			Type <- "classification"
		} else Type <- "regression"
		.args. <- list(levels = levels(response),
		n = c(intial = NROW(train), final = NROW(train)),
		type = Type, na.action = "na.pass",
		mlearning.call = match.call(), method = "mlSvm")
	}
	dots$scale <- scale
	dots$type <- type
	dots$kernel <- kernel
	dots$class.weigths <- classwt
	#dots$probability <- TRUE
	
	## Return a mlearning object
	structure(e1071:::svm.default(x = sapply(train, as.numeric), y = response,
		scale = scale, type = type, kernel = kernel, class.weights = classwt,
		probability = TRUE, ...), formula = .args.$formula, train = train,
		response = response, levels = .args.$levels, n = .args.$n, args = dots,
		optim = .args.$optim, numeric.only = TRUE, type = .args.$type,
		pred.type = c(class = "class", membership = "raw"),
		summary = "summary", na.action = .args.$na.action,
		mlearning.call = .args.$mlearning.call, method = .args.$method,
		algorithm = "support vector machine",
		class = c("mlSvm", "mlearning", "svm"))
}

predict.mlSvm <- function(object, newdata,
type = c("class", "membership", "both"), method = c("direct", "cv"),
na.action = na.exclude, ...)
{
	if (!inherits(object, "mlSvm"))
		stop("'object' must be a 'mlSvm' object")
	
	## If method == "cv", delegate to cvpredict()
	method <- as.character(method)[1]
	if (method == "cv") {
		if (!missing(newdata))
			stop("cannot handle new data with method = 'cv'")
		return(cvpredict(object = object, type = type, ...))	
	}
	
	## Recalculate newdata according to formula...
	if (missing(newdata)) { # Use train
		newdata <- attr(object, "train")
	} else if (attr(object, "optim")) { # Use optimized approach
		## Just keep vars similar as in train
		vars <- names(attr(object, "train"))
		if (!all(vars %in% names(newdata)))
			stop("One or more missing variables in newdata")
		newdata <- newdata[, vars]
	} else { # Use model.frame
		newdata <- model.frame(formula = attr(object, "formula"),
			data = newdata, na.action = na.pass)[, names(attr(object, "train"))]
	}
	## Only numerical predictors
	newdata <- sapply(as.data.frame(newdata), as.numeric)
	
	## Determine how many data and perform na.action
	n <- NROW(newdata)
	newdata <- match.fun(na.action)(newdata)
	ndrop <- attr(newdata, "na.action")
	attr(newdata, "na.action") <- NULL
		
	## Delegate to the e1071 predict.svm method
	if (as.character(type)[1] == "class") proba <- FALSE else proba <- TRUE
	class(object) <- class(object)[-(1:2)]
	if (attr(object, "type") == "regression")
		return(predict(object, newdata = newdata, ...))
	
	## This is for classification
	res <- predict(object, newdata = newdata,
		probability = proba, ...)
	proba <- attr(res, "probabilities")	
	
	## Rework results according to what we want
	switch(as.character(type)[1],
		class = .expandFactor(factor(as.character(res), levels = levels(object)),
			n, ndrop),
		membership = .expandMatrix(.membership(proba, levels = levels(object)),
			n, ndrop),
		both = list(class = .expandFactor(factor(as.character(res),
			levels = levels(object)), n, ndrop),
			membership = .expandMatrix(.membership(proba, levels = levels(object)),
			n, ndrop)),
		stop("unrecognized 'type' (must be 'class', 'membership' or 'both')"))
}

## NaiveBayes from e1071 package
mlNaiveBayes <- function (...)
	UseMethod("mlNaiveBayes")

mlNaiveBayes.formula <- function(formula, data, laplace = 0, ...,
subset, na.action)
	mlearning(formula, data = data, method = "mlNaiveBayes", model.args =
		list(formula  = formula, data = substitute(data),
		subset = substitute(subset)), call = match.call(), laplace = laplace,
		..., subset = subset, na.action = substitute(na.action))

mlNaiveBayes.default <- function (train, response, laplace = 0, ...)
{
	if (!is.factor(response))
		stop("only factor response (classification) accepted for mlNaiveBayes")

	dots <- list(...)
	.args. <- dots$.args.
	dots$.args. <- NULL
	dots$laplace <- laplace
	if (!length(.args.)) .args. <- list(levels = levels(response),
		n = c(intial = NROW(train), final = NROW(train)),
		type = "classification", na.action = "na.pass",
		mlearning.call = match.call(), method = "mlNaiveBayes")
	
	## Return a mlearning object
	structure(e1071:::naiveBayes.default(x = train, y = response,
		laplace = laplace, ...), formula = .args.$formula, train = train,
		response = response, levels = .args.$levels, n = .args.$n, args = dots,
		optim = .args.$optim, numeric.only = FALSE, type = .args.$type,
		pred.type = c(class = "class", membership = "raw"),
		summary = NULL, na.action = .args.$na.action,
		mlearning.call = .args.$mlearning.call, method = .args.$method,
		algorithm = "naive Bayes classifier",
		class = c("mlNaiveBayes", "mlearning", "naiveBayes"))
}
	
## NaiveBayes from RWeka package
## TODO: keep this for mlearningWeka package!
#mlNaiveBayesWeka <- function (...)
#	UseMethod("mlNaiveBayesWeka")
#
#mlNaiveBayesWeka.formula <- function(formula, data, ..., subset, na.action)
#	mlearning(formula, data = data, method = "mlNaiveBayesWeka", model.args =
#		list(formula  = formula, data = substitute(data),
#		subset = substitute(subset)), call = match.call(),
#		..., subset = subset, na.action = substitute(na.action))
#
#mlNaiveBayesWeka.default <- function (train, response, ...)
#{
#	if (!is.factor(response))
#		stop("only factor response (classification) accepted for mlNaiveBayesWeka")
#
#	.args. <- dots <- list(...)$.args.
#	if (!length(.args.)) .args. <- list(levels = levels(response),
#		n = c(intial = NROW(train), final = NROW(train)),
#		type = "classification", na.action = "na.pass",
#		mlearning.call = match.call(), method = "mlNaiveBayesWeka")
#	
#	wekaArgs <- list(control = .args.$control)
#	
#	## If response is not NULL, add it to train
#	if (length(response)) {
#		formula <- .args.$formula
#		if (!length(formula)) response.label <- "Class" else
#			response.label <- all.vars(formula)[1]
#		data <- data.frame(response, train)
#		names(data) <- c(response.label, colnames(train))
#		wekaArgs$data <- data
#		wekaArgs$formula <- as.formula(paste(response.label, "~ ."))
#	} else { # Unsupervised classification
#		wekaArgs$data <- train
#		wekaArgs$formula <- ~ . 
#	}
#	
#	WekaClassifier <- make_Weka_classifier("weka/classifiers/bayes/NaiveBayes")
#	
#	## Return a mlearning object
#	structure(do.call(WekaClassifier, wekaArgs), formula = .args.$formula,
#		train = train, response = response, levels = .args.$levels, n = .args.$n,
#		args = dots, optim = .args.$optim, numeric.only = FALSE,
#		type = .args.$type, pred.type = c(class = "class", membership = "probability"),
#		summary = "summary", na.action = .args.$na.action,
#		mlearning.call = .args.$mlearning.call, method = .args.$method,
#		algorithm = "Weka naive Bayes classifier",
#		class = c("mlNaiveBayesWeka", "mlearning", "Weka_classifier"))
#}
