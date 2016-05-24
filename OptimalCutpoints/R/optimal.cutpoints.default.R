optimal.cutpoints.default <-
function(X, status, tag.healthy, methods, data, direction = c("<", ">"), categorical.cov = NULL, pop.prev = NULL, control = control.cutpoints(), ci.fit = FALSE, conf.level = 0.95, trace = FALSE, ...) {
	if(missing(methods) || is.null(methods)) {
		stop("'methods' argument required.", call.=FALSE)
	}
	if(any(!(methods %in% c("CB","MCT","MinValueSp","MinValueSe","ValueSp","ValueSe","MinValueSpSe", "MaxSp", "MaxSe", "MaxSpSe",
		"MaxProdSpSe","ROC01","SpEqualSe","Youden","MaxEfficiency","Minimax","MaxDOR","MaxKappa",
		"MinValueNPV","MinValuePPV","ValueNPV","ValuePPV","MinValueNPVPPV","PROC01","NPVEqualPPV","MaxNPVPPV", "MaxSumNPVPPV", "MaxProdNPVPPV","ValueDLR.Negative","ValueDLR.Positive","MinPvalue","ObservedPrev","MeanPrev","PrevalenceMatching")))) {
		stop ("You have entered an invalid method.", call. = FALSE)
	}
	if (missing(data)|| is.null(data)) {
		stop("'data' argument required.", call. = FALSE)
	}
	if (missing(X)|| is.null(X)) {
		stop("'X' argument required.", call. = FALSE)
	}
	if (missing(status)|| is.null(status)) {
		stop("'status' argument required.", call. = FALSE)
	}
	if (missing(tag.healthy)|| is.null(tag.healthy)) {
		stop("'tag.healthy' argument required.", call. = FALSE)
	}
	if (is.logical(ci.fit) == FALSE) {
		stop("'ci.fit' must be a logical-type argument.", call. = FALSE)
	}	
	if (conf.level < 0 | conf.level > 1 | length(conf.level) != 1) {
		stop("'conf.level' must be a single number between 0 and 1.", call. = FALSE)
	}
	if (is.logical(trace) == FALSE) {
		stop("'trace' must be a logical-type argument.", call. = FALSE)
	}
	if (is.null(pop.prev) & ci.fit == TRUE & !control$ci.PV %in% c("Exact","Quadratic","Wald","AgrestiCoull","RubinSchenker")) {	 
		warning(paste("Predictive Vaues CI: ``",control$ci.PV,"'' method is not valid when prevalence is estimated from the sample.\n", sep = ""), call. = FALSE)
	}
	if (!is.null(pop.prev) & ci.fit == TRUE & !control$ci.PV %in% c("Transformed","NotTransformed","GartNam")) {
		warning(paste("Predictive Values CI: \"",control$ci.PV,"\" method is not valid when prevalence is not estimated from the sample.\n", sep = ""), call. = FALSE)
	}
	direction <- match.arg(direction)
	if(!all(c(X,status,categorical.cov) %in% names(data))) {
		stop("Not all needed variables are supplied in 'data'.", call. = FALSE)
	}
	# NA's deleted
	data <- na.omit(data[,c(X,status,categorical.cov)])
	# A data frame with the results is created:
	res <- vector("list", length(methods))
	names(res) <- methods
	# Categorical covariate levels:	 	
	if(!is.null(categorical.cov)) {
		if(!is.factor(data[, categorical.cov])) data[, categorical.cov] <- factor(data[, categorical.cov])
		data[, categorical.cov] <- droplevels(data[, categorical.cov])
		levels.cat <- levels(data[, categorical.cov])			
		for (i in 1: length(methods)) {
			res[[i]] <- vector("list", length(levels.cat))
			names(res[[i]]) <- levels.cat
		}		
	} else {
		levels.cat = 1
		res[[1]] <- vector("list", 1)
		names(res[[1]]) <- "Global"
	}
	pop.prev.new <- vector(length=length(levels(data[, categorical.cov])))
	if(is.null(pop.prev)) pop.prev <- NA
	if (!is.null(categorical.cov) & length(pop.prev) != 1 & length(pop.prev) != length(levels(data[, categorical.cov]))) {
		stop("You have entered different values for prevalence which \n do not coincide with categorical covariate levels.", call. = FALSE)
	} else if (!is.null(categorical.cov) & length(pop.prev) == 1) {
		pop.prev.new <- rep(pop.prev, length(levels(data[, categorical.cov])))				 
	} else if (is.null(categorical.cov) & length(pop.prev) > 1) {
 		warning("You have entered several values for prevalence. \n The first value has been selected.", call. = FALSE, immediate. = TRUE)	 		
 		pop.prev.new <- pop.prev[1]			
	} else {
		pop.prev.new <- pop.prev
	}
	# Each method is called up:
	for(i in 1:length(levels.cat)) {
		if(trace) {
			if(length(levels.cat) > 1) {
				text <- paste("Level: ", levels.cat[i], sep = "")
				cat(text)
				cat("\nAnalysing ...\n\n")
			}
		}
		data.m <- if(length(levels.cat) != 1) data[data[,categorical.cov] == levels.cat[i], ] else data
		if (is.na(pop.prev.new[i])) {
			pop.prev.new[i] <- calculate.sample.prev(data.m, status, tag.healthy)
		}
		validate.prevalence(pop.prev.new[i])
		measures.acc <- calculate.accuracy.measures(data = data.m, marker = X, status = status, tag.healthy = tag.healthy, direction = direction, pop.prev = pop.prev.new[i], control = control, conf.level = conf.level, ci.fit = ci.fit)																
		for (j in 1: length(methods)) {  		 
			if(trace) {
				text <- paste("Method: ", methods[j],sep = "")
				cat(text)
				cat("\nAnalysing ...\n\n")
			}			
			res[[j]][[i]] <- eval(parse(text = paste("function.", methods[j], sep = "")))(data = data.m, marker = X, status = status, tag.healthy = tag.healthy, direction = direction, pop.prev = pop.prev.new[i], control = control, conf.level = conf.level, ci.fit = ci.fit, measures.acc = measures.acc)
		}
	}
	res$methods <-  methods
	if(length(levels.cat) != 1) res$levels.cat  <-  levels.cat
	res$call <- match.call()
	res$data <- data
	class(res) <- "optimal.cutpoints"
	invisible(res)
	res
}
