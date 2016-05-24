## Linear Predictor Score fitting
## Author : Sylvain Mareschal <maressyl@gmail.com>
LPS <- function(data, coeff, response, k, threshold, formula, method="fdr", ...) {
	# Checks
	if(!is.matrix(coeff) || !is.numeric(coeff))                         stop("'coeff' must be a numeric matrix")
	if(is.null(rownames(coeff)))                                        stop("'coeff' must have row names")
	if(is.null(colnames(coeff)))                                        stop("'coeff' must have col names")
	if(!is.data.frame(data) && (!is.matrix(data) || !is.numeric(data))) stop("'data' must be a data.frame or a numeric matrix")
	
	# formula, k, threshold
	if(!missing(k)) {
		# Build K formula
		type <- "k"
		k <- as.integer(k)
		formula <- LPS.k.form(responseName=".", coeff=coeff[,1], k=k)
	} else if(!missing(threshold)) {
		# Build P formula
		type <- "p"
		threshold <- as.double(threshold)
		formula <- LPS.p.form(responseName=".", p.values=coeff[,"p.value"], threshold=threshold, method=method)
	} else if(!missing(formula)) {
		# Extract 'k' or 'threshold' from 'formula'
		if(length(formula) == 3 && (is.name(formula[[3]]) || is.atomic(formula[[3]]))) {
			# Single scalar in variables
			if(formula[[3]] == ".") {
				# '.' for all genes (pseudo K formula)
				type <- "k"
				k <- nrow(coeff)
				formula <- update.formula(formula, LPS.k.form(responseName=".", coeff=coeff[,1], k=k))
			} else {
				# Gene selection
				value <- eval(formula[[3]], parent.frame())
				if(is.integer(value)) {
					# Update K formula
					type <- "k"
					k <- value
					formula <- update.formula(formula, LPS.k.form(responseName=".", coeff=coeff[,1], k=k))
				} else if(is.double(value)) {
					# Update P formula
					type <- "p"
					threshold <- value
					formula <- update.formula(formula, LPS.p.form(responseName=".", p.values=coeff[,"p.value"], threshold=threshold, method=method))
				}
			}
		} else {
			# Custom formula
			type <- "custom"
		}
	} else stop("Please provide a gene list or selection criterion via 'k', 'threshold' or 'formula'")
	
	# Extract data
	if(missing(response)) {
		# Model frame (slower)
		mf <- model.frame(formula=formula, data=data, ...)
		response <- mf[,1L]
	} else {
		# Subsetting rows (matrix)
		argList <- list(...)
		if("subset" %in% names(argList) && !is.null(argList[["subset"]])) {
			# Subset data
			if(is.null(rownames(data)) && is.character(argList[["subset"]])) stop("'subset' is a character vector but 'data' rows are unnamed")
			data <- data[ argList[["subset"]] ,]
			
			# Subset response
			if(is.null(names(response)) && is.character(argList[["subset"]])) stop("'subset' is a character vector but 'response' is unnamed")
			response <- response[ argList[["subset"]] ]
		}
		
		# No model.frame
		mf <- NULL
	}
	
	# Response factor
	if(is.character(response) || is.numeric(response) || is.logical(response)) response <- factor(response)
	if(length(levels(response)) != 2) stop("The response element of the 'formula' must have 2 levels")
	
	# Coefficient selection (keeping name, formula already expanded)
	columns <- var.formula(formula, backtick=FALSE)
	coeffSel <- coeff[,1][ columns , drop=FALSE ]
	tSel <- coeff[,"t"][ columns , drop=FALSE ]
	
	# Expression matrix
	expr <- as.matrix(data[ , columns , drop=FALSE ])
	
	# NAs and subsetting (model.frame)
	if(!is.null(mf)) {
		expr <- expr[ sub("\\.[0-9]+", "", rownames(mf)) , , drop=FALSE ]
	}
	
	# Expression matrixes
	M1 <- expr[ response == levels(response)[1] , , drop=FALSE ]
	M2 <- expr[ response == levels(response)[2] , , drop=FALSE ]
	
	# LPS computation
	LPS <- list(
		apply(t(M1) * coeffSel, 2, sum),
		apply(t(M2) * coeffSel, 2, sum)
	)
	
	# Model object
	object <- list(
		coeff = coeffSel,
		t = tSel,
		classes = levels(response),
		scores = lapply(LPS, "attributes<-", NULL),
		means = sapply(LPS, mean),
		sds = sapply(LPS, sd)
	)
	object$ovl <- OVL(means=object$means, sds=object$sds)
	if(type == "k") object$k <- k
	if(type == "p") {
		object$p.threshold <- threshold
		object$p.method <- method
	}
	class(object) <- "LPS"
	
	return(object)
}
