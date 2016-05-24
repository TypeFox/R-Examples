## Linear Predictor Score fast coefficient computation
## Author : Sylvain Mareschal <maressyl@gmail.com>
LPS.coeff <- function(data, response, formula=~1, type=c("t", "limma"), p.value=TRUE, log=FALSE, weighted=FALSE, ...) {
	# Checks
	type <- match.arg(type)
	if(!is(formula, "formula"))                                         stop("'formula' must be a formula object")
	if(!is.data.frame(data) && (!is.matrix(data) || !is.numeric(data))) stop("'data' must be a data.frame or a numeric matrix")
	
	# Extract data
	if(missing(response)) {
		# Model frame (slower)
		mf <- model.frame(formula=formula, data=data, ...)
		response <- mf[,1L]
	} else {
		# Subsetting rows (matrix)
		argList <- list()
		if("subset" %in% names(argList) && !is.null(argList[["subset"]])) {
			data <- data[ argList[["subset"]] ,]
			response <- response[ argList[["subset"]] ]
		}
		
		# No model.frame
		mf <- NULL
	}
	
	# Response factor
	if(is.character(response) || is.numeric(response) || is.logical(response)) response <- factor(response)
	if(length(levels(response)) != 2) stop("The response element of the 'formula' must have 2 levels")
	
	# Column subsetting
	if(grepl("~ *1 *$", deparse(formula, nlines=1))) {
		if(is.matrix(data)) {
			# All columns for matrix 'data'
			columns <- TRUE
		} else {
			# All columns except from non numeric ones for data.frame 'data'
			columns <- -which(!sapply(data, is.numeric))
			if(length(columns) == 0) columns <- TRUE
		}
	} else {
		# Only columns in formula
		if(grepl("~ *\\. *$", deparse(formula, nlines=1)) && is.matrix(data)) stop("'.' in 'formula' is not compatible with a matrix 'data'")
		columns <- attr(terms(formula, data=data), "term.labels")
		columns <- sub("`(.+)`", "\\1", columns)
		
	}
	
	# Expression matrix
	expr <- as.matrix(data[ , columns , drop=FALSE ])
	
	# NAs and subsetting (model.frame)
	if(!is.null(mf)) {
		expr <- expr[ sub("\\.[0-9]+", "", rownames(mf)) , , drop=FALSE ]
	}
	
	# LIMMA T
	if(type == "limma") {
		# LIMMA fit
		design <- cbind(base=1, diff=as.integer(response == levels(response)[2]))
		fit <- limma::lmFit(t(expr), design)
		fit <- limma::eBayes(fit)
		out <- fit$t[,"diff"]
		
		# P-value
		if(isTRUE(p.value)) { out <- cbind(t=out, p.value=fit$p.value[,"diff"])
		} else              { out <- matrix(out, dimnames=list(names(out), "t"))
		}
	} else {
		# Expression matrixes
		M1 <- expr[ response == levels(response)[1] , , drop=FALSE ]
		M2 <- expr[ response == levels(response)[2] , , drop=FALSE ]
		
		# Probe means
		s1 <- apply(M1, 2, sum, na.rm=TRUE)
		s2 <- apply(M2, 2, sum, na.rm=TRUE)
		n1 <- apply(!is.na(M1), 2, sum)
		n2 <- apply(!is.na(M2), 2, sum)
		m1 <- s1 / n1
		m2 <- s2 / n2
		
		# T statistic (alternative="two.sided", paired=FALSE, var.equal=FALSE)
		v1 <- apply(M1, 2, var, na.rm=TRUE)
		v2 <- apply(M2, 2, var, na.rm=TRUE)
		out <- (m1 - m2) / sqrt(v1/n1 + v2/n2)
		
		# Student's p-value
		if(isTRUE(p.value)) {
			dfr <- ((v1/n1 + v2/n2)^2) / (((v1/n1)^2)/(n1-1) + ((v2/n2)^2)/(n2-1))
			pval <- pt(q=abs(out), df=dfr, lower.tail=FALSE)*2
			out <- cbind(t=out, p.value=pval)
		} else {
			out <- matrix(out, dimnames=list(names(out), "t"))
		}
	}
	
	# Log transformation
	if(isTRUE(log)) {
		log.t <- double(nrow(out))
		log.t[ out[,"t"] >  1 ] <-  log( out[ out[,"t"] >  1 , "t" ])
		log.t[ out[,"t"] < -1 ] <- -log(-out[ out[,"t"] < -1 , "t" ])
		out <- cbind(log.t, out)
	}
	
	# Mean weighting
	if(isTRUE(weighted)) {
		out <- cbind(
			out[,1] / apply(expr, 2, mean, na.rm=TRUE),
			out
		)
		colnames(out)[1] <- sprintf("weighted.%s", colnames(out)[2])
	}
	
	return(out)
}
