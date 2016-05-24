summary.gsaResult <- function(object,
		mode = c("summary", "table"),
		orderBy = c("adjustedPValues", "rawPValues", "geneSetName"),
		significantOnly = FALSE,
		signLevel = object$signLevel,
		...){
	
	if(class(object) != "gsaResult"){
		stop("'object' mut be of class gsaResult.")
	}
	mode <- match.arg(mode)

	if(mode == "summary"){
		cat("\n    Analysis: ", object$analysis$name, "\n", sep ="")

		cat("\n    ", length(object$rawPValues),
			" gene set(s) tested:", "\n", sep ="")

		cat("      - ", sum(object$rawPValues < object$signLevel),
			" gene set(s) with raw p-value < ",
			object$signLevel, "\n", sep ="")
		cat("      - min p-value: ",
			names(object$res.all)[which.min(object$rawPValues)], " (", min(object$rawPValues), ")", "\n\n", sep ="")

		cat("    ", "Correction for multiple testing: ",
			object$adjustmentMethod, "\n", sep ="")

		cat("      - ", sum(object$adjustedPValues < object$signLevel),
			" gene set(s) with adjusted p-value < ",
			object$signLevel, "\n\n", sep ="")
	
	}else{
		print(createSummaryTable(object,
			orderBy = orderBy,
			significantOnly = significantOnly,
			signLevel = signLevel))
	}

	return(invisible(NULL))
}


createSummaryTable <- function(object,
		orderBy = c("adjustedPValues", "rawPValues", "geneSetName"),
		significantOnly = FALSE,
		signLevel = object$signLevel){

	if(class(object) != "gsaResult"){
		stop("'object' mut be of class gsaResult.")
	}

	orderBy <- match.arg(orderBy)

	#determine intersect size
	intersects <- sapply(object$res.all, function(i){
			return(length(i$geneSetValues$intersectGeneSetCoreSet))
		})

	#if all intersect sizes are 0 no fisher test has been performed and
	#the sizes must not be printed
	if(sum(intersects) == 0){
		res <- data.frame(geneSetName = names(object$res.all),
			adjustedPValues = object$adjustedPValues,
			rawPValues = object$rawPValues,
			geneSetSize = sapply(object$res.all, function(i){return(length(i$geneSet))}),
			stringsAsFactors=FALSE)
	}else{
		res <- data.frame(geneSetName = names(object$res.all),
			adjustedPValues = object$adjustedPValues,
			rawPValues = object$rawPValues,
			geneSetSize = sapply(object$res.all, function(i){return(length(i$geneSet))}),
			intersectSize = intersects,
			stringsAsFactors=FALSE)
	}

	#order table according to parameter
	indOrdered <- order(res[[orderBy]])
	res <- res[indOrdered,]

	if(significantOnly){
		ind <- which(res[["adjustedPValues"]] < signLevel)
		res <- res[ind,]
	}

	rownames(res) <- NULL

	return(res)
}