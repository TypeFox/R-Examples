
.apply.MPT.restrictions <- function(tree, restrictions) {
	safeDeparse <- function(expr){
		ret <- paste(deparse(expr), collapse="")
		#rm whitespace
		ret <- gsub("[[:space:]][[:space:]]+", " ", ret)
		gsub("^expression\\(", "", gsub("[[:space:]][[:space:]]", " ", gsub("\\)$", "", ret)))
	}
	for (c1 in 1:length(tree)){
		for (c2 in 1:length(tree[[c1]])) {
			for (c3 in 1:length(restrictions)) {
				tree[[c1]][c2] <- parse(text = gsub(paste("\\<",restrictions[[c3]][1], "\\>", sep =""),restrictions[[c3]][2],safeDeparse(tree[[c1]][c2])))[1]
			}
		}
	}
	return(tree)
}


.check.restrictions <- function(restrictions.filename, model) {
	restrictions.tmp <- .read.MPT.restrictions(restrictions.filename)
	restricted.parameter <- vapply(restrictions.tmp, "[", "", i = 1)
	orig.tree <- model
	orig.params <- .find.MPT.params(model)
	#browser()
	if (any(!restricted.parameter %in% orig.params)) {
		warning(paste("Restricted parameter(s)",  paste(restricted.parameter[!restricted.parameter %in% orig.params], collapse = ", "), "is/are not in the original model."))
		restrictions <- restrictions.tmp[restricted.parameter %in% orig.params]
		restrictions.tmp <- restrictions
	} else restrictions <- restrictions.tmp
	restricted.parameter <- vapply(restrictions.tmp, "[", "", i = 2)
	restricted.parameter[!is.na(suppressWarnings(as.numeric(restricted.parameter)))] <- orig.params[1]
	restricted.parameter <- vapply(strsplit(restricted.parameter, "\\*"), "[", "", i = 1)
	restricted.parameter <- gsub(" ", "", restricted.parameter)
	#browser()
	orig.tree <- model
	orig.params <- .find.MPT.params(model)
	#browser()
	if (any(!restricted.parameter %in% orig.params)) {
		warning(paste("Restricted parameter(s)",  paste(restricted.parameter[!restricted.parameter %in% orig.params], collapse = ", "), "is/are not in the original model."))
		restrictions <- restrictions.tmp[restricted.parameter %in% orig.params]
	} else restrictions <- restrictions.tmp
	restrictions
}
	
