# OAPPLY ########### tapply of ovariables applies a function to each cell of a ragged array, that is to each (non-empty) group of
############ values given by a unique combination of the levels of certain factors.
### parameters (other parameters are as in generic tapply):
### X an ovariable
### cols overrides INDEX by choosing INDEX as all marginals NOT given in cols (character vector) argument

oapply = function(X, INDEX = NULL, FUN = NULL, cols = NULL, #use_plyr = FALSE, 
		drop_na = TRUE, use_aggregate = TRUE, ..., simplify = TRUE) {
	if(!use_aggregate) out <- X@output
	marginals <- colnames(X@output)[X@marginal]
	if (is.data.frame(INDEX)) INDEX <- colnames(INDEX)
	if (is.null(INDEX) & is.null(cols)) stop("No INDEX nor cols defined!\n")
	if (!is.null(cols)) INDEX <- marginals[!marginals %in% cols]
	if (length(INDEX) == 0) {
		warning("Zero length INDEX while oapplying. All columns except relevant Result removed.")
		res <- FUN(result(X))
		X@output <- data.frame(res)
		colnames(X@output) <- paste(X@name, "Result", sep = "")
		X@marginal <- FALSE
		return(X)
	}
	if (use_aggregate) {
		out <- aggregate(result(X), X@output[INDEX], FUN, ...)
		colnames(out)[ncol(out)] <- paste(X@name, "Result", sep = "")
	#} else if (use_plyr) {
	#	if (is.null(INDEX)) stop("Unable to determine index name, please use character input.")
	#	out <- ddply(
	#		X@output, 
	#		INDEX,
	#		oapplyf(FUN),
	#		rescol = paste(X@name, "Result", sep = ""),
	#		datvars = var, 
	#		...,
	#		.drop = TRUE
	#	)
	} else {
		# Old implementation
		out <- tapply(
			X = out[[paste(X@name, "Result", sep = "")]], 
			INDEX = X@output[INDEX],
			FUN = FUN,
			...,
			simplify = simplify
		)
		if (length(out) == 0) stop("0 length return from tapply!\n")
		if (is.list(out[1])) { # function returned array
			rows <- tapply(1:nrow(X@output), X@output[INDEX], I)
			if (length(dim(out[[1]])) == 2) {
				out <- lapply(out, t)
			}
			out <- lapply(lapply(out, as.table), as.data.frame)
			for (i in 1:length(out)) {
				out[[i]]$Row <- rows[[i]]
			}
			out <- do.call(rbind, out)
			
			temp <- X@output[!colnames(X@output) %in% paste(X@name, "Result", sep = "")]
			temp$Row <- 1:nrow(temp)
			
			out <- merge(temp, out)
			out <- out[colnames(out) != "Row"]
		}
		else { # function returned single value
			out <- as.data.frame(as.table(out))
		}
		nas <- is.na(out$Freq)
		if (any(nas)) {
			out <- out[!nas,]
			warning(paste(sum(nas), "NAs removed. Consider using na.rm = TRUE if this seems unusual or drop_na = FALSE if you do not want to remove NAs automatically."))
		}
		
		colnames(out)[colnames(out) == "Freq"] <- paste(X@name, "Result", sep = "")
	}
	X@output <- out
	X@marginal <- colnames(out) %in% marginals # Marginals can be easily corrected here disrequiring CheckMarginals
	return(X)
}

oapplyf <- function(fun) {
	if (is.character(fun)) fun <- get(fun)
	out <- function(dat, rescol, datvars, ...) {
		# Take first entry of each index (since they should contain only one distinct value)
		out <- data.frame(dat[[datvars[1]]][1])
		if (length(datvars > 1)) {
			for (i in 2:length(datvars)) {
				out[[i]] <- dat[[datvars[i]]][1]
			}
		}
		out <- data.frame(out, fun(dat[[rescol]], ...))
		colnames(out) <- c(datvars, rescol)
		return(out)
	}
	return(out)
}

# A memory-saving function for oapply when there is exactly one row for each unique combination.
# All non-marginal indices are removed.
ooapply <- function(
	X, # An ovariable
	cols, # Names of index columns to aggregrate over
	FUN = "sum", # A function to used in aggregation. Only "sum", "mean", "min", "max" and "prod" are available atm.
	... # For compatibility.
) {

	rescol <- paste(X@name, "Result", sep = "")
	X <- unkeep(X, # Unkeep all columns except critical marginals and the result.
		cols = setdiff(colnames(X@output)[!X@marginal], rescol),
		prevresults = TRUE,
		sources = TRUE 
	) 
	keeps <- colnames(X@output)[X@marginal & !colnames(X@output) %in% cols] # Marginals to keep
	if(any(colnames(X@output)[X@marginal] %in% cols)) {
		ro <- unique(X@output[cols[1]]) # data.frame with all combinations of marginal locations
		if(length(cols) == 1) colu <- keeps else colu <- c(cols[2:length(cols)], keeps)
		for(j in colu) {
			ro <- merge(unique(X@output[j]), ro)
		}
		ro <- ro[ncol(ro):1]
		nro <- nrow(ro)
		ro <- merge(ro, X@output, all.x = TRUE)
		ro <- ro[do.call(order, ro[1:(ncol(ro) - 1)]) , ]
		res <- ro[[rescol]] # Result column in the right order.
		if(length(res) != nro) stop("The numbers of rows don't match.\n")
		if(FUN == "prod") out <- 1
		if(FUN %in% c("sum", "mean")) out <- 0
		if(FUN == "min") out <- Inf
		if(FUN == "max") out <- -Inf
		res[is.na(res)] <- out
		block <- unique(ro[keeps]) # All combinations of locations of marginals to keep
		keepn <- nrow(block)
		if(FUN == "mean") res <- res * keepn / nro
		for(i in 1:(nro / keepn)) { # Loop across all combinations of locations of marginals not to keep
			addi <- res[((i - 1) * keepn + 1):((i - 1) * keepn + keepn)] 
			if(FUN == "prod") out <- out * addi
			if(FUN %in% c("sum", "mean")) out <- out + addi
			if(FUN == "min") out <- pmin(out, addi)
			if(FUN == "max") out <- pmax(out, addi)
		}
		out <- data.frame(block, Result = out)
		colnames(out)[colnames(out) == "Result"] <- rescol
		X@output <- out
		X@marginal <- colnames(X@output) %in% keeps
	}
	return(X)
}
