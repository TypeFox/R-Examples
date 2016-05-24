#######################################
# CollapseMarginal
##########################################
# Collapses marginals by applying sums, means or samples
# Also loses all non-marginal columns except the relevant Result
# Parse-able table should have columns Index, Function and Probs
# Probs can be left out for equal weights sampling
# If Function is not given mean is assumed
#####################################

CollapseTableParser <- function(CTable, env = .GlobalEnv){ # CTable is a data.frame
	for (i in unique(as.character(CTable$Variable))) {
		temp <- CTable[CTable$Variable == i,]
		cols <- temp[["Index"]]
		probs <- strsplit(as.character(temp[["Probs"]]), ",")
		probs <- lapply(probs, as.numeric)
		fun <- temp[["Function"]]
		out <- list(cols = cols, probs = probs, fun = fun)
		assign(paste("Col", i, sep = ""), out, envir = env)
	}
}

CheckCollapse <- function(variable, indent = 0, verbose = TRUE, ...) {
	if (exists(paste("Col", variable@name, sep = ""))) {
		if (verbose) cat(rep("-", indent), "Processing", variable@name, "marginal collapses", "...")
		Col <- get(paste("Col", variable@name, sep = ""))
		variable <- CollapseMarginal(variable, Col$cols, Col$fun, Col$probs, ...)
		if (verbose) cat(" done!\n")
	}
	return(variable)
}

CollapseMarginal <- function(variable, cols, fun = "mean", probs = NULL, ...) { # cols is a character vector, while probs is a list
	if (length(fun) == 0) fun <- "mean"
	# If no probabilities given use function
	# Also if given funtion is sample then equal weights are used and this section is skipped
	if (is.na(fun)) stop("No function given to collapse with!\n")
	if (length(probs) == 0 & fun != "sample") {
		fun <- get(fun)
		out <- oapply(variable, FUN = fun, cols = cols, na.rm = TRUE)
		return(out)
	}
	
	# Else use sample with option of given probabilities
	out <- variable@output
	marginals <- colnames(out)[variable@marginal]
	
	if (!is.list(probs) & is.numeric(probs)) probs <- list(probs)
	if (!is.null(probs) & length(probs) != length(cols)) stop("Number of columns does not match number of probability vectors given!\n")
	if ("Iter" %in% colnames(out)) {
		N <- max(out$Iter)
	} else {
		N <- get("N", openv)
	}
	for (i in 1:length(cols)) {
		b <- probs[[i]]
		locs <- levels(out[[cols[i]]])
		if (is.null(b)) b <- rep(1, length(locs))
		if (any(is.na(b))) b <- rep(1, length(locs)) # dont see why NA would turn up here, but hey lets just be sure
		if (length(b) != length(locs)) {
			stop(paste("Number of locations does not match number of probabilities given for ", cols[i], "!\n", sep = ""))
		}
		selection <- data.frame(
			Iter = 1:N, 
			sample(
				locs, 
				size = N, 
				replace = TRUE, 
				prob = b
			)
		)
		colnames(selection)[2] <- cols[i]
		out <- merge(selection, out)
	}
	variable@output <- out
	variable@marginal <- colnames(out) %in% c(marginals, "Iter") & ! colnames(out) %in% cols
	return(variable)
}