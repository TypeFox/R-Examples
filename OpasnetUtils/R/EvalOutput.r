# EvalOutput #################### evaluates the output slot of ovariables
#### Runs CheckMarginals as well

EvalOutput <- function(variable, fillna = FALSE, indent = 0, verbose = FALSE, ...) { # ... for e.g na.rm 
	if (verbose) {
		cat(rep("-", indent), "Evaluating", variable@name, "...")
		t1 <- Sys.time()
	}
	ComputeDependencies(variable@dependencies, fillna = fillna, indent = indent + 1, verbose = verbose, new_code = TRUE, ...)
	variable <- ddata_apply(variable, ...)
	if (nrow(variable@data) > 0) {
		colnames(variable@data)[colnames(variable@data) %in% "Result"] <- paste(variable@name, "Result", sep = "")
		rescol <- paste(variable@name, "Result", sep = "")
		if (!is.numeric(variable@data[[rescol]]) & !is.null(variable@data[[rescol]])) {
			a <- interpret(variable@data, rescol = rescol, ...) 
		} else a <- variable@data
	} else a <- variable@data
	b <- variable@formula(variable@dependencies, indent = indent, verbose = verbose, ...)
	tempmarginals <- character()
	if (class(b)[1]=="ovariable") {
		if (length(b@marginal) > 0) {
			tempmarginals <- c(
				tempmarginals, 
				colnames(b@output)[b@marginal], 
				paste(variable@name, "Source", sep = "") # CheckMarginal expects complete a marginal if it exists at all
			) 
		}
		b <- b@output
	}
	if (is.numeric(b) & nrow(a) == 0) {
		if (verbose) cat("\n")
		stop(paste("No proper data nor formula defined for ", variable@name, "! (Numeric formula return and 0 rows data)\n", sep = ""))
	}
	if (is.numeric(b)) {
		colnames(a)[colnames(a) == rescol] <- paste(variable@name, "Result", sep = "")
		a[,paste(variable@name, "Source", sep = "")] <- "Data"
		variable@output <- a
		if (verbose) {
			td <- Sys.time() - t1
			cat(paste(" done(", round(td, 2), " ", attributes(td)$units, ")!\n", sep = ""))
		}
	}
	else if (nrow(a) == 0) {
		colnames(b)[
			colnames(b) %in% "Result"
		] <- paste(variable@name, "Result", sep = "")
		b[,paste(variable@name, "Source", sep = "")] <- "Formula"
		variable@output <- b
		if (length(tempmarginals) > 1) variable@marginal <- colnames(variable@output) %in% tempmarginals
		if (verbose) {
			td <- Sys.time() - t1
			cat(paste(paste(rep("-", indent), collapse = ""), " done(", round(td, 2), " ", attributes(td)$units, ")!\n", sep = ""))
		}
	}
	else {
		colnames(a)[colnames(a) == rescol] <- "FromData"
		colnames(b)[colnames(b) %in% c(paste(variable@name, "Result", sep = ""), "Result")] <- "FromFormula" # *
		# <variablename> prefix not necessitated for "Result" column of formula output
		temp <- melt(
			merge(a, b, all = TRUE, ...), # Will cause problems if dependencies contain non-marginal indices that match with -
			# marginal indeces in data. Or maybe not.
			measure.vars = c("FromData", "FromFormula"),
			variable.name = paste(variable@name, "Source", sep = ""),
			value.name = paste(variable@name, "Result", sep = ""),
			...
		)
		levels(
			temp[[paste(variable@name, "Source", sep = "")]]
		) <- gsub("^From", "", 
			levels(
				temp[[paste(variable@name, "Source", sep = "")]]
			)
		)
		variable@output <- temp
		if (length(tempmarginals) > 1) variable@marginal <- colnames(variable@output) %in% tempmarginals
		if (verbose) {
			td <- Sys.time() - t1
			cat(paste(paste(rep("-", indent), collapse = ""), " done(", round(td, 2), " ", attributes(td)$units, ")!\n", sep = ""))
		}
	}
	#if (verbose) cat(rep("-", indent), " done!\n")
	#if (verbose) cat(" done!\n")
	variable <- CheckMarginals(variable, indent = indent, verbose = verbose, ...)
	if (fillna) {
		ret <- tryCatch(variable@output <-  fillna(variable@output, 1:ncol(variable@output)[variable@marginal]), error = function(e) return(NULL))
		if (is.null(ret)) warning("Unable to FillNA.")
	}
	return(variable)
}