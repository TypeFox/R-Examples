# Fetch used to be named Fetch2
Fetch2 <- function(...) Fetch(...)

# Fetch #################### loads all given dependencies to global memory

Fetch <- function(dependencies, evaluate = FALSE, indent = 0, verbose = TRUE, ...) {
	if (nrow(dependencies) > 0) {
		for (i in 1:nrow(dependencies)) {
			if(!exists(as.character(dependencies$Name[i]))) {
				testkey <- if (is.null(dependencies$Key[i])) TRUE else is.na(dependencies$Key[i]) | dependencies$Key[i] == ""
				testid <- if (is.null(dependencies$Ident[i])) TRUE else is.na(dependencies$Ident[i]) | dependencies$Ident[i] == ""
				if (testkey & testid) {
					stop(paste("No key nor ident given for dependent variable: ", dependencies$Name[i], "!", sep = ""))
				}
				if (!testkey) {
					objects.get(dependencies$Key[i]) # Key is the R-tools session identifier (shown at the end of the url)
				}
				if (testkey & !testid) {
					ident <- strsplit(as.character(dependencies$Ident[i]), "/")[[1]] # Ident should be in format <page_id>/<code_name>
					objects.latest(ident[1], ident[2])
				}
				if (evaluate) assign(
						as.character(dependencies$Name[i]), 
						EvalOutput(get(as.character(dependencies$Name[i])), ...), 
						envir = as.environment(find(as.character(dependencies$Name[i])))
				)
				if (verbose) cat("\n", rep("-", indent), as.character(dependencies$Name[i]), "fetched successfully!\n")
			}
		}
	}
} # no need to return anything since variables are written in global memory by objects.get