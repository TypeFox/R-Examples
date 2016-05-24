contrastCoefficients <- function(..., contrast.definitions, data=parent.frame(), normalize=FALSE){
	# Get contrasts definitions from contrast.definitions or ...
	if (missing(contrast.definitions)) contrast.definitions <- list(...)
	# Get levels from data or the parent.frame()
	factorlevels <- getXLevels(as.list(data))
	# Preapre return list (by default with the assigned names)
	coefvectors <- vector("list", length(contrast.definitions))
	names(coefvectors) <- names(contrast.definitions)
	# Create a vector for each contrast defined as a formula
	# Leave the rest unchanged
	ignored <- seq(along=contrast.definitions) # Factors that wil be ignored (all by default)
	for (f in seq(along=contrast.definitions)){
		# Work with valid formulas of factors that exist
		if (is.call(contrast.definitions[[f]]) && (contrast.definitions[[f]][1] == "~"())){
			if (length(contrast.definitions[[f]]) < 3) stop("Malformed formula of factor contrasts.")
			names(coefvectors)[f] <- fname <- as.character(contrast.definitions[[f]][2])
			if (fname %in% names(factorlevels)){
				ignored <- ignored[ignored != f]
				# Code factor levels as k-th dimensional dummy variables in f_envir
				f_envir <- as.data.frame(diag(length(factorlevels[[fname]])))
				names(f_envir) <- factorlevels[[fname]]
				coefvectors[[f]] <- eval(languageEl(contrast.definitions[[f]], 3), envir=f_envir)
				names(coefvectors[[f]]) <- factorlevels[[fname]]
				# Normalize value if requested and suitable
				if (normalize && is.numeric(coefvectors[[f]])){
					coefvectors[[f]] <- coefvectors[[f]] / sqrt(sum(coefvectors[[f]]^2))
				}
			}
		}
	}
	# Delete ignored, and group vectors of the same name
	coefvectors[ignored] <- NULL
	coefmatrices <- sapply(unique(names(coefvectors)),
		function(n) as.matrix(as.data.frame(coefvectors[names(coefvectors)==n])),
		simplify=FALSE, USE.NAMES=TRUE)
	# Return combined list
	c(coefmatrices, contrast.definitions[ignored])
}
