## Returns a factor from one variable or a list of variable
## some code taken from the tapply function

group <- function (INDEX, factor=TRUE) {

	if (!is.factor(INDEX)) {
		if (!is.list(INDEX)) { INDEX <- list(INDEX) }

		il <- sapply(INDEX, length)

		if (min(il)!=max(il)) {
			stop(" [!] all factors must have the same length")
		}

		if (factor) {
			group <- rep("", length(INDEX[[1]]))
			for (i in seq_along(INDEX)) {
		        	index <- as.factor(INDEX[[i]])
		        	group <- paste(group, index)
				group[is.na(index)] <- NA
			}
			group <- as.factor(group)
		}
		else {
			## fromt apply
			nI <- length(INDEX)

			namelist <- vector("list", nI)
			names(namelist) <- names(INDEX)
			extent <- integer(nI)
			one <- 1L

			group <- rep.int(one, length(INDEX[[1]]))
			ngroup <- one

			for (i in seq_along(INDEX)) {
		        	index <- as.factor(INDEX[[i]])
		        	namelist[[i]] <- levels(index)
		        	extent[i] <- nlevels(index)
		        	group <- group + ngroup * (as.integer(index) - one)
		        	ngroup <- ngroup * nlevels(index)
			}
		}
	} else {
		## Eliminate the unused levels
		uf <- unique(INDEX)
		fl <- levels(INDEX)[levels(INDEX) %in% uf]
		group <- factor(INDEX, levels=fl)
	}

	return(group)
}

