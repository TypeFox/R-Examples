`mergeSGP` <- 
function(list_1,
	list_2) {

	### Merge lists

	if (is.null(names(list_1))) return(list_2)
	if (!is.null(names(list_2))) {
		for (j in c("Coefficient_Matrices", "Cutscores", "Goodness_of_Fit", "Knots_Boundaries", "Linkages", "SGPercentiles", "SGProjections", "Simulated_SGPs", "Error_Reports")) {
			list_1[[j]] <- c(list_1[[j]], list_2[[j]])[!duplicated(names(c(list_1[[j]], list_2[[j]])))]
		}

		### SGPercentiles, SGProjections, Simulated_SGPs

		for (j in c("SGPercentiles", "SGProjections", "Simulated_SGPs")) {
			if (all(names(list_2[[j]]) %in% names(list_1[[j]]))) {
				for (k in names(list_2[[j]])) { # merging list_2 in with list_1, so use it here
					if (!identical(list_1[[j]][[k]], list_2[[j]][[k]])) { # keeps it from copying first set of results
						list_1[[j]][[k]] <- rbindlist(list(list_1[[j]][[k]], list_2[[j]][[k]]), fill=TRUE)
					}
				}
			}
		}

		### Goodness_of_Fit, Knots_Boundaries

		for (j in c("Goodness_of_Fit", "Knots_Boundaries")) {
			for (k in names(list_2[[j]])) {
				if (!identical(list_1[[j]][[k]], list_2[[j]][[k]])) {
 					names.list <- c(unique(names(list_1[[j]][[k]])), unique(names(list_2[[j]][[k]]))) # Get list of (unique) names first.
					list_1[[j]][[k]] <- c(list_1[[j]][[k]], list_2[[j]][[k]][!names(list_2[[j]][[k]]) %in% names(list_1[[j]][[k]])]) #new.elements
					if (any(duplicated(names.list))) {
						dups <- names.list[which(duplicated(names.list))]
						for (l in seq(dups)) {
							if (!identical(list_1[[j]][[k]][[dups[l]]], list_2[[j]][[k]][[dups[l]]])) { # could be same matrices, different @Version (???)
								x <- length(list_1[[j]][[k]])+1
								list_1[[j]][[k]][[x]] <- list_2[[j]][[k]][[dups[l]]]
								names(list_1[[j]][[k]]) <- c(names(list_1[[j]][[k]])[-x], dups[l])
							}
						}
					}
				}
			}
		} # j in c("Goodness_of_Fit", "Knots_Boundaries")

		### Coefficient_Matrices

		j <- "Coefficient_Matrices"
		for (k in names(list_2[[j]])) {
			if (!grepl("SIMEX", k)) {
				if (!identical(list_1[[j]][[k]], list_2[[j]][[k]])) {
					list_1[[j]][[k]] <- unique.splineMatrix(c(list_1[[j]][[k]], list_2[[j]][[k]]))
				}
			} else {
				for (grd_ord in names(list_2[[j]][[k]])) {
					for (lambda in names(list_2[[j]][[k]][[grd_ord]])) {
						if (!identical(list_1[[j]][[k]][[grd_ord]][[lambda]], list_2[[j]][[k]][[grd_ord]][[lambda]])) {
							list_1[[j]][[k]][[grd_ord]][[lambda]] <- unique.splineMatrix(c(list_1[[j]][[k]][[grd_ord]][[lambda]], list_2[[j]][[k]][[grd_ord]][[lambda]]))
						}
					}
				}
			}
		} # j <- "Coefficient_Matrices"
	}
	list_1[which(names(list_1) != "Panel_Data")]
} ### END mergeSGP
