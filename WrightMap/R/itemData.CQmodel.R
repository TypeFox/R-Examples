itemData.CQmodel <- function(thresholds, item.table = NULL, interactions = NULL, step.table = NULL, item.type = "default", throld = 0.5, 
	...) {

	unpack.GIN <- function(GIN) {
		if (class(GIN) == "matrix") 
			return(GIN)
		else {
			return(do.call(cbind, lapply(GIN, unpack.GIN)))
		}
	}

	unpack.names <- function(GIN, sofar = "") {
		if (class(GIN) == "matrix") {
			my.names <- c(1:ncol(GIN))
		} else my.names <- names(GIN)


		if (length(my.names) == 1) {
			names <- sofar
		} else if (length(sofar) == 1) {
			names <- my.names
		} else names <- c(outer(my.names, sofar, paste))

		if (class(GIN) == "matrix") 
			return(names)
		else return(unpack.names(GIN[[1]], names))
	}

	model <- thresholds

	if (throld == 0.5 && !is.null(model$GIN) && is.null(item.table) && (item.type != "deltas")) {
		#print("false")
		throlds <- unpack.GIN(model$GIN)
		names <- unpack.names(model$GIN)
		colnames(throlds) <- names

		message("Using GIN table for threshold parameters")
	} else {
		RMP <- model$RMP

		if (item.type != "thresholds") {
			throlds <- make.deltas(model, item.table = item.table, interactions = interactions, step.table = step.table)

		} else {
			throlds <- make.thresholds(model, item.table = item.table, interactions = interactions, step.table = step.table, throld = throld)


		}

	}
	
	return(throlds)
}
