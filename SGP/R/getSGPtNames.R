`getSGPtNames` <- 
function(sgp.iter,
	SGPt=NULL,
	type) {

	if (is.null(SGPt)) {
		return(NULL)
	}

	if (type %in% c("sgp.percentiles", "sgp.percentiles.equated", "sgp.percentiles.baseline")) {
		tmp.list <- as.list(paste(c("TIME", "TIME_LAG"), tail(sgp.iter[["sgp.panel.years"]], 1), tail(sgp.iter[["sgp.content.areas"]], 1), sep="."))
		names(tmp.list) <- c("TIME", "TIME_LAG")
		return(tmp.list)
	}

	if (type=="sgp.projections") {
		return(paste("DATE", tail(sgp.iter[["sgp.projection.panel.years"]], 1), tail(sgp.iter[["sgp.projection.content.areas"]], 1), sep="."))
	}

	if (type=="sgp.projections.baseline") {
		return(paste("DATE", tail(sgp.iter[["sgp.projection.baseline.panel.years"]], 1), tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1), sep="."))
	}

	if (type %in% c("sgp.projections.lagged", "sgp.projections.lagged.baseline")) {
		return(paste("DATE", tail(sgp.iter[["sgp.panel.years"]], 1), tail(sgp.iter[["sgp.content.areas"]], 1), sep="."))
	}
} ## END getSGPtNames
