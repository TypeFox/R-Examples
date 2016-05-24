`getPercentileTableNames` <- 
function(sgp_object,
	content_areas=NULL,
	state=NULL,
	years=NULL,
	sgp.type,
	sgp.percentiles.equated=FALSE,
	use.cohort.for.baseline.when.missing=NULL) {

        if (is.null(use.cohort.for.baseline.when.missing)) {
                if (!is.null(state) && is.null(SGP::SGPstateData[[state]][["SGP_Configuration"]][["use.cohort.for.baseline.when.missing"]])) {
                        use.cohort.for.baseline.when.missing <- FALSE
                } 
                if (!is.null(state) && !is.null(SGP::SGPstateData[[state]][["SGP_Configuration"]][["use.cohort.for.baseline.when.missing"]])) {
                        use.cohort.for.baseline.when.missing <- SGP::SGPstateData[[state]][["SGP_Configuration"]][["use.cohort.for.baseline.when.missing"]]
                }
        }

	if (sgp.type %in% c("sgp.percentiles", "sgp.percentiles.baseline")) {
		tmp.sgp.names <- as.character(names(sgp_object@SGP$SGPercentiles))
		tmp.baseline.names <- grep("BASELINE", names(sgp_object@SGP$SGPercentiles), value=TRUE)
		if (sgp.type=="sgp.percentiles") tmp.names <- setdiff(tmp.sgp.names, tmp.baseline.names)
		if (sgp.type=="sgp.percentiles.baseline") {
			tmp.names <- tmp.baseline.names
			if (use.cohort.for.baseline.when.missing) {
				tmp.content_areas.diff <- setdiff(unique(sapply(strsplit(tmp.sgp.names, "[.]"), function(x) paste(x[1:2], collapse="."))), 
					unique(sapply(strsplit(tmp.baseline.names, "[.]"), function(x) paste(x[1:2], collapse="."))))
				if (length(tmp.content_areas.diff) > 0) {
					if (!is.null(years)) tmp.content_areas.diff <- tmp.content_areas.diff[sapply(strsplit(tmp.content_areas.diff, "[.]"), function(x) x[2] %in% years)]
					if (!is.null(content_areas)) tmp.content_areas.diff <- tmp.content_areas.diff[sapply(strsplit(tmp.content_areas.diff, "[.]"), function(x) x[2] %in% content_areas)]
				}
				if (length(tmp.content_areas.diff) > 0) {
					tmp.names <- c(tmp.names, unlist(lapply(tmp.content_areas.diff, function(x) tmp.sgp.names[grep(x, tmp.sgp.names)])))
					message(c("\tNOTE: Cohort referenced SGPs being used for baseline referenced SGPs for content areas and years:\n\t\t", 
						paste(unlist(lapply(tmp.content_areas.diff, function(x) tmp.sgp.names[grep(x, tmp.sgp.names)])), collapse=",\n\t\t")))
				}
			} ### END if (use.cohort.for.baseline.when.missing)
		}
		if (length(tmp.names) > 0 & !is.null(years)) tmp.names <- tmp.names[sapply(tmp.names, function(x) getTableNameYear(x) %in% years)]
		if (length(tmp.names) > 0 & !is.null(content_areas)) tmp.names <- tmp.names[sapply(strsplit(tmp.names, "[.]"), function(x) x[1] %in% content_areas)]
		if (sgp.percentiles.equated) tmp.names <- grep("EQUATED", tmp.names, value=TRUE, invert=TRUE)
		return(tmp.names)
	}

	if (sgp.type %in% c("sgp.projections", "sgp.projections.baseline")) {
		tmp.target.names <- grep("TARGET_SCALE_SCORES", names(sgp_object@SGP$SGProjections), value=TRUE)
		tmp.lagged.names <- grep("LAGGED", names(sgp_object@SGP$SGProjections), value=TRUE)
		tmp.sgp.projections.names <- setdiff(names(sgp_object@SGP$SGProjections), c(tmp.target.names, tmp.lagged.names))
		if (!is.null(tmp.sgp.projections.names)) {
		tmp.baseline.names <- grep("BASELINE", tmp.sgp.projections.names, value=TRUE)
		if (sgp.type=="sgp.projections") tmp.names <- setdiff(tmp.sgp.projections.names, tmp.baseline.names)
		if (sgp.type=="sgp.projections.baseline") {
			tmp.names <- tmp.baseline.names
			if (use.cohort.for.baseline.when.missing) {
				tmp.content_areas.diff <- setdiff(unique(sapply(strsplit(tmp.sgp.projections.names, "[.]"), function(x) paste(x[1:2], collapse="."))), 
					unique(sapply(strsplit(tmp.baseline.names, "[.]"), function(x) paste(x[1:2], collapse="."))))
				if (length(tmp.content_areas.diff) > 0) {
					if (!is.null(years)) tmp.content_areas.diff <- tmp.content_areas.diff[sapply(strsplit(tmp.content_areas.diff, "[.]"), function(x) x[2] %in% years)]
					if (!is.null(content_areas)) tmp.content_areas.diff <- tmp.content_areas.diff[sapply(strsplit(tmp.content_areas.diff, "[.]"), function(x) x[2] %in% content_areas)]
				}
				if (length(tmp.content_areas.diff) > 0) {
					tmp.names <- c(tmp.names, tmp.content_areas.diff)
					message(c("\tNOTE: Cohort referenced projections being used for baseline referenced projections for content areas and years: ", paste(tmp.content_areas.diff, collapse=", ")))
				}
			} ### END if (use.cohort.for.baseline.when.missing)
		}
		if (length(tmp.names) > 0) tmp.names <- tmp.names[sapply(lapply(tmp.names, function(x) grep("LEVEL", names(sgp_object@SGP[['SGProjections']][[x]]))), function(x) length(x)>0)]
		if (length(tmp.names) > 0 & !is.null(years)) tmp.names <- tmp.names[sapply(tmp.names, function(x) getTableNameYear(x) %in% years)]
		if (length(tmp.names) > 0 & !is.null(content_areas)) tmp.names <- tmp.names[sapply(strsplit(tmp.names, "[.]"), function(x) x[1] %in% content_areas)]
		return(tmp.names)
		} else return(NULL)
	}

	if (sgp.type %in% c("sgp.projections.lagged", "sgp.projections.lagged.baseline")) {
		tmp.target.names <- grep("TARGET_SCALE_SCORES", names(sgp_object@SGP$SGProjections), value=TRUE)
		tmp.lagged.names <- setdiff(grep("LAGGED", names(sgp_object@SGP$SGProjections), value=TRUE), tmp.target.names)
		tmp.baseline.names <- intersect(tmp.lagged.names, grep("BASELINE", names(sgp_object@SGP$SGProjections), value=TRUE))
		if (sgp.type=="sgp.projections.lagged") tmp.names <- setdiff(tmp.lagged.names, tmp.baseline.names)
		if (sgp.type=="sgp.projections.lagged.baseline") {
			tmp.names <- tmp.baseline.names
			if (use.cohort.for.baseline.when.missing) {
				tmp.content_areas.diff <- setdiff(unique(sapply(strsplit(tmp.lagged.names, "[.]"), function(x) paste(x[1:2], collapse="."))), 
					unique(sapply(strsplit(tmp.baseline.names, "[.]"), function(x) paste(x[1:2], collapse="."))))
				if (length(tmp.content_areas.diff) > 0) {
					if (!is.null(years)) tmp.content_areas.diff <- tmp.content_areas.diff[sapply(strsplit(tmp.content_areas.diff, "[.]"), function(x) x[2] %in% years)]
					if (!is.null(content_areas)) tmp.content_areas.diff <- tmp.content_areas.diff[sapply(strsplit(tmp.content_areas.diff, "[.]"), function(x) x[2] %in% content_areas)]
				}
				if (length(tmp.content_areas.diff) > 0) {
					tmp.names <- c(tmp.names, unlist(lapply(tmp.content_areas.diff, function(x) tmp.lagged.names[grep(x, tmp.lagged.names)])))
					message(c("\tNOTE: Cohort referenced lagged.projections being used for baseline referenced lagged projections for content areas and years: ", paste(unlist(lapply(tmp.content_areas.diff, function(x) tmp.lagged.names[grep(x, tmp.lagged.names)])), collapse=", ")))
				}
			} ### END if (use.cohort.for.baseline.when.missing)
		}
		if (length(tmp.names) > 0) tmp.names <- tmp.names[sapply(lapply(tmp.names, function(x) grep("LEVEL", names(sgp_object@SGP[['SGProjections']][[x]]))), function(x) length(x)>0)]
		if (length(tmp.names) > 0 & !is.null(years)) tmp.names <- tmp.names[sapply(tmp.names, function(x) getTableNameYear(x) %in% years)]
		if (length(tmp.names) > 0 & !is.null(content_areas)) tmp.names <- tmp.names[sapply(strsplit(tmp.names, "[.]"), function(x) x[1] %in% content_areas)]
		return(tmp.names)
	}
} ## END getPercentileTableNames
