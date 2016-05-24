`getSGPBaselineConfig` <- 
function(sgp_object, 
	content_areas,
	grades, 
	sgp.baseline.panel.years,
	sgp.percentiles.baseline.max.order,
	calculate.simex.baseline=NULL) {

	sgp.baseline.config <- tmp.sgp.baseline.config <- .content_areas <- .years <- .grades <- .sgp.grade.sequences <- list()

	if (is.null(content_areas)) {
		.content_areas <- unique(sgp_object@Data["VALID_CASE"][["CONTENT_AREA"]])
	} else {
		.content_areas <- content_areas
	}

	if (is.null(sgp.baseline.panel.years)) {
		.years <- head(sort(unique(sgp_object@Data[SJ("VALID_CASE", .content_areas)][["YEAR"]])), 5)
	} else {
		.years <- sgp.baseline.panel.years
	}

	.baseline.max.order <- min(sgp.percentiles.baseline.max.order, length(.years)-2)

	for (i in .content_areas) {
		.grades <- sort(type.convert(unique(sgp_object@Data[SJ("VALID_CASE", i)][["GRADE"]])))
		tmp.sgp.grade.sequences <- lapply(.grades[-1], function(x) tail(.grades[.grades <= x], (.baseline.max.order+1)))
		tmp.sgp.baseline.grade.sequences <- sapply(tmp.sgp.grade.sequences, function(x) x[(tail(x,1)-x) <= length(.years)-2])
		if (!is.null(grades)) {
			tmp.sgp.baseline.grade.sequences <- tmp.sgp.baseline.grade.sequences[sapply(tmp.sgp.baseline.grade.sequences, function(x) tail(x, 1) %in% grades)]
		}
	
		sgp.baseline.grade.sequences <- list()
		for (a in seq_along(tmp.sgp.baseline.grade.sequences)) {
			sgp.baseline.grade.sequences[[a]] <-
				eval(parse(text=paste("list(", paste("tail(tmp.sgp.baseline.grade.sequences[[", a, "]],", length(tmp.sgp.baseline.grade.sequences[[a]]):2, ")", collapse=", "), ")")))
		}
	
		sgp.baseline.grade.sequences <- unlist(sgp.baseline.grade.sequences, recursive=FALSE)
		sgp.baseline.grade.sequences.lags <- lapply(sgp.baseline.grade.sequences, diff)

		tmp.sgp.baseline.config[[as.character(i)]] <- 
			list(
				sgp.baseline.content.areas=i, 
				sgp.baseline.panel.years=.years,
				sgp.baseline.grade.sequences=sgp.baseline.grade.sequences,
				sgp.baseline.grade.sequences.lags=sgp.baseline.grade.sequences.lags,
				sgp.baseline.calculate.simex.baseline=calculate.simex.baseline)
	}

	for (a in seq_along(tmp.sgp.baseline.config)) {
		# Set tmp.length only once to the first list length.  Overwrites elements if subsequent list length differ from the first
		if (a == 1) tmp.length <- length(tmp.sgp.baseline.config[[a]][["sgp.baseline.grade.sequences"]])
		for (b in 1:length(tmp.sgp.baseline.config[[a]][["sgp.baseline.grade.sequences"]])) {
			sgp.baseline.config[[b+(a-1)*tmp.length]] <- tmp.sgp.baseline.config[[a]]
			sgp.baseline.config[[b+(a-1)*tmp.length]][["sgp.baseline.grade.sequences"]] <- unlist(tmp.sgp.baseline.config[[a]][["sgp.baseline.grade.sequences"]][b])
			sgp.baseline.config[[b+(a-1)*tmp.length]][["sgp.baseline.grade.sequences.lags"]] <- unlist(tmp.sgp.baseline.config[[a]][["sgp.baseline.grade.sequences.lags"]][b])
			sgp.baseline.config[[b+(a-1)*tmp.length]][["sgp.baseline.content.areas"]] <- 
				rep(sgp.baseline.config[[b+(a-1)*tmp.length]][["sgp.baseline.content.areas"]], length(sgp.baseline.config[[b+(a-1)*tmp.length]][["sgp.baseline.grade.sequences"]]))
			sgp.baseline.config[[b+(a-1)*tmp.length]][["sgp.baseline.calculate.simex.baseline"]] <- tmp.sgp.baseline.config[[a]][["sgp.baseline.calculate.simex.baseline"]]
			if ("YEAR_WITHIN" %in% names(sgp_object@Data)) {
				sgp.baseline.config[[b+(a-1)*tmp.length]][['sgp.baseline.panel.years.within']] <- rep("LAST_OBSERVATION", length(content_areas))
			}
		}
	}

	checkConfig(sgp.baseline.config, "Baseline")
} ## END getSGPBaselineConfig
