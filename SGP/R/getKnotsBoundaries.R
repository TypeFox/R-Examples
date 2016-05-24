`getKnotsBoundaries` <- 
function(sgp.iter,
	state,
	sgp.iter.type="sgp.percentiles",
	sgp.year.baseline=NULL) {

	kb <- list()


	### Identify content areas needs for knots/boundaries

	if (sgp.iter.type=="sgp.percentiles") {
		tmp.content_areas <- unique(sgp.iter[['sgp.content.areas']])
	}

	if (sgp.iter.type=="sgp.percentiles.baseline") {
		tmp.content_areas <- unique(sgp.iter[['sgp.baseline.content.areas']])
	}

	if (sgp.iter.type=="sgp.projections") {
		if (identical(sgp.iter[['sgp.projection.sequence']] %in% names(SGP::SGPstateData[[state]][['SGP_Configuration']][['content_area.projection.sequence']]), TRUE)) {
			tmp.content_areas <- unique(unlist(sapply(sgp.iter[['sgp.projection.sequence']], 
				function(x) unique(SGP::SGPstateData[[state]][['SGP_Configuration']][['content_area.projection.sequence']][[x]]))))
		} else {
			tmp.content_areas <- unique(sgp.iter[['sgp.projection.content.areas']])
		}
	}

	if (sgp.iter.type=="sgp.projections.lagged") {
		if (identical(sgp.iter[['sgp.projection.sequence']] %in% names(SGP::SGPstateData[[state]][['SGP_Configuration']][['content_area.projection.sequence']]), TRUE)) {
			tmp.content_areas <- unique(unlist(sapply(unique(sgp.iter[['sgp.projection.sequence']]), 
				function(x) unique(SGP::SGPstateData[[state]][['SGP_Configuration']][['content_area.projection.sequence']][[x]]))))
		} else {
			tmp.content_areas <- unique(sgp.iter[['sgp.projection.content.areas']])
		}
	}

	if (sgp.iter.type=="sgp.projections.baseline") {
		if (identical(sgp.iter[['sgp.projection.sequence']] %in% names(SGP::SGPstateData[[state]][['SGP_Configuration']][['content_area.projection.sequence']]), TRUE)) {
			tmp.content_areas <- unique(unlist(sapply(unique(sgp.iter[['sgp.projection.sequence']]), 
				function(x) unique(SGP::SGPstateData[[state]][['SGP_Configuration']][['content_area.projection.sequence']][[x]]))))
		} else {
			tmp.content_areas <- unique(sgp.iter[['sgp.projection.baseline.content.areas']])
		}
	}

	if (sgp.iter.type=="sgp.projections.lagged.baseline") {
		if (identical(sgp.iter[['sgp.projection.sequence']] %in% names(SGP::SGPstateData[[state]][['SGP_Configuration']][['content_area.projection.sequence']]), TRUE)) {
			tmp.content_areas <- unique(unlist(sapply(unique(sgp.iter[['sgp.projection.sequence']]), 
				function(x) unique(SGP::SGPstateData[[state]][['SGP_Configuration']][['content_area.projection.sequence']][[x]]))))
		} else {
			tmp.content_areas <- unique(sgp.iter[['sgp.projection.baseline.content.areas']])
		}
	}


	### Create label(s) for Knots_Boundaries

	if (sgp.iter.type=="sgp.percentiles") {
		content_area.label <- tail(sgp.iter[["sgp.content.areas"]], 1)
		if (!is.null(sgp.year.baseline)) tmp.year <- sgp.year.baseline else tmp.year <- tail(sgp.iter[["sgp.panel.years"]], 1)
	}

	if (sgp.iter.type=="sgp.percentiles.baseline") {
		content_area.label <- tail(sgp.iter[["sgp.baseline.content.areas"]], 1)
		if (!is.null(sgp.year.baseline)) tmp.year <- sgp.year.baseline else tmp.year <- tail(sgp.iter[["sgp.panel.years"]], 1)
	}

	if (sgp.iter.type=="sgp.projections") {
		content_area.label <- tail(sgp.iter[["sgp.projection.content.areas"]], 1)
		if (!is.null(sgp.year.baseline)) tmp.year <- sgp.year.baseline else tmp.year <- tail(sgp.iter[["sgp.projection.panel.years"]], 1)
	}

	if (sgp.iter.type=="sgp.projections.lagged") {
		content_area.label <- tail(sgp.iter[["sgp.content.areas"]], 1)
		if (!is.null(sgp.year.baseline)) tmp.year <- sgp.year.baseline else tmp.year <- tail(sgp.iter[["sgp.panel.years"]], 1)
	}

	if (sgp.iter.type=="sgp.projections.baseline") {
		content_area.label <- tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1)
		if (!is.null(sgp.year.baseline)) tmp.year <- sgp.year.baseline else tmp.year <- tail(sgp.iter[["sgp.projection.baseline.panel.years"]], 1)
	}

	if (sgp.iter.type=="sgp.projections.lagged.baseline") {
		content_area.label <- tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1)
		if (!is.null(sgp.year.baseline)) tmp.year <- sgp.year.baseline else tmp.year <- tail(sgp.iter[["sgp.panel.years"]], 1)
	}

	### Create Knots_Boundaries list

	for (j in tmp.content_areas) {
		for (i in grep(j, names(SGP::SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]]), value=TRUE)) {
			kb[["Knots_Boundaries"]][[paste(content_area.label, tmp.year, sep=".")]][[i]] <- SGP::SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]][[i]]
		}
	}

	return(kb[["Knots_Boundaries"]])
} ## END getKnotsBoundaries
