`checkConfig` <-
function(my.config,
	config.type="Standard") {

	if (config.type=="Standard") {
		
		for (i in seq_along(my.config)) {
			my.config[[i]][['sgp.content.areas']] <- my.config[[i]][['sgp.content.areas']][!is.na(my.config[[i]][['sgp.content.areas']])]
			my.config[[i]][['sgp.panel.years']] <- my.config[[i]][['sgp.panel.years']][!is.na(my.config[[i]][['sgp.panel.years']])]
			my.config[[i]][['sgp.panel.years']] <- as.character(my.config[[i]][['sgp.panel.years']])
			my.config[[i]][['sgp.grade.sequences']] <- my.config[[i]][['sgp.grade.sequences']][!is.na(my.config[[i]][['sgp.grade.sequences']])]
			my.config[[i]][['sgp.grade.sequences']] <- as.character(my.config[[i]][['sgp.grade.sequences']])
			if (length(my.config[[i]][['sgp.projection.grade.sequences']]) > 0) {
				my.config[[i]][['sgp.projection.grade.sequences']] <- 
					my.config[[i]][['sgp.projection.grade.sequences']][!is.na(my.config[[i]][['sgp.projection.grade.sequences']])]
				tmp.min <- min(length(my.config[[i]][['sgp.content.areas']]), length(my.config[[i]][['sgp.projection.grade.sequences']]))
				my.config[[i]][['sgp.projection.grade.sequences']] <- tail(as.character(my.config[[i]][['sgp.projection.grade.sequences']]), tmp.min)
			} else my.config[[i]][['sgp.projection.grade.sequences']] <- NA

			if (length(my.config[[i]][['sgp.content.areas']]) != length(my.config[[i]][['sgp.grade.sequences']])) {
				tmp.min <- min(length(my.config[[i]][['sgp.content.areas']]), length(my.config[[i]][['sgp.grade.sequences']]))
				my.config[[i]][['sgp.content.areas']] <- tail( my.config[[i]][['sgp.content.areas']], tmp.min)
				my.config[[i]][['sgp.grade.sequences']] <- tail( my.config[[i]][['sgp.grade.sequences']], tmp.min)
			}
		}
		return(my.config)
	}

	if (config.type=="Baseline") {

		if (!all(unlist(sapply(lapply(my.config, names),
			function(x) x %in% c("sgp.baseline.content.areas", "sgp.baseline.panel.years", "sgp.baseline.grade.sequences", "sgp.baseline.grade.sequences.lags", "sgp.baseline.panel.years.within", "sgp.baseline.exclude.sequences"))))) {
				stop("Please specify an appropriate list of SGP function labels (sgp.baseline.config).  See help page for details.")
		}       

		for (i in seq_along(my.config)) {
			my.config[[i]][['sgp.baseline.content.areas']] <- my.config[[i]][['sgp.baseline.content.areas']][!is.na(my.config[[i]][['sgp.baseline.content.areas']])]
			my.config[[i]][['sgp.baseline.panel.years']] <- my.config[[i]][['sgp.baseline.panel.years']][!is.na(my.config[[i]][['sgp.baseline.panel.years']])]
			my.config[[i]][['sgp.baseline.grade.sequences']] <- my.config[[i]][['sgp.baseline.grade.sequences']][!is.na(my.config[[i]][['sgp.baseline.grade.sequences']])]
			my.config[[i]][['sgp.baseline.panel.years']] <- as.character(my.config[[i]][['sgp.baseline.panel.years']])
			my.config[[i]][['sgp.baseline.grade.sequences']] <- as.character(my.config[[i]][['sgp.baseline.grade.sequences']])

			if (length(my.config[[i]][['sgp.baseline.content.areas']]) != length(my.config[[i]][['sgp.baseline.grade.sequences']])) {
				tmp.min <- min(length(my.config[[i]][['sgp.baseline.content.areas']]), length(my.config[[i]][['sgp.baseline.grade.sequences']]))
				my.config[[i]][['sgp.baseline.content.areas']] <- tail( my.config[[i]][['sgp.baseline.content.areas']], tmp.min)
				my.config[[i]][['sgp.baseline.grade.sequences']] <- tail( my.config[[i]][['sgp.baseline.grade.sequences']], tmp.min)
			}
		}
		return(my.config)
	}
} ### END checkConfig
