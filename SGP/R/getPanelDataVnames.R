`getPanelDataVnames` <- 
function(sgp.type,
	sgp.iter,
	sgp.data.names,
	scale.score.variable.name=NULL) {

	if (is.null(scale.score.variable.name)) scale.score.variable.name <- "SCALE_SCORE"

	if (sgp.type=="sgp.percentiles") {
		if ("YEAR_WITHIN" %in% sgp.data.names) {
			return(c("ID", paste("GRADE", tail(sgp.iter[["sgp.panel.years"]], length(sgp.iter[["sgp.grade.sequences"]])), 
				tail(sgp.iter[["sgp.content.areas"]], length(sgp.iter[["sgp.grade.sequences"]])), 
				tail(sgp.iter[["sgp.panel.years.within"]], length(sgp.iter[["sgp.grade.sequences"]])), sep="."), 
				paste(scale.score.variable.name, tail(sgp.iter[["sgp.panel.years"]], length(sgp.iter[["sgp.grade.sequences"]])),
				tail(sgp.iter[["sgp.content.areas"]], length(sgp.iter[["sgp.grade.sequences"]])), 
				tail(sgp.iter[["sgp.panel.years.within"]], length(sgp.iter[["sgp.grade.sequences"]])), sep=".")))
		}

		return(c("ID", paste("GRADE", tail(sgp.iter[["sgp.panel.years"]], length(sgp.iter[["sgp.grade.sequences"]])), 
			tail(sgp.iter[["sgp.content.areas"]], length(sgp.iter[["sgp.grade.sequences"]])), sep="."), 
			paste(scale.score.variable.name, tail(sgp.iter[["sgp.panel.years"]], length(sgp.iter[["sgp.grade.sequences"]])),
			tail(sgp.iter[["sgp.content.areas"]], length(sgp.iter[["sgp.grade.sequences"]])), sep=".")))
	}

	if (sgp.type=="sgp.percentiles.baseline") {
		if ("YEAR_WITHIN" %in% sgp.data.names) {
			return(c("ID", paste("GRADE", tail(sgp.iter[["sgp.panel.years"]], length(sgp.iter[["sgp.baseline.grade.sequences"]])), 
				tail(sgp.iter[["sgp.content.areas"]], length(sgp.iter[["sgp.baseline.grade.sequences"]])), 
				tail(sgp.iter[["sgp.panel.years.within"]], length(sgp.iter[["sgp.baseline.grade.sequences"]])), sep="."), 
				paste(scale.score.variable.name, tail(sgp.iter[["sgp.panel.years"]], length(sgp.iter[["sgp.baseline.grade.sequences"]])),
				tail(sgp.iter[["sgp.content.areas"]], length(sgp.iter[["sgp.baseline.grade.sequences"]])), 
				tail(sgp.iter[["sgp.panel.years.within"]], length(sgp.iter[["sgp.baseline.grade.sequences"]])), sep=".")))
		}

		return(c("ID", paste("GRADE", tail(sgp.iter[["sgp.panel.years"]], length(sgp.iter[["sgp.baseline.grade.sequences"]])), 
			tail(sgp.iter[["sgp.content.areas"]], length(sgp.iter[["sgp.baseline.grade.sequences"]])), sep="."), 
			paste(scale.score.variable.name, tail(sgp.iter[["sgp.panel.years"]], length(sgp.iter[["sgp.baseline.grade.sequences"]])),
			tail(sgp.iter[["sgp.content.areas"]], length(sgp.iter[["sgp.baseline.grade.sequences"]])), sep=".")))
	}

	if (sgp.type=="sgp.percentiles.to.return") {
		special.names.to.return <- ("YEAR_WITHIN")
		if (any(special.names.to.return %in% sgp.data.names)) {
			tmp.list <- list()
			for (i in intersect(special.names.to.return, sgp.data.names)) {
				tmp.list <- c(tmp.list, list(TMP=i))
				names(tmp.list)[[length(tmp.list)]] <- 
					paste(i, tail(sgp.iter[["sgp.panel.years"]], 1), tail(sgp.iter[["sgp.content.areas"]], 1), tail(sgp.iter[['sgp.panel.years.within']], 1), sep=".")
			}
			return(tmp.list)
		} else {
			return(NULL)
		}
	}

	if (sgp.type=="sgp.projections") {
		if ("YEAR_WITHIN" %in% sgp.data.names) {
			return(c("ID", paste("GRADE", tail(sgp.iter[["sgp.projection.panel.years"]], length(sgp.iter[["sgp.projection.grade.sequences"]])), 
				tail(sgp.iter[["sgp.projection.content.areas"]], length(sgp.iter[["sgp.projection.grade.sequences"]])), 
				head(sgp.iter[["sgp.panel.years.within"]], length(sgp.iter[["sgp.projection.grade.sequences"]])), sep="."), 
				paste(scale.score.variable.name, tail(sgp.iter[["sgp.projection.panel.years"]], length(sgp.iter[["sgp.projection.grade.sequences"]])),
				tail(sgp.iter[["sgp.projection.content.areas"]], length(sgp.iter[["sgp.projection.grade.sequences"]])), 
				head(sgp.iter[["sgp.panel.years.within"]], length(sgp.iter[["sgp.projection.grade.sequences"]])), sep=".")))
		} else {
			return(c("ID", paste("GRADE", tail(sgp.iter[["sgp.projection.panel.years"]], length(sgp.iter[["sgp.projection.grade.sequences"]])), 
				tail(sgp.iter[["sgp.projection.content.areas"]], length(sgp.iter[["sgp.projection.grade.sequences"]])), sep="."), 
				paste(scale.score.variable.name, tail(sgp.iter[["sgp.projection.panel.years"]], length(sgp.iter[["sgp.projection.grade.sequences"]])),
				tail(sgp.iter[["sgp.projection.content.areas"]], length(sgp.iter[["sgp.projection.grade.sequences"]])), sep=".")))
		}
	}

	if (sgp.type=="sgp.projections.baseline") {
		if ("YEAR_WITHIN" %in% sgp.data.names) {
			return(c("ID", paste("GRADE", tail(sgp.iter[["sgp.projection.baseline.panel.years"]], length(sgp.iter[["sgp.projection.baseline.grade.sequences"]])), 
				tail(sgp.iter[["sgp.projection.baseline.content.areas"]], length(sgp.iter[["sgp.projection.baseline.grade.sequences"]])), 
				head(sgp.iter[["sgp.panel.years.within"]], length(sgp.iter[["sgp.projection.baseline.grade.sequences"]])), sep="."), 
				paste(scale.score.variable.name, tail(sgp.iter[["sgp.projection.baseline.panel.years"]], length(sgp.iter[["sgp.projection.baseline.grade.sequences"]])),
				tail(sgp.iter[["sgp.projection.baseline.content.areas"]], length(sgp.iter[["sgp.projection.baseline.grade.sequences"]])), 
				head(sgp.iter[["sgp.panel.years.within"]], length(sgp.iter[["sgp.projection.baseline.grade.sequences"]])), sep=".")))
		} else {
			return(c("ID", paste("GRADE", tail(sgp.iter[["sgp.projection.baseline.panel.years"]], length(sgp.iter[["sgp.projection.baseline.grade.sequences"]])), 
				tail(sgp.iter[["sgp.projection.baseline.content.areas"]], length(sgp.iter[["sgp.projection.baseline.grade.sequences"]])), sep="."), 
				paste(scale.score.variable.name, tail(sgp.iter[["sgp.projection.baseline.panel.years"]], length(sgp.iter[["sgp.projection.baseline.grade.sequences"]])),
				tail(sgp.iter[["sgp.projection.baseline.content.areas"]], length(sgp.iter[["sgp.projection.baseline.grade.sequences"]])), sep=".")))
		}
	}

	if (sgp.type=="sgp.projections.lagged") {
		if ("YEAR_WITHIN" %in% sgp.data.names) {
			return(c("ID", paste("GRADE", head(tail(sgp.iter[["sgp.panel.years"]], length(sgp.iter[["sgp.grade.sequences"]])), -1), 
				head(tail(sgp.iter[["sgp.content.areas"]], length(sgp.iter[["sgp.grade.sequences"]])), -1), 
				head(sgp.iter[["sgp.panel.years.within"]], length(sgp.iter[["sgp.projection.grade.sequences"]])), sep="."), 
				paste(scale.score.variable.name, head(tail(sgp.iter[["sgp.panel.years"]], length(sgp.iter[["sgp.grade.sequences"]])), -1), 
				head(tail(sgp.iter[["sgp.content.areas"]], length(sgp.iter[["sgp.grade.sequences"]])), -1), 
				head(sgp.iter[["sgp.panel.years.within"]], length(sgp.iter[["sgp.projection.grade.sequences"]])), sep=".")))
		} else {
			return(c("ID", paste("GRADE", head(tail(sgp.iter[["sgp.panel.years"]], length(sgp.iter[["sgp.grade.sequences"]])), -1), 
				head(tail(sgp.iter[["sgp.content.areas"]], length(sgp.iter[["sgp.grade.sequences"]])), -1), sep="."), 
				paste(scale.score.variable.name, head(tail(sgp.iter[["sgp.panel.years"]], length(sgp.iter[["sgp.grade.sequences"]])), -1), 
				head(tail(sgp.iter[["sgp.content.areas"]], length(sgp.iter[["sgp.grade.sequences"]])), -1), sep=".")))
		}
	}
} ## END getPanelDataVnames
