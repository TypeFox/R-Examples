`piecewiseTransform` <-
function(scale_score,
	state,
	content_area,
	year,
	grade,
	output.digits=1,
	sgp.projections.equated=NULL,
	new.cutscores=NULL,
	equating.method="equipercentile") {

	if (all(is.na(scale_score))) return(scale_score)

	### Test to deal with assessment transition scenario

	if (!is.null(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]])) {
		equate.year <- SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][["Year"]]
		if (year < equate.year)  {
			tmp.test <- "Transformed_Achievement_Level_Cutscores"
		} else {
			if (!is.null(new.cutscores) && length(new.cutscores) > 0) {
				tmp.test <- "NOT_NULL"
			} else {
				tmp.test <- NULL
			}
		}
	} else {
		tmp.test <- NULL
	}


	if (is.null(sgp.projections.equated) | !is.null(tmp.test)) {
		if ((content_area %in% names(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Transformed_Achievement_Level_Cutscores"]]) &&
			grade %in% matrix(unlist(strsplit(unlist(lapply(SGP::SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]][names(SGP::SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]]) %in% content_area], names)), "_"))) || !is.null(tmp.test))) {

			if (!is.null(new.cutscores)) {
				tmp.new.cuts <- new.cutscores
			} else {
				if (!is.null(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]])) {
					tmp.new.cuts <- SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[tmp.test]][[content_area]]
				} else {
					tmp.new.cuts <- SGP::SGPstateData[[state]][["Student_Report_Information"]][["Transformed_Achievement_Level_Cutscores"]][[content_area]]
				}
			}

			my.knots_boundaries.label <- getMyLabel(state, content_area, year, "Knots_Boundaries")
			tmp.loss.hoss <- SGP::SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]][[my.knots_boundaries.label]][[paste("loss.hoss_", grade, sep="")]]
			scale_score[scale_score < tmp.loss.hoss[1]] <- tmp.loss.hoss[1]; scale_score[scale_score > tmp.loss.hoss[2]] <- tmp.loss.hoss[2]
			my.content_area <- getMyLabel(state, content_area, year)
			tmp.old.cuts <- c(tmp.loss.hoss[1], SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]][[my.content_area]][[paste("GRADE_", grade, sep="")]],
				tmp.loss.hoss[2])
			tmp.index <- findInterval(scale_score, tmp.old.cuts, rightmost.closed=TRUE)
			tmp.diff <- diff(tmp.new.cuts)/diff(tmp.old.cuts)
			round(tmp.new.cuts[tmp.index] + (scale_score - tmp.old.cuts[tmp.index]) * (diff(tmp.new.cuts)/diff(tmp.old.cuts))[tmp.index], digits=output.digits)
		} else {
			as.numeric(scale_score)
		}
	} else {
		if (!is.na(content_area) & !is.na(grade)) {
			sgp.projections.equated[['Linkages']][[paste(content_area, sgp.projections.equated[['Year']], sep=".")]][[paste("GRADE", grade, sep="_")]][[toupper(equating.method)]][['OLD_TO_NEW']][['interpolated_function']](scale_score)
		} else {
			as.numeric(scale_score)
		}
	}
} ## END piecewiseTransform Function
