`transformScaleScore` <-
function(tmp.data,
	state,
	content_areas,
	linkages,
	slot.data,
	equating.method="equipercentile",
	assessment.transition.type.cutscores=NULL) {

	TRANSFORMED_SCALE_SCORE <- SCALE_SCORE <- TEMP_SCALE_SCORE <- SCALE_SCORE_EQUATED <- CONTENT_AREA <- CONTENT_AREA_LABELS <- YEAR <- GRADE <- GRADE_NUMERIC <- ID <- NULL
	CUTSCORES <- CUTSCORES_ORIGINAL <- GRADE_FOR_CUTSCORES <- NULL

	### Create relevant variables

	Cutscores <- list()


	### Utility functions

	get.min.max.grade <- function(Cutscores) {

		if ("GRADE_NUMERIC" %in% names(Cutscores)) {
			tmp.grades.numeric <- range(sort(type.convert(subset(Cutscores, !GRADE %in% c("GRADE_LOWER", "GRADE_UPPER"))[['GRADE_NUMERIC']])))
			tmp.grades <- sort(subset(Cutscores, GRADE_NUMERIC %in% tmp.grades.numeric)[['GRADE']])
		} else {
			tmp.grades <- sort(type.convert(subset(Cutscores, !GRADE %in% c("GRADE_LOWER", "GRADE_UPPER"))[['GRADE']]))
		}
		return(c(tmp.grades[1], rev(tmp.grades)[1]))
	}


	### Return Data and Cutscores based upon whether scale score transition

	if (!is.null(linkages)) {

		### Define variables

		year.for.equate <- tail(sort(sapply(strsplit(names(linkages), "[.]"), '[', 2)), 1)
		assessment.transition.type <- c(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][['Vertical_Scale']],
			SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[paste('Vertical_Scale', year.for.equate, sep=".")]])
		if (is.null(assessment.transition.type.cutscores)) {
			assessment.transition.type.cutscores <- assessment.transition.type
		}

		for (i in content_areas) {
			Cutscores[[i]] <- createLongCutscores(state=state, content_area=i, assessment.transition.type=assessment.transition.type.cutscores)
			Cutscores[[i]][, CUTSCORES_ORIGINAL:=CUTSCORES]
		}


		### Transform Cutscores

		for (content_area.iter in content_areas) {
			for (grade.iter in unique(Cutscores[[content_area.iter]][['GRADE']])) {
				Cutscores[[content_area.iter]][CONTENT_AREA==content_area.iter & GRADE==grade.iter & (is.na(YEAR) | YEAR < year.for.equate),
					CUTSCORES:=linkages[[paste(content_area.iter, year.for.equate, sep=".")]][[paste("GRADE", grade.iter, sep="_")]][[toupper(equating.method)]][['OLD_TO_NEW']][["interpolated_function"]](CUTSCORES)]
				tmp.min.max <- get.min.max.grade(Cutscores[[content_area.iter]])
				if (grade.iter=="GRADE_UPPER") {
					Cutscores[[content_area.iter]][CONTENT_AREA=="PLACEHOLDER" & GRADE=="GRADE_UPPER" & (is.na(YEAR) | YEAR < year.for.equate),
						CUTSCORES:=linkages[[paste(content_area.iter, year.for.equate, sep=".")]][[paste("GRADE", tmp.min.max[2], sep="_")]][[toupper(equating.method)]][['OLD_TO_NEW']][["interpolated_function"]](CUTSCORES)]
				}
				if (grade.iter=="GRADE_LOWER") {
					Cutscores[[content_area.iter]][CONTENT_AREA=="PLACEHOLDER" & GRADE=="GRADE_LOWER" & (is.na(YEAR) | YEAR < year.for.equate),
						CUTSCORES:=linkages[[paste(content_area.iter, year.for.equate, sep=".")]][[paste("GRADE", tmp.min.max[1], sep="_")]][[toupper(equating.method)]][['OLD_TO_NEW']][["interpolated_function"]](CUTSCORES)]
				}
			}
		}

		#############################################################
		### Vertical to Vertical scale transition
		#############################################################

		if (identical(toupper(assessment.transition.type), c("YES", "YES"))) {

			### Create TRANSFORMED_SCALE_SCORE

			tmp.data[, TRANSFORMED_SCALE_SCORE:=SCALE_SCORE_EQUATED]
		}


		#######################################################
		### Non-Vertical to Vertical scale transition
		#######################################################

		if (identical(toupper(assessment.transition.type), c("NO", "YES"))) {

			### Create TRANSFORMED_SCALE_SCORE

			tmp.data[!is.na(CONTENT_AREA_LABELS) & YEAR <= year.for.equate, GRADE_FOR_CUTSCORES:=tail(mixedsort(sort(GRADE)), 1), by=list(CONTENT_AREA_LABELS, ID)]
			tmp.data[!is.na(CONTENT_AREA_LABELS) & YEAR < year.for.equate,
				TRANSFORMED_SCALE_SCORE:=piecewiseTransform(
					SCALE_SCORE,
					state,
					CONTENT_AREA_LABELS,
					YEAR,
					GRADE,
					new.cutscores=sort(Cutscores[[CONTENT_AREA_LABELS[1]]][list(CONTENT_AREA_LABELS[1], rev(sort(unique(Cutscores[[CONTENT_AREA_LABELS[1]]]$YEAR), na.last=FALSE))[2], GRADE_FOR_CUTSCORES[1])][['CUTSCORES']])),
						by=list(CONTENT_AREA_LABELS, YEAR, GRADE, GRADE_FOR_CUTSCORES)]
			tmp.data[!is.na(CONTENT_AREA_LABELS) & YEAR >= year.for.equate, TRANSFORMED_SCALE_SCORE:=SCALE_SCORE_EQUATED]
		}


		#######################################################
		### Vertical to Non-Vertical scale transition
		#######################################################

		if (identical(toupper(assessment.transition.type), c("YES", "NO"))) {

			### Create TRANSFORMED_SCALE_SCORE

			tmp.data[!is.na(CONTENT_AREA_LABELS) & YEAR >= year.for.equate, GRADE_FOR_CUTSCORES:=head(mixedsort(sort(GRADE)), 1), by=list(CONTENT_AREA_LABELS, ID)]
			tmp.data[!is.na(CONTENT_AREA_LABELS) & YEAR >= year.for.equate,
				TRANSFORMED_SCALE_SCORE:=piecewiseTransform(
					SCALE_SCORE,
					state,
					CONTENT_AREA_LABELS,
					YEAR,
					GRADE,
					new.cutscores=sort(Cutscores[[CONTENT_AREA_LABELS[1]]][list(CONTENT_AREA_LABELS[1], year.for.equate, GRADE_FOR_CUTSCORES[1])][['CUTSCORES']])),
						by=list(CONTENT_AREA_LABELS, YEAR, GRADE)]
			tmp.data[!is.na(CONTENT_AREA_LABELS) & YEAR < year.for.equate, TRANSFORMED_SCALE_SCORE:=SCALE_SCORE_EQUATED]
		}


		###############################################################
		### Non-Vertical to Non-Vertical scale transition
		###############################################################

		if (identical(toupper(assessment.transition.type), c("NO", "NO"))) {

			### Create TRANSFORMED_SCALE_SCORE

			tmp.data[!is.na(CONTENT_AREA_LABELS) & YEAR <= year.for.equate, GRADE_FOR_CUTSCORES:=tail(mixedsort(sort(GRADE)), 1), by=list(CONTENT_AREA_LABELS, ID)]
			tmp.data[!is.na(CONTENT_AREA_LABELS) & YEAR < year.for.equate,
				TRANSFORMED_SCALE_SCORE:=piecewiseTransform(
					SCALE_SCORE,
					state,
					CONTENT_AREA_LABELS,
					YEAR,
					GRADE,
					new.cutscores=sort(Cutscores[[CONTENT_AREA_LABELS[1]]][list(CONTENT_AREA_LABELS[1], rev(sort(unique(Cutscores[[CONTENT_AREA_LABELS[1]]]$YEAR), na.last=FALSE))[2], GRADE_FOR_CUTSCORES[1])][['CUTSCORES']])),
						by=list(CONTENT_AREA_LABELS, YEAR, GRADE)]

			tmp.data[!is.na(CONTENT_AREA_LABELS) & YEAR >= year.for.equate, GRADE_FOR_CUTSCORES:=head(mixedsort(sort(GRADE)), 1), by=list(CONTENT_AREA_LABELS, ID)]
			tmp.data[!is.na(CONTENT_AREA_LABELS) & YEAR >= year.for.equate,
				TRANSFORMED_SCALE_SCORE:=piecewiseTransform(
					SCALE_SCORE,
					state,
					CONTENT_AREA_LABELS,
					YEAR,
					GRADE,
					new.cutscores=sort(Cutscores[[CONTENT_AREA_LABELS[1]]][list(CONTENT_AREA_LABELS[1], year.for.equate, GRADE_FOR_CUTSCORES[1])][['CUTSCORES']])),
						by=list(CONTENT_AREA_LABELS, YEAR, GRADE)]
		}

		### Return data

		return(list(Data=tmp.data, Cutscores=Cutscores, sgp.projections.equated=list(Year=year.for.equate, Linkages=linkages, Assessment_Transition_Type=assessment.transition.type)))
	} else {
		for (i in content_areas) {
			Cutscores[[i]] <- createLongCutscores(state, i)
			Cutscores[[i]][, CUTSCORES_ORIGINAL:=CUTSCORES]
		}
		tmp.data[, TRANSFORMED_SCALE_SCORE:=piecewiseTransform(SCALE_SCORE, state, CONTENT_AREA_LABELS, as.character(YEAR), as.character(GRADE)), by=list(CONTENT_AREA_LABELS, YEAR, GRADE)]
		return(list(Data=tmp.data, Cutscores=Cutscores, sgp.projections.equated=NULL))
	}
} ### END transformScaleScore function
