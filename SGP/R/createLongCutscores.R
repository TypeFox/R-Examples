`createLongCutscores` <-
function(state,
	content_area,
	add.GRADE_NUMERIC=FALSE,
	assessment.transition.type=NULL) {

	GRADE <- GRADE_NUMERIC <- CUTSCORES <- YEAR <- CUTLEVEL <- YEAR_LAG <- CONTENT_AREA <- SCALE_SCORE <- first.content_area <- last.content_area <- NULL

	### Create relevant variables

	tmp.cutscore.list <- list()
	content_area.argument <- content_area
	if (!is.null(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Content_Areas_Domains"]])) {
		first.content_area <- head(names(which(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Content_Areas_Domains"]]==content_area)),1)
		last.content_area <- tail(names(which(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Content_Areas_Domains"]]==content_area)),1)
		content_area <- unique(names(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Content_Areas_Domains"]])[
			SGP::SGPstateData[[state]][["Student_Report_Information"]][["Content_Areas_Domains"]]==content_area])
	}
	if (length(content_area)>1) multiple.content_areas <- TRUE else multiple.content_areas <- FALSE


	### Utility functions

	get.long.cutscores <- function(content_area, transformed.cutscores=NULL, subset.year=NULL) {

		cutscore.list <- list()
		for (i in grep(paste("^", content_area, "$", sep=""), sapply(strsplit(names(SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]]), '[.]'), '[', 1))) {
			cutscores.content_area <- unlist(strsplit(names(SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]])[i], '[.]'))[1]
			grades <- as.character(matrix(unlist(strsplit(names(SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]][[i]]), "_")), ncol=2, byrow=TRUE)[,2])
			cutscores.year <- as.character(unlist(strsplit(names(SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]])[i], "[.]"))[2])

			if (!is.null(transformed.cutscores) & content_area %in% names(transformed.cutscores)) {
				cutscores.iter <- seq(length(transformed.cutscores[[content_area]])-2)
				cutscores <- tail(transformed.cutscores[[content_area]], -1)
				loss <- transformed.cutscores[[content_area]][1]
				hoss <- tail(transformed.cutscores[[content_area]], 1)
			} else {
				cutscores.iter <- seq(length(SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]][[i]][[1]]))
				cutscores <- as.data.table(matrix(unlist(SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]][[i]]),
					ncol=length(SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]][[i]][[1]]), byrow=TRUE))
				if (names(SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]])[i] %in% names(SGP::SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]])) {
					loss.hoss.label <- names(SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]])[i]
				} else {
					loss.hoss.label <- cutscores.content_area
				}
				loss <- sapply(SGP::SGPstateData[[state]][['Achievement']][['Knots_Boundaries']][[loss.hoss.label]][
						grep("loss.hoss", names(SGP::SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]][[loss.hoss.label]]))], '[', 1)
				loss <- as.numeric(loss[sapply(strsplit(names(loss), "_"), '[', 2) %in% grades])
				hoss <- sapply(SGP::SGPstateData[[state]][['Achievement']][['Knots_Boundaries']][[loss.hoss.label]][
						grep("loss.hoss", names(SGP::SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]][[loss.hoss.label]]))], '[', 2)
				hoss <- as.numeric(hoss[sapply(strsplit(names(hoss), "_"), '[', 2) %in% grades])
			}

			for (j in cutscores.iter) {
				cutscore.list[[paste(i, j, sep="_")]] <- data.table(
					GRADE=grades,
					CONTENT_AREA=cutscores.content_area,
					CUTLEVEL=as.character(j),
					CUTSCORES=cutscores[[j]],
					YEAR=cutscores.year)
				cutscore.list[[paste(i, j, sep="_")]] <- subset(cutscore.list[[paste(i, j, sep="_")]],
					GRADE %in% SGP::SGPstateData[[state]][["Student_Report_Information"]][["Grades_Reported"]][[content_area]])
			}

			cutscore.list[[paste(i, "LOSS", sep="_")]] <- data.table(
				GRADE=grades,
				CONTENT_AREA=cutscores.content_area,
				CUTLEVEL="LOSS",
				CUTSCORES=loss,
				YEAR=cutscores.year)
			cutscore.list[[paste(i, "LOSS", sep="_")]] <- subset(cutscore.list[[paste(i, "LOSS", sep="_")]],
				GRADE %in% SGP::SGPstateData[[state]][["Student_Report_Information"]][["Grades_Reported"]][[content_area]])

			cutscore.list[[paste(i, "HOSS", sep="_")]] <- data.table(
				GRADE=grades,
				CONTENT_AREA=cutscores.content_area,
				CUTLEVEL="HOSS",
				CUTSCORES=hoss,
				YEAR=cutscores.year)
			cutscore.list[[paste(i, "HOSS", sep="_")]] <- subset(cutscore.list[[paste(i, "HOSS", sep="_")]],
				GRADE %in% SGP::SGPstateData[[state]][["Student_Report_Information"]][["Grades_Reported"]][[content_area]])
		}

		### Add GRADE_LOWER/GRADE_UPPER
		long.cutscores <- rbindlist(cutscore.list)
		if (all(long.cutscores$GRADE == "EOCT")) {
			if (content_area == last.content_area) {
				extension.cutscores <-
					# data.table(CONTENT_AREA="PLACEHOLDER", GRADE="GRADE_UPPER", subset(long.cutscores, !is.na(as.numeric(CUTLEVEL)))[,list(CUTSCORES=CUTSCORES+CUTSCORES*0.15), by=list(YEAR, CUTLEVEL)])
					# data.table(CONTENT_AREA="PLACEHOLDER", GRADE="GRADE_UPPER", long.cutscores[,list(CUTSCORES=CUTSCORES+CUTSCORES*0.15), by=list(YEAR, CUTLEVEL)])
					data.table(CONTENT_AREA="PLACEHOLDER", GRADE="GRADE_UPPER", long.cutscores[,list(CUTSCORES=extendrange(rnorm(rep(CUTSCORES, 1000), CUTSCORES, CUTSCORES/500), f=0.075)[2]), by=list(YEAR, CUTLEVEL)])
			} else {
				extension.cutscores <-
					data.table(CONTENT_AREA="PLACEHOLDER", GRADE="GRADE_LOWER", long.cutscores[,list(CUTSCORES=extendrange(rnorm(rep(CUTSCORES, 1000), CUTSCORES, CUTSCORES/500), f=0.075)[1]), by=list(YEAR, CUTLEVEL)])

			}
		} else {
			extension.cutscores <-
				data.table(CONTENT_AREA="PLACEHOLDER", GRADE=c("GRADE_LOWER", "GRADE_UPPER"), long.cutscores[,list(CUTSCORES=extendrange(CUTSCORES, f=0.15)), by=list(YEAR, CUTLEVEL)])
		}
		long.cutscores <- rbindlist(list(long.cutscores, setcolorder(extension.cutscores, names(cutscore.list[[1]]))))
		setkeyv(long.cutscores, c("GRADE", "CONTENT_AREA"))

		if (length(sort(long.cutscores$YEAR)) > 0 & !is.null(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Earliest_Year_Reported"]][[content_area]])) {
			long.cutscores <- subset(long.cutscores, as.numeric(unlist(sapply(strsplit(as.character(long.cutscores$YEAR), "_"), function(x) x[1]))) >=
				as.numeric(sapply(strsplit(as.character(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Earliest_Year_Reported"]][[content_area]]), "_"), function(x) x[1])))
		}

		if (!is.null(subset.year)) if (identical(subset.year, NA)) long.cutscores <- long.cutscores[is.na(YEAR)] else long.cutscores <- long.cutscores[YEAR==subset.year]

		if (multiple.content_areas) {
			if (content_area %in% c(first.content_area, last.content_area)) {
				if (content_area == first.content_area) long.cutscores <- long.cutscores[GRADE!="GRADE_UPPER"]
				if (content_area == last.content_area) long.cutscores <- long.cutscores[GRADE!="GRADE_LOWER"]
			} else long.cutscores <- long.cutscores[!GRADE %in% c("GRADE_LOWER", "GRADE_UPPER")]
		}

		return(data.table(long.cutscores, key=c("GRADE", "CONTENT_AREA")))
	} ### END get.long.cutscores


	###############################################################################################################################
	### Create long cutscores based upon whether an assessment transition has occurred
	###############################################################################################################################`

	if (is.null(assessment.transition.type)) {
		if (any(content_area %in% names(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Transformed_Achievement_Level_Cutscores"]])) &
			!all(content_area %in% names(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Transformed_Achievement_Level_Cutscores"]]))) {
				stop("Not all content areas have Transformed Achievement Level Cutscores available in SGP::SGPstateData.
					Please augment the SGP::SGPstateData set with your data or contact the SGP package maintainer to have your data added to the SGP package.")
		}
		tmp.transformed.cutscores <- SGP::SGPstateData[[state]][["Student_Report_Information"]][["Transformed_Achievement_Level_Cutscores"]]
	} else tmp.transformed.cutscores <- NULL

	if (multiple.content_areas) {
		long.cutscores <- data.table(rbindlist(lapply(content_area, function(x) get.long.cutscores(x, tmp.transformed.cutscores))), key=c("GRADE", "CONTENT_AREA"))
	} else long.cutscores <- get.long.cutscores(content_area, tmp.transformed.cutscores)

	### Add GRADE_NUMERIC

	if (add.GRADE_NUMERIC) {
		if (!is.null(SGP::SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[content_area.argument]])) {
			grades.content_areas.reported.in.state <- data.table(
					GRADE=as.character(SGP::SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[content_area.argument]]),
					YEAR_LAG=c(1, SGP::SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[content_area.argument]]),
					CONTENT_AREA=SGP::SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[content_area.argument]])
		} else {
			grades.content_areas.reported.in.state <- data.table(
					GRADE=as.numeric(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Grades_Reported"]][[content_area.argument]]),
					YEAR_LAG=c(1, diff(as.numeric(SGP::SGPstateData[[state]][["Student_Report_Information"]][["Grades_Reported"]][[content_area.argument]]))),
					CONTENT_AREA=content_area.argument)
		}
		grades.content_areas.reported.in.state[,GRADE_NUMERIC:=
			(as.numeric(grades.content_areas.reported.in.state$GRADE[2])-1)+c(0, cumsum(tail(grades.content_areas.reported.in.state$YEAR_LAG, -1)))]
		grades.content_areas.reported.in.state[,GRADE:=as.character(grades.content_areas.reported.in.state$GRADE)]
		setkeyv(grades.content_areas.reported.in.state, c("GRADE", "CONTENT_AREA"))
		long.cutscores <- grades.content_areas.reported.in.state[long.cutscores]
		long.cutscores[,YEAR_LAG:=NULL]
	}

	return(unique(data.table(long.cutscores, key=c("CONTENT_AREA", "YEAR", "GRADE", "CUTLEVEL"))))
} ## END createLongCutscores Function
