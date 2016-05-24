`getAchievementLevel` <-
function(sgp_data,
	state=NULL,
	year=NULL,
	content_area=NULL,
	grade=NULL,
	achievement.level.name="ACHIEVEMENT_LEVEL",
	scale.score.name="SCALE_SCORE") {

	CONTENT_AREA <- YEAR <- GRADE <- STATE <- TMP_ACH_LEVEL <- SCALE_SCORE <- ACHIEVEMENT_LEVEL <- NULL

	###  Utility functions

	get.cutscore.label <- function(state, year, content_area) {
		tmp.cutscore.names <- names(SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]])
		tmp.cutscore.years <- sapply(strsplit(tmp.cutscore.names[grep(content_area, tmp.cutscore.names)], "[.]"), function(x) x[2])
		if (any(!is.na(tmp.cutscore.years))) {
			if (year %in% tmp.cutscore.years) {
				return(paste(content_area, year, sep="."))
			} else {
				if (year==sort(c(year, tmp.cutscore.years))[1]) {
					return(content_area)
				} else {
					return(paste(content_area, sort(tmp.cutscore.years)[which(year==sort(c(year, tmp.cutscore.years)))-1], sep="."))
				}
			}
		} else {
			return(content_area)
		}
	}

	get.achievement_level.label <- function(state, year) {
		tmp.achievement_level.names <- grep("Achievement_Levels", names(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]]), value=TRUE)
		tmp.achievement_level.years <- sapply(strsplit(tmp.achievement_level.names, "[.]"), function(x) x[2])
		if (any(!is.na(tmp.achievement_level.years))) {
			if (year %in% tmp.achievement_level.years) {
				return(paste("Achievement_Levels", year, sep="."))
			} else {
				if (year==sort(c(year, tmp.achievement_level.years))[1]) {
					return("Achievement_Levels")
				} else {
					return(paste("Achievement_Levels", sort(tmp.achievement_level.years)[which(year==sort(c(year, tmp.achievement_level.years)))-1], sep="."))
				}
			}
		} else {
			return("Achievement_Levels")
		}
	}

	getAchievementLevel_INTERNAL <- function(state, content_area, year, grade, scale_score) {
		if (is.null(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]])) {
			tmp.levels <- seq_along(SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]][!is.na(SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Proficient"]])])
			tmp.labels <- SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]][!is.na(SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Proficient"]])]
		} else {
			tmp.levels <- seq_along(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[get.achievement_level.label(state, year)]][["Labels"]][
				!is.na(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[get.achievement_level.label(state, year)]][["Proficient"]])])
			tmp.labels <- SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[get.achievement_level.label(state, year)]][["Labels"]][
				!is.na(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[get.achievement_level.label(state, year)]][["Proficient"]])]
		}
		as.character(factor(findInterval(scale_score, SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]][[get.cutscore.label(state, year, content_area)]][[paste("GRADE_", grade, sep="")]])+1,
			levels=tmp.levels, labels=tmp.labels))
	}

	if ("STATE" %in% names(sgp_data) & !is.null(SGP::SGPstateData[[state]][["Achievement"]][["Cutscore_Information"]])) {
		cutscore.states <- SGP::SGPstateData[[state]][["Achievement"]][["Cutscore_Information"]][["Cutscore_States"]]
		cutscore.subjects <- unique(sapply(names(SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]]), function(x) strsplit(x, "[.]")[[1]][1], USE.NAMES=FALSE))
		if (any(unique(sgp_data[['STATE']]) %in% cutscore.states)) {
			sgp_data[which(STATE %in% cutscore.states & CONTENT_AREA %in% cutscore.subjects), ACHIEVEMENT_LEVEL := paste("Level", findInterval(SCALE_SCORE,
				SGP::SGPstateData[[state]][["Achievement"]][["Cutscores"]][[paste(CONTENT_AREA[1], STATE[1], sep=".")]][[paste("GRADE", GRADE[1], sep="_")]])+1L),
				by=c("STATE", "CONTENT_AREA", "GRADE")]
		} else {
			sgp_data[,ACHIEVEMENT_LEVEL:=as.character(NA)]
		}
	} else {
		if (is.null(year)) year <- sort(unique(sgp_data[['YEAR']]))
		if (is.null(content_area)) content_area <- sort(unique(sgp_data[['CONTENT_AREA']][sgp_data[['YEAR']] %in% year]))
		if (is.null(grade)) grade <- sort(unique(sgp_data[['GRADE']][sgp_data[['YEAR']] %in% year & sgp_data[['CONTENT_AREA']] %in% content_area]))

		setkeyv(sgp_data, c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE"))
			sgp_data[sgp_data[CJ("VALID_CASE", content_area, year, grade), which=TRUE, nomatch=0], achievement.level.name :=
			sgp_data[CJ("VALID_CASE", content_area, year, grade), nomatch=0][, getAchievementLevel_INTERNAL(state, CONTENT_AREA, YEAR, GRADE, get(scale.score.name)),
				by=list(CONTENT_AREA, YEAR, GRADE)][["V1"]], with=FALSE]
	}
	return(sgp_data)
} ### END getAchievementLevel Function
