`getTargetAchievementLevels` <- 
function(state,
	status.type="CATCH_UP_KEEP_UP") {

	tmp.list <- list()

	if (status.type=="CATCH_UP_KEEP_UP") {
		if (!is.null(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]])) {
			tmp.index <- grep("Achievement_Levels", names(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]]))
			levels.that.are.proficient <- 
				sort(unlist(sapply(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][tmp.index], function(x) x[['Labels']]))[
				unlist(sapply(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][tmp.index], function(x) x[['Proficient']]))=="Proficient"])
			levels.that.are.not.proficient <- 
				sort(unlist(sapply(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][tmp.index], function(x) x[['Labels']]))[
				unlist(sapply(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][tmp.index], function(x) x[['Proficient']]))=="Not Proficient"])
		} else {
			levels.that.are.proficient <- SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]][
				which(SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Proficient"]]=="Proficient")]
			levels.that.are.not.proficient <- SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]][
				which(SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Proficient"]]=="Not Proficient")]
		}
		return(list(NO=as.character(levels.that.are.not.proficient), YES=as.character(levels.that.are.proficient)))
	}

	if (status.type=="MOVE_UP_STAY_UP") {
		if (!is.null(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]])) {
			levels.that.are.advanced <- levels.that.are.not.advanced <- list()
			for (i in grep("Achievement_Levels", names(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]]), value=TRUE)) {
				levels.that.are.advanced[[i]] <- SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[i]][['Labels']][
					tail(which(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[i]][["Proficient"]]=="Proficient"), -1)]
				levels.that.are.not.advanced[[i]] <- SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[i]][['Labels']][
					c(which(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[i]][["Proficient"]]=="Not Proficient"),
					head(which(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[i]][["Proficient"]]=="Proficient"), 1))]
			}
		} else {
			levels.that.are.advanced <- SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]][
				tail(which(SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Proficient"]]=="Proficient"), -1)]
			levels.that.are.not.advanced <- SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]][
				c(which(SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Proficient"]]=="Not Proficient"),
				head(which(SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Proficient"]]=="Proficient"), 1))]
		}
		return(list(NO=as.character(unlist(levels.that.are.not.advanced)), YES=as.character(unlist(levels.that.are.advanced))))
	}
} ### END getTargetAchievementLevels
