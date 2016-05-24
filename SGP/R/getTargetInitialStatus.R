`getTargetInitialStatus` <- 
function(achievement_level,
	state,
	state.iter=NULL,
	status.type="CATCH_UP_KEEP_UP") {

		if (!is.null(SGP::SGPstateData[[state]][['Achievement']][['Cutscore_Information']])) {
			tmp.state.level <- which(sapply(lapply(SGP::SGPstateData[[state]][["Achievement"]][["Cutscore_Information"]][['State_Levels']], '[[', 1), function(x) state.iter %in% x))
			levels.that.are.proficient <- SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]][
				which(SGP::SGPstateData[[state]][["Achievement"]][["Cutscore_Information"]][["State_Levels"]][[tmp.state.level]][['Levels']]=="Proficient")]
			levels.that.are.not.proficient <- SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]][
				which(SGP::SGPstateData[[state]][["Achievement"]][["Cutscore_Information"]][["State_Levels"]][[tmp.state.level]][['Levels']]=="Not Proficient")]
			levels.that.are.advanced <- SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]][
				tail(which(SGP::SGPstateData[[state]][["Achievement"]][["Cutscore_Information"]][["State_Levels"]][[tmp.state.level]][["Levels"]]=="Proficient"), -1)]
			levels.that.are.not.advanced <- SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]][
				c(which(SGP::SGPstateData[[state]][["Achievement"]][["Cutscore_Information"]][["State_Levels"]][[tmp.state.level]][['Levels']]=="Not Proficient"),
				head(which(SGP::SGPstateData[[state]][["Achievement"]][["Cutscore_Information"]][["State_Levels"]][[tmp.state.level]][['Levels']]=="Proficient"), 1))]
		} else {
			if (!is.null(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]])) {
				levels.that.are.advanced <- levels.that.are.not.advanced <- list()
				tmp.index <- grep("Achievement_Levels", names(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]]))
				levels.that.are.proficient <- sort(unlist(sapply(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][tmp.index], function(x) x[['Labels']]))[
					unlist(sapply(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][tmp.index], function(x) x[['Proficient']]))=="Proficient"])
				levels.that.are.not.proficient <- sort(unlist(sapply(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][tmp.index], function(x) x[['Labels']]))[
					unlist(sapply(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][tmp.index], function(x) x[['Proficient']]))=="Not Proficient"])
				for (i in grep("Achievement_Levels", names(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]]), value=TRUE)) {
					levels.that.are.advanced[[i]] <- SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[i]][['Labels']][
						tail(which(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[i]][["Proficient"]]=="Proficient"), -1)]
					levels.that.are.not.advanced[[i]] <- SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[i]][['Labels']][
						c(which(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[i]][["Proficient"]]=="Not Proficient"),
						head(which(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[i]][["Proficient"]]=="Proficient"), 1))]
				}
				levels.that.are.advanced <- unlist(levels.that.are.advanced)
				levels.that.are.not.advanced <- unlist(levels.that.are.not.advanced)
			} else {
				levels.that.are.proficient <- SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]][
					which(SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Proficient"]]=="Proficient")]
				levels.that.are.not.proficient <- SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]][
					which(SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Proficient"]]=="Not Proficient")]
				levels.that.are.advanced <- SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]][
					tail(which(SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Proficient"]]=="Proficient"), -1)]
				levels.that.are.not.advanced <- SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Labels"]][
					c(which(SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Proficient"]]=="Not Proficient"),
					head(which(SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Proficient"]]=="Proficient"), 1))]
			}
		}

		if (status.type=="CATCH_UP_KEEP_UP") {
			achievement_level[achievement_level %in% levels.that.are.proficient] <- 2
			achievement_level[achievement_level %in% levels.that.are.not.proficient] <- 1
			return(factor(achievement_level, levels=1:2, labels=c("Catching Up", "Keeping Up")))
		}

		if (status.type=="MOVE_UP_STAY_UP") {
			achievement_level[achievement_level %in% levels.that.are.not.proficient] <- NA
			achievement_level[achievement_level %in% levels.that.are.advanced] <- 2
			achievement_level[achievement_level %in% levels.that.are.not.advanced] <- 1
			return(factor(factor(achievement_level, levels=1:2, labels=c("Moving Up", "Staying Up"))))
		}
} ### END getTargetInitialStatus
