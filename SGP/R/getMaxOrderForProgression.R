`getMaxOrderForProgression` <-
function(year,
	content_area,
	state,
	sgp.projections.equated) {

	if (is.null(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[content_area]])) {
		return(SGP::SGPstateData[[state]][["SGP_Configuration"]][["max.order.for.projection"]]) ## Returns NULL if it doesn't exist
	} else {
		if (is.null(sgp.projections.equated)) {
			tmp <- as.numeric(tail(unlist(strsplit(as.character(year), "_")), 1)) -
				as.numeric(tail(unlist(strsplit(as.character(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[content_area]]), "_")), 1))
			if (tmp < 0) return(SGP::SGPstateData[[state]][["SGP_Configuration"]][["max.order.for.projection"]])
			if (tmp > 0) return(min(c(as.numeric(tmp), SGP::SGPstateData[[state]][["SGP_Configuration"]][["max.order.for.projection"]])))
			if (tmp==0 &&
				!identical(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][["Baseline_Projections_in_Transition_Year"]], TRUE)) {
					message(paste("\tNOTE: Based upon state scale changes in ", year, ". student growth projections are not possible. No student growth projections will be generated", sep="")) 
			}
		} else {
			return(SGP::SGPstateData[[state]][["SGP_Configuration"]][["max.order.for.projection"]]) ## Returns NULL if it doesn't exist
		}
	}
} ### getMaxOrderForProgression
