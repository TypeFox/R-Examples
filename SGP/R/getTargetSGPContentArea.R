`getTargetSGPContentArea` <- 
function(grade,
	content_area,
	state,
	max.sgp.target.years.forward,
	my.sgp.target.content_area) {

	if (is.null(SGP::SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[content_area]])) {
		return(content_area)
	} else {
		tmp.content_areas.by.grade <- paste(SGP::SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[content_area]], 
						SGP::SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[content_area]], sep=".")
		tmp.index <- match(paste(content_area, grade, sep="."), tmp.content_areas.by.grade)
		if (is.na(tmp.index)) {
			message(paste("\tNOTE: '", content_area, "', GRADE '", grade, "' combination is not in current @Data. Will return ", content_area, " for '", my.sgp.target.content_area, "'.", sep=""))
			return(content_area)
		} else {
			return(SGP::SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[content_area]][
				min(tmp.index+max.sgp.target.years.forward, length(SGP::SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[content_area]]))])
		}
	}
} ### END getTargetSGPContentArea
