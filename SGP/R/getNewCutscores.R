`getNewCutscores` <- 
function(content_area,
	content_area_labels,
	year,
	grade,
	Cutscores) {

	### Utility function

	get.cutscore.year <- function(year, cutscore.years) {
		if (year %in% cutscore.years) return(year)
		if (year==sort(c(year, cutscore.years))[1]) return(as.character(NA))
		return(sort(cutscore.years)[which(year==sort(c(year, cutscore.years)))-1])
	}

	### Define variables

	cutscore.year <- get.cutscore.year(year, unique(Cutscores[[content_area]][['YEAR']]))

	if (is.na(content_area) | is.na(grade)) {
		return(NULL)
	} else {
		return(sort(Cutscores[[content_area]][list(content_area_labels, cutscore.year, grade)][["CUTSCORES"]]))
	}
} ### END getNewCutscores function
