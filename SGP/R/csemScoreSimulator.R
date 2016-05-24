csemScoreSimulator <- 
function(
	scale_scores, 
	grade, 
	content_area, 
	year, 
	state, 
	variable=NULL, 
	distribution=NULL, 
	round.digits=NULL) {

	GRADE <- CONTENT_AREA <- YEAR <- SIM <- NULL

	### Define relevant variables

	if (is.null(round.digits)) round.digits <- 2
	if (is.null(distribution)) distribution <- "Normal"
	if (!is.null(state)) {
		min.max <- SGP::SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]][[content_area]][[paste("loss.hoss_", grade, sep="")]]
	} else {
		min.max <- range(scale_scores, na.rm=TRUE)
	}
	Interpolation_Function <- function(scale_score, variance) return(splinefun(scale_score, variance, method="natural"))

	### Create scale score dependent CSEMs

	if (!is.null(state)) {
		if ("YEAR" %in% names(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["CSEM"]])) {
			Interpolation_Data <- subset(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["CSEM"]], GRADE==grade & CONTENT_AREA==content_area & YEAR==year)
		} else {
			Interpolation_Data <- subset(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["CSEM"]], GRADE==grade & CONTENT_AREA==content_area)
		}
		tmp.omega <- Interpolation_Function(Interpolation_Data[['SCALE_SCORE']], Interpolation_Data[['SCALE_SCORE_CSEM']])(scale_scores)
	}
	if (!is.null(variable)) tmp.omega <- variable

	if (distribution=="Skew-Normal") {
		tmp.scores <- data.table(SIM=round(rsn(length(scale_scores), xi=scale_scores, omega=tmp.omega,
			alpha=tan((pi/2)*((min.max[1]+min.max[2]) - 2*scale_scores)/(min.max[2]-min.max[1]))), digits=round.digits))
	} else {
		tmp.scores <- data.table(SIM=round(rnorm(length(scale_scores), scale_scores, tmp.omega), digits=round.digits))
	}
	tmp.scores[SIM < min.max[1], SIM:=min.max[1]]
	tmp.scores[SIM > min.max[2], SIM:=min.max[2]]
	return(tmp.scores[['SIM']])
} ### END csemScoreSimulator
