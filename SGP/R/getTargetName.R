`getTargetName` <- 
function(state,
	target.type="sgp.projections.lagged",
	target.level="CATCH_UP_KEEP_UP",
	target.years=3,
	target.label="SGP_TARGET",
	projection.unit.label="YEAR",
	projection_group.iter=NULL) {

	if (target.level=="CATCH_UP_KEEP_UP") target.level <- NULL
	if (target.type %in% c("sgp.projections", "sgp.projections.baseline")) target.period <- "CURRENT" else target.period <- NULL
	if (target.type %in% c("sgp.projections.baseline", "sgp.projections.lagged.baseline")) sgp.type <- "BASELINE" else sgp.type <- NULL
	if (!is.null(projection_group.iter) && projection_group.iter %in% names(SGP::SGPstateData[[state]][['SGP_Configuration']][['grade.projection.sequence']])) {
		tmp.target.years <- SGP::SGPstateData[[state]][['SGP_Configuration']][['max.forward.projection.sequence']][[projection_group.iter]]
		if (!is.null(tmp.target.years)) target.years <- tmp.target.years
	}

	return(paste(c(target.label, sgp.type, target.level, target.years, projection.unit.label, target.period), collapse="_"))
} ### END getTargetName
