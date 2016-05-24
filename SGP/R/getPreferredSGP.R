`getPreferredSGP` <- 
function(tmp.data,
	state,
	type="COHORT_REFERENCED") {

	YEAR <- SGP_NORM_GROUP <- VALID_CASE <- CONTENT_AREA <- ID <- PREFERENCE <- NULL

	if (type=="TARGET") {
		tmp.sgp.norm.group.variables <- c("SGP_PROJECTION_GROUP", "PREFERENCE")
	}
	if (type=="BASELINE") {
		tmp.sgp.norm.group.variables <- c("YEAR", "SGP_NORM_GROUP_BASELINE", "PREFERENCE")
		tmp.message <- "\tNOTE: Multiple Baseline SGPs exist for individual students. Unique Baseline SGPs will be created using SGP Norm Group Preference Table for "
	} 
	if (type=="COHORT_REFERENCED") {
		tmp.sgp.norm.group.variables <- c("YEAR", "SGP_NORM_GROUP", "PREFERENCE")
		tmp.message <- "\tNOTE: Multiple SGPs exist for individual students. Unique SGPs will be created using SGP Norm Group Preference Table for "
	}

	if (type=="TARGET") {
		if (!is.null(SGP::SGPstateData[[state]][['SGP_Progression_Preference']])) {
			setkeyv(SGP::SGPstateData[[state]][['SGP_Progression_Preference']], tmp.sgp.norm.group.variables)
			setkeyv(tmp.data, tmp.sgp.norm.group.variables[1])
		} else {
			stop("\tNOTE: Multiple Projections exist for individual students. Please examine results in @SGP[['SGProjections']].")
		}
		tmp.data <- data.table(SGP::SGPstateData[[state]][['SGP_Progression_Preference']][,tmp.sgp.norm.group.variables,with=FALSE][tmp.data], 
			key=c(getKey(tmp.data), "PREFERENCE"))
	} else {
		if (!is.null(SGP::SGPstateData[[state]][['SGP_Norm_Group_Preference']])) {
			message(paste(tmp.message, state, ".", sep=""))
			setkeyv(SGP::SGPstateData[[state]][['SGP_Norm_Group_Preference']], tmp.sgp.norm.group.variables[1:2])
			setkeyv(tmp.data, tmp.sgp.norm.group.variables[1:2])
		} else {
			stop("\tNOTE: Multiple SGPs exist for individual students. Please examine results in @SGP[['SGPercentiles']].")
		}
		tmp.data <- data.table(SGP::SGPstateData[[state]][['SGP_Norm_Group_Preference']][,tmp.sgp.norm.group.variables,with=FALSE][tmp.data], 
			key=c(getKey(tmp.data), "PREFERENCE"))
	}
	setkeyv(tmp.data, getKey(tmp.data))
	tmp.data <- tmp.data[!duplicated(tmp.data)][,PREFERENCE:=NULL]
	setkeyv(tmp.data, getKey(tmp.data))
	return(tmp.data)
} ### END getPreferredSGP
