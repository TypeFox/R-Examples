`getTargetSGP` <- 
function(sgp_object,
	content_areas, 
	state, 
	years, 
	target.type, 
	target.level,
	max.sgp.target.years.forward=3,
	subset.ids=NULL,
	return.lagged.status=TRUE) {

	VALID_CASE <- ID <- CONTENT_AREA <- YEAR <- FIRST_OBSERVATION <- LAST_OBSERVATION <- STATE <- SGP_PROJECTION_GROUP <- NULL

	### Utility functions

	getTargetSGP_INTERNAL <- function(tmp_object_1, state, state.iter, projection_group.iter, target.type, target.level, year_within) {

		if (year_within) {
			###  Assumes that any "canonical progression" will use the LAST_OBSERVATION for all (or at least the most recent) prior(s) in straight progressions
			if (target.type %in% c("sgp.projections", "sgp.projections.baseline")) {
				tmp_object_1[, LAST_OBSERVATION := 1L]; year.within.key <- "LAST_OBSERVATION"
			}
			###  lagged progressions would still be based on the FIRST_OBSERVATION score (used to produce SGP)
			if (target.type %in% c("sgp.projections.lagged", "sgp.projections.lagged.baseline")) {
				tmp_object_1[, FIRST_OBSERVATION := 1L]; year.within.key <- "FIRST_OBSERVATION"
			}
			setkeyv(tmp_object_1, c("VALID_CASE", "CONTENT_AREA", "YEAR", "ID", year.within.key))
			setkeyv(sgp_object@Data, c("VALID_CASE", "CONTENT_AREA", "YEAR", "ID", year.within.key))
			tmp_object_1 <- data.table(sgp_object@Data[,c(key(tmp_object_1), "YEAR_WITHIN"), with=FALSE], key=key(tmp_object_1))[tmp_object_1]
			jExp_Key <- c('VALID_CASE', 'CONTENT_AREA', 'YEAR', 'ID', 'YEAR_WITHIN')
		} else {
			setkeyv(tmp_object_1, c("VALID_CASE", "CONTENT_AREA", "YEAR", "ID"))
			jExp_Key <- c('VALID_CASE', 'CONTENT_AREA', 'YEAR', 'ID')
		}

		if (target.type %in% c("sgp.projections", "sgp.projections.baseline")) {
			if (year_within) {
				tmp_object_1 <- data.table(sgp_object@Data[,c(key(tmp_object_1), "ACHIEVEMENT_LEVEL"), with=FALSE], key=key(tmp_object_1))[tmp_object_1]
				setkeyv(sgp_object@Data, getKey(sgp_object))
			} else 	tmp_object_1 <- data.table(sgp_object@Data[,c(key(tmp_object_1), "ACHIEVEMENT_LEVEL"), with=FALSE], key=key(tmp_object_1))[tmp_object_1]
		}

		tmp_object_1[, paste(target.level, "STATUS_INITIAL", sep="_") := 
			getTargetInitialStatus(tmp_object_1[[grep("ACHIEVEMENT", names(tmp_object_1), value=TRUE)]], state, state.iter, target.level), with=FALSE]
		tmp_object_1 <- tmp_object_1[!is.na(get(paste(target.level, "STATUS_INITIAL", sep="_")))]
	
		## Find min/max of targets based upon CATCH_UP_KEEP_UP_STATUS_INITIAL status

		if (dim(tmp_object_1)[1] > 0) {
			num.years.available <- length(grep("LEVEL_[123456789]", names(tmp_object_1)))
			if (projection_group.iter %in% names(SGP::SGPstateData[[state]][['SGP_Configuration']][['grade.projection.sequence']])) {
				num.years.to.get <- min(SGP::SGPstateData[[state]][['SGP_Configuration']][['max.forward.projection.sequence']][[projection_group.iter]], num.years.available)
				if (!is.null(SGP::SGPstateData[[state]][['SGP_Configuration']][['max.forward.projection.sequence']][[projection_group.iter]])) {
					num.years.to.get.label <- SGP::SGPstateData[[state]][['SGP_Configuration']][['max.forward.projection.sequence']][[projection_group.iter]]
				} else {
					num.years.to.get.label <- max.sgp.target.years.forward
				}
			} else {
				num.years.to.get <- min(max.sgp.target.years.forward, num.years.available)
				num.years.to.get.label <- max.sgp.target.years.forward
			}
			if (target.type %in% c("sgp.projections.lagged", "sgp.projections.lagged.baseline")) num.years.to.get <- num.years.to.get+1
			
			tmp.level.variables <- 
				paste(grep(paste(sgp.projections.projection.unit.label, "_[", paste(seq(num.years.to.get), collapse=""), "]", sep=""), names(tmp_object_1), value=TRUE), collapse=", ")
	
			jExpression <- parse(text=paste("{catch_keep_move_functions[[unclass(", target.level, "_STATUS_INITIAL)]](", tmp.level.variables, ", na.rm=TRUE)}", sep=""))
			tmp_object_2 <- tmp_object_1[, eval(jExpression), keyby = jExp_Key]

			if (target.type %in% c("sgp.projections.baseline", "sgp.projections.lagged.baseline")) baseline.label <- "_BASELINE" else baseline.label <- NULL
			if (target.type %in% c("sgp.projections", "sgp.projections.baseline")) projection.label <- "_CURRENT" else projection.label <- NULL
			if (target.level=="MOVE_UP_STAY_UP") target.level.label <- "_MOVE_UP_STAY_UP" else target.level.label <- NULL
	
			setnames(tmp_object_2, "V1", 
				paste("SGP_TARGET", baseline.label, target.level.label, "_",  num.years.to.get.label, "_", sgp.projections.projection.unit.label, projection.label, sep=""))

			if (target.type %in% c("sgp.projections.lagged", "sgp.projections.lagged.baseline") && return.lagged.status) {
				tmp_object_2[,c("ACHIEVEMENT_LEVEL_PRIOR", grep("STATUS_INITIAL", names(tmp_object_1), value=TRUE)) := 
					list(tmp_object_1[["ACHIEVEMENT_LEVEL_PRIOR"]], tmp_object_1[[grep("STATUS_INITIAL", names(tmp_object_1), value=TRUE)]]), with=FALSE]
			}
			return(tmp_object_2[,SGP_PROJECTION_GROUP:=projection_group.iter])
		} else {
			return(NULL)
		}
	} ### getTargetSGP_INTERNAL

	
	### Define variables

	tmp.sgpTarget.list <- list()

	catch_keep_move_functions <- c(min, max)

	if (!is.null(SGP::SGPstateData[[state]][["SGP_Configuration"]][["sgp.projections.projection.unit.label"]])) {
		sgp.projections.projection.unit.label <- SGP::SGPstateData[[state]][["SGP_Configuration"]][["sgp.projections.projection.unit.label"]]
	} else {
		sgp.projections.projection.unit.label <- "YEAR"
	}

	### Loop over different states (usually just 1 state)

	tmp.names <- getPercentileTableNames(sgp_object, content_areas, state, years, target.type)
	if (length(tmp.names)==0) return(NULL)
	tmp.list <- list()

	if ("STATE" %in% names(sgp_object@Data)) {
		tmp.unique.states <- sort(unique(unlist(sapply(tmp.names, function(x) unique(sgp_object@SGP[['SGProjections']][[x]][['STATE']])))))
	} else {
		tmp.unique.states <- state
	}

	for (state.iter in tmp.unique.states) {
		level.to.get <- getTargetSGPLevel(state, state.iter, target.level)

		### Calculate Targets

		for (i in tmp.names) {
			cols.to.get.names <- names(sgp_object@SGP[["SGProjections"]][[i]])[
				c(grep(paste("LEVEL_", level.to.get, sep=""), names(sgp_object@SGP[["SGProjections"]][[i]])), grep("SGP_PROJECTION_GROUP", names(sgp_object@SGP[["SGProjections"]][[i]])))]
			if (target.type %in% c("sgp.projections.lagged", "sgp.projections.lagged.baseline")) cols.to.get.names <- c("ACHIEVEMENT_LEVEL_PRIOR", cols.to.get.names)
			if ("STATE" %in% names(sgp_object@Data)) cols.to.get.names <- c("STATE", cols.to.get.names)
			cols.to.get <- match(cols.to.get.names, names(sgp_object@SGP[["SGProjections"]][[i]]))

			if ("STATE" %in% names(sgp_object@Data)) {
				tmp.list[[i]] <- data.table(
					CONTENT_AREA=unlist(strsplit(i, "[.]"))[1],
					YEAR=getTableNameYear(i),
					sgp_object@SGP[["SGProjections"]][[i]][,c(1,cols.to.get), with=FALSE])[STATE==state.iter]
			} else {
				tmp.list[[i]] <- data.table(
					CONTENT_AREA=unlist(strsplit(i, "[.]"))[1],
					YEAR=getTableNameYear(i),
					sgp_object@SGP[["SGProjections"]][[i]][,c(1,cols.to.get), with=FALSE])
			}
		}

		if (!is.null(subset.ids)) {
			tmp_object_1 <- data.table(rbindlist(tmp.list, fill=TRUE), VALID_CASE="VALID_CASE", key="ID")[subset.ids, nomatch=0]
		} else {
			tmp_object_1 <- data.table(rbindlist(tmp.list, fill=TRUE), VALID_CASE="VALID_CASE")
		}

		if (!"SGP_PROJECTION_GROUP" %in% names(tmp_object_1)) tmp_object_1[,SGP_PROJECTION_GROUP:=CONTENT_AREA]

		for (projection_group.iter in unique(tmp_object_1[['SGP_PROJECTION_GROUP']])) {
			tmp.sgpTarget.list[[paste(state.iter, projection_group.iter, sep=".")]] <-
				getTargetSGP_INTERNAL(tmp_object_1[SGP_PROJECTION_GROUP==projection_group.iter], state, state.iter, projection_group.iter, target.type, target.level,
					year_within="YEAR_WITHIN" %in% names(sgp_object@Data))
		}
	} ### END for state.iter

	return(data.table(rbindlist(tmp.sgpTarget.list, fill=TRUE), key=getKey(sgp_object)))
} ### END getTargetSGP
