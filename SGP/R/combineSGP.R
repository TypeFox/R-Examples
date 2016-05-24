`combineSGP` <-
function(
	sgp_object,
	state=NULL,
	years=NULL,
	content_areas=NULL,
	sgp.percentiles=TRUE,
	sgp.percentiles.baseline=TRUE,
	sgp.projections=TRUE,
	sgp.projections.baseline=TRUE,
	sgp.projections.lagged=TRUE,
	sgp.projections.lagged.baseline=TRUE,
	sgp.target.scale.scores=FALSE,
	sgp.target.scale.scores.only=FALSE,
	sgp.target.content_areas=NULL,
	max.sgp.target.years.forward=3,
	update.all.years=FALSE,
	sgp.config=NULL,
	sgp.percentiles.equated=FALSE,
	SGPt=NULL,
	parallel.config=NULL) {

	started.at <- proc.time()
	message(paste("Started combineSGP", date()))

	ID <- CONTENT_AREA <- YEAR <- GRADE <- YEAR_INTEGER_TMP <- ACHIEVEMENT_LEVEL <- CATCH_UP_KEEP_UP_STATUS_INITIAL <- MOVE_UP_STAY_UP_STATUS_INITIAL <- VALID_CASE <- NULL
	MOVE_UP_STAY_UP_STATUS <- CATCH_UP_KEEP_UP_STATUS <- ACHIEVEMENT_LEVEL_PRIOR <- target.type <- SGP_PROJECTION_GROUP <- NULL

	tmp.messages <- NULL

	### Create slot.data from sgp_object@Data

	slot.data <- copy(sgp_object@Data)


	### Create state (if missing) from sgp_object (if possible)

	if (is.null(state)) {
		tmp.name <- toupper(gsub("_", " ", deparse(substitute(sgp_object))))
		state <- getStateAbbreviation(tmp.name, "combineSGP")
	}

	if (is.null(state)) {
		tmp.name <- toupper(gsub("_", " ", deparse(substitute(sgp_object))))
		tmp.name.position <- sapply(c(datasets::state.name, "AOB", "DEMONSTRATION"), function(x) regexpr(toupper(x), tmp.name))
		if (any(tmp.name.position!=-1)) {
			state <- c(datasets::state.abb, "AOB", "DEMO")[which(names(sort(tmp.name.position[tmp.name.position!=-1])[1])==c(datasets::state.name, "AOB", "DEMONSTRATION"))]
		} else {
			tmp.messages <- c(tmp.messages, "\tNOTE: argument 'state' required for target SGP calculation. Target SGPs will not be calculated.\n")
			sgp.projections.lagged <- sgp.projections.lagged.baseline <- FALSE
		}
	}

	### Create SGP_TARGET_CONTENT_AREA in certain cases

	if (is.null(sgp.target.content_areas) & any(sapply(SGP::SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]], function(x) length(unique(x))) > 1)) {
		sgp.target.content_areas <- TRUE
		tmp.messages <- c(tmp.messages, "\tNOTE: Multple content areas detected for student growth targets. 'sgp.target.content_areas set to TRUE.\n")
	}

	### Check to see if max.sgp.target.years.forward if configured in SGPstateData

	if (!is.null(SGP::SGPstateData[[state]][['SGP_Configuration']][['max.sgp.target.years.forward']])) {
		max.sgp.target.years.forward <- SGP::SGPstateData[[state]][['SGP_Configuration']][['max.sgp.target.years.forward']]
	}

	if (!is.null(SGP::SGPstateData[[state]][['SGP_Configuration']][['sgp.projections.projection.unit.label']])) {
		projection.unit.label <- SGP::SGPstateData[[state]][['SGP_Configuration']][['sgp.projections.projection.unit.label']]
	} else {
		projection.unit.label <- "YEAR"
	}


	### Setup for equated SGPs and scale score targets

	if (!is.null(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][["Year"]])) {
		year.for.equate <- tail(sort(unique(sgp_object@Data$YEAR)), 1)
		if (SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][["Year"]]!=year.for.equate) {
			sgp.percentiles.equated <- FALSE
			if (sgp.target.scale.scores) sgp.projections.equated <- NULL
		} else {
			sgp.percentiles.equated <- TRUE
			if (sgp.target.scale.scores) sgp.projections.equated <- list(Year=year.for.equate, Linkages=sgp_object@SGP[['Linkages']])
		}
	} else {
		if (sgp.percentiles.equated) {
			message("\tNOTE: 'sgp.percentiles.equated' has been set to TRUE but no meta-data exists in SGPstateData associated with the assessment transition. Equated/linked SGP analyses require meta-data embedded in 'SGPstateData' to correctly work. Contact package administrators on how such data can be added to the package.")
			sgp.percentiles.equated <- FALSE
			sgp.target.scale.scores <- FALSE
		}
		if (sgp.target.scale.scores) sgp.projections.equated <- NULL
	}


	### Utility functions

	get.target.arguments <- function(system.type, target.type=NULL, projection.unit.label) {
		tmp.list <- list()
		if (is.null(system.type)) {
			if (identical(target.type, c("sgp.projections", "sgp.projections.lagged"))) system.type <- "Cohort Referenced"
			if (identical(target.type, c("sgp.projections.baseline", "sgp.projections.lagged.baseline"))) system.type <- "Baseline Referenced"
			if (identical(target.type, c("sgp.projections", "sgp.projections.baseline", "sgp.projections.lagged", "sgp.projections.lagged.baseline"))) {
				system.type <- "Cohort and Baseline Referenced"
			}
		}

		if (!is.null(target.type)) {
			if (identical(target.type, "sgp.projections.lagged")) system.type <- "Cohort Referenced"
			if (identical(target.type, "sgp.projections.lagged.baseline")) system.type <- "Baseline Referenced"
			if (identical(target.type, c("sgp.projections.lagged", "sgp.projections.lagged.baseline"))) system.type <- "Cohort and Baseline Referenced"
		}

		if (identical(system.type, "Cohort Referenced")) {
			tmp.list[['target.type']] <- intersect(target.type, c("sgp.projections", "sgp.projections.lagged"))
			tmp.list[['my.sgp']] <- "SGP"
			tmp.list[['my.sgp.target']] <- paste("SGP_TARGET", max.sgp.target.years.forward, projection.unit.label, sep="_")
			tmp.list[['my.sgp.target.content_area']] <- paste("SGP_TARGET", max.sgp.target.years.forward, projection.unit.label, "CONTENT_AREA", sep="_")
			tmp.list[['my.sgp.target.move.up.stay.up']] <- paste("SGP_TARGET_MOVE_UP_STAY_UP", max.sgp.target.years.forward, projection.unit.label, sep="_")
			if (sgp.target.scale.scores) tmp.list[['sgp.target.scale.scores.types']] <- intersect(target.type, c("sgp.projections", "sgp.projections.lagged"))
		}
		if (identical(system.type, "Baseline Referenced")) {
			tmp.list[['target.type']] <- intersect(target.type, c("sgp.projections.baseline", "sgp.projections.lagged.baseline"))
			tmp.list[['my.sgp']] <- "SGP_BASELINE"
			tmp.list[['my.sgp.target']] <- paste("SGP_TARGET_BASELINE", max.sgp.target.years.forward, projection.unit.label, sep="_")
			tmp.list[['my.sgp.target.content_area']] <- paste("SGP_TARGET_BASELINE", max.sgp.target.years.forward, projection.unit.label, "CONTENT_AREA", sep="_")
			tmp.list[['my.sgp.target.move.up.stay.up']] <- paste("SGP_TARGET_BASELINE_MOVE_UP_STAY_UP", max.sgp.target.years.forward, projection.unit.label, sep="_")
			if (sgp.target.scale.scores) tmp.list[['sgp.target.scale.scores.types']] <- intersect(target.type, c("sgp.projections.baseline", "sgp.projections.lagged.baseline"))
		}
		if (identical(system.type, "Cohort and Baseline Referenced")) {
			tmp.list[['target.type']] <- intersect(target.type, c("sgp.projections", "sgp.projections.baseline", "sgp.projections.lagged", "sgp.projections.lagged.baseline"))
			tmp.list[['my.sgp']] <- c("SGP", "SGP_BASELINE")[c(sgp.percentiles, sgp.percentiles.baseline)]
			tmp.list[['my.sgp.target']] <- c(paste("SGP_TARGET", max.sgp.target.years.forward, projection.unit.label, sep="_"),
				paste("SGP_TARGET_BASELINE", max.sgp.target.years.forward, projection.unit.label, sep="_"))
			tmp.list[['my.sgp.target.content_area']] <- c(paste("SGP_TARGET", max.sgp.target.years.forward, projection.unit.label, "CONTENT_AREA", sep="_"),
				paste("SGP_TARGET_BASELINE", max.sgp.target.years.forward, projection.unit.label, "CONTENT_AREA", sep="_"))
			tmp.list[['my.sgp.target.move.up.stay.up']] <- c(paste("SGP_TARGET_MOVE_UP_STAY_UP", max.sgp.target.years.forward, projection.unit.label, sep="_"),
				paste("SGP_TARGET_BASELINE_MOVE_UP_STAY_UP", max.sgp.target.years.forward, projection.unit.label, sep="_"))
			if (sgp.target.scale.scores) tmp.list[['sgp.target.scale.scores.types']] <-
				intersect(target.type, c("sgp.projections", "sgp.projections.baseline", "sgp.projections.lagged", "sgp.projections.lagged.baseline"))
		}

		tmp.list[['target.level']] <- c("CATCH_UP_KEEP_UP", "MOVE_UP_STAY_UP")
		if (!is.null(SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Proficient"]]) &&
			length(which(SGP::SGPstateData[[state]][["Achievement"]][["Levels"]][["Proficient"]]=="Proficient")) <= 1) {
			tmp.list[['target.level']] <- "CATCH_UP_KEEP_UP"
		}
		if (!is.null(SGP::SGPstateData[[state]][["SGP_Configuration"]][['sgp.target.types']]) &&
			length(grep("MUSU", SGP::SGPstateData[[state]][["SGP_Configuration"]][['sgp.target.types']]))==0) {
				tmp.list[['target.level']] <- "CATCH_UP_KEEP_UP"
		}

		return(tmp.list)
	} ### END get.target.arguments

	catch_keep_move_functions <- c(min, max)

	"%w/o%" <- function(x,y) x[!x %in% y]

	getTargetData <- function(tmp.target.data, projection_group.iter, tmp.target.level.names) {
		tmp.data <- tmp.target.data[SGP_PROJECTION_GROUP==projection_group.iter, c("ID", "CONTENT_AREA", "YEAR", tmp.target.level.names), with=FALSE]
		return(tmp.data[apply(tmp.data[,tmp.target.level.names, with=FALSE],1,function(x)any(!is.na(x)))])
	}


	############################################################################
	### Check update.all.years
	############################################################################

	if (update.all.years) {
		variables.to.null.out <- c(
			"SGP", "SGP_SIMEX", "SGP_LEVEL", "SGP_STANDARD_ERROR", "SCALE_SCORE_PRIOR", "SCALE_SCORE_PRIOR_STANDARDIZED", "SGP_BASELINE", "SGP_LEVEL_BASELINE",
			"SGP_TARGET", "SGP_TARGET_MU", "SGP_TARGET_MU_BASELINE", "SGP_TARGET_MOVE_UP_STAY_UP", "SGP_TARGET_MOVE_UP_STAY_UP_BASELINE", "ACHIEVEMENT_LEVEL_PRIOR",
			"CATCH_UP_KEEP_UP_STATUS_INITIAL", "SGP_TARGET_BASELINE", "CATCH_UP_KEEP_UP_STATUS", "CATCH_UP_KEEP_UP_STATUS_BASELINE",
			"MOVE_UP_STATUS", "MOVE_UP_STAY_UP_STATUS", "MOVE_UP_STAY_UP_STATUS_BASELINE",
			"SGP_NORM_GROUP", "SGP_NORM_GROUP_BASELINE", "SGP_BASELINE_STANDARD_ERROR", "SGP_NORM_GROUP_SCALE_SCORES", "SGP_NORM_GROUP_BASELINE_SCALE_SCORES",
			grep("SGP_ORDER", names(slot.data), value=TRUE), grep("SGP_BASELINE_ORDER", names(slot.data), value=TRUE),
			grep("PERCENTILE_CUT", names(slot.data), value=TRUE), grep("CONFIDENCE_BOUND", names(slot.data), value=TRUE),
			paste("SGP_TARGET", max.sgp.target.years.forward, projection.unit.label, sep="_"),
			paste("SGP_TARGET_MOVE_UP_STAY_UP", max.sgp.target.years.forward, projection.unit.label, sep="_"),
			paste("SGP_TARGET", max.sgp.target.years.forward, projection.unit.label, "CURRENT", sep="_"),
			paste("SGP_TARGET_MOVE_UP_STAY_UP", max.sgp.target.years.forward, projection.unit.label, "CURRENT", sep="_"),
			paste("SGP_TARGET_BASELINE", max.sgp.target.years.forward, projection.unit.label, sep="_"),
			paste("SGP_TARGET_BASELINE_MOVE_UP_STAY_UP", max.sgp.target.years.forward, projection.unit.label, sep="_"),
			paste("SGP_TARGET_BASELINE", max.sgp.target.years.forward, projection.unit.label, "CURRENT", sep="_"),
			paste("SGP_TARGET_BASELINE_MOVE_UP_STAY_UP", max.sgp.target.years.forward, projection.unit.label, "CURRENT", sep="_"))

		for (tmp.variables.to.null.out in intersect(names(slot.data), variables.to.null.out)) {
			slot.data[,tmp.variables.to.null.out:=NULL, with=FALSE]
		}
	}


	############################################################################
	### sgp.percentiles: Merge Cohort Referenced SGPs with student data
	############################################################################

	## Determine names of Cohort Referenced SGPs

	if (!sgp.target.scale.scores.only && length(tmp.names <- getPercentileTableNames(sgp_object, content_areas, state, years, "sgp.percentiles", sgp.percentiles.equated)) == 0 && sgp.percentiles) {
		tmp.messages <- c(tmp.messages, "\tNOTE: No cohort referenced SGP results available in SGP slot. No cohort referenced SGP results will be merged.\n")
		sgp.percentiles <- FALSE
	}

	if (sgp.percentiles & !sgp.target.scale.scores.only) {

		tmp.list <- list()
		for (i in tmp.names) {
		tmp.list[[i]] <- data.table(
					CONTENT_AREA=unlist(strsplit(i, "[.]"))[1],
					YEAR=getTableNameYear(i),
					sgp_object@SGP[["SGPercentiles"]][[i]])
		}

		tmp.data <- data.table(rbindlist(tmp.list, fill=TRUE), VALID_CASE="VALID_CASE", key=key(slot.data))

		if (any(duplicated(tmp.data))) {
			tmp.data <- getPreferredSGP(tmp.data, state)
		}

		variables.to.merge <- names(tmp.data) %w/o% key(slot.data)
		tmp.index <- slot.data[tmp.data[,key(slot.data), with=FALSE], which=TRUE]
		slot.data[tmp.index, variables.to.merge := tmp.data[, variables.to.merge, with=FALSE], with=FALSE, nomatch=0]

		setkeyv(slot.data, getKey(slot.data))
	}


	###################################################################################
	### sgp.percentiles.baseline: Merge baseline referenced SGPs with student data
	###################################################################################

	## Determine names of Baseline Referenced SGPs

	if (!sgp.target.scale.scores.only && length(tmp.names <- getPercentileTableNames(sgp_object, content_areas, state, years, "sgp.percentiles.baseline"))==0 && sgp.percentiles.baseline) {
		 tmp.messages <- c(tmp.messages, "\tNOTE: No baseline referenced SGP results available in SGP slot. No baseline referenced SGP results will be merged.\n")
		 sgp.percentiles.baseline <- FALSE
	}

	if (sgp.percentiles.baseline & !sgp.target.scale.scores.only) {

		tmp.list <- list()
		for (i in tmp.names) {
			tmp.list[[i]] <- data.table(
				CONTENT_AREA=unlist(strsplit(i, "[.]"))[1],
				YEAR=getTableNameYear(i),
				sgp_object@SGP[["SGPercentiles"]][[i]])

			if (is.na(unlist(strsplit(i, "[.]"))[3])) { ### If cohort referenced SGP are to be included in baseline SGP (e.g., Georgia)
				setnames(tmp.list[[i]], "SGP", "SGP_BASELINE")
				if ("SGP_LEVEL" %in% names(tmp.list[[i]])) setnames(tmp.list[[i]], "SGP_LEVEL", "SGP_LEVEL_BASELINE")
				if ("SGP_NORM_GROUP" %in% names(tmp.list[[i]])) setnames(tmp.list[[i]], "SGP_NORM_GROUP", "SGP_NORM_GROUP_BASELINE")
				if ("SGP_NORM_GROUP_SCALE_SCORES" %in% names(tmp.list[[i]])) setnames(tmp.list[[i]], "SGP_NORM_GROUP_SCALE_SCORES", "SGP_NORM_GROUP_BASELINE_SCALE_SCORES")
				if ("SGP_SIMEX" %in% names(tmp.list[[i]])) setnames(tmp.list[[i]], "SGP_SIMEX", "SGP_SIMEX_BASELINE")
			}
		}

		tmp.data <- data.table(rbindlist(tmp.list, fill=TRUE), VALID_CASE="VALID_CASE", key=key(slot.data))

		if (any(duplicated(tmp.data))) {
			tmp.data <- getPreferredSGP(tmp.data, state, type="BASELINE")
		}

		variables.to.merge <- names(tmp.data) %w/o% key(slot.data)
		tmp.index <- slot.data[tmp.data[,key(slot.data), with=FALSE], which=TRUE]
		slot.data[tmp.index, variables.to.merge := tmp.data[, variables.to.merge, with=FALSE], with=FALSE, nomatch=0]

		setkeyv(slot.data, getKey(slot.data))
	}


	######################################################################################
	### Create SGP targets (Cohort and Baseline referenced) and merge with student data
	######################################################################################


	if (!sgp.target.scale.scores.only && length(getPercentileTableNames(sgp_object, content_areas, state, years, "sgp.projections"))==0 && sgp.projections) {
		tmp.messages <- c(tmp.messages, "\tNOTE: No SGP projections available in SGP slot. No current year student growth projection targets will be produced.\n")
		sgp.projections <- FALSE;
	}
	if (!sgp.target.scale.scores.only && length(getPercentileTableNames(sgp_object, content_areas, state, years, "sgp.projections.baseline"))==0 && sgp.projections.baseline) {
		tmp.messages <- c(tmp.messages, "\tNOTE: No SGP baseline projections available in SGP slot. No current year baseline student growth projection targets will be produced.\n")
		sgp.projections.baseline <- FALSE;
	}
	if (!sgp.target.scale.scores.only && length(getPercentileTableNames(sgp_object, content_areas, state, years, "sgp.projections.lagged"))==0 && sgp.projections.lagged) {
		tmp.messages <- c(tmp.messages, "\tNOTE: No SGP lagged projections available in SGP slot. No student growth projection targets will be produced.\n")
		sgp.projections.lagged <- FALSE;
	}
	if (!sgp.target.scale.scores.only && length(getPercentileTableNames(sgp_object, content_areas, state, years, "sgp.projections.lagged.baseline"))==0 && sgp.projections.lagged.baseline) {
		tmp.messages <- c(tmp.messages, "\tNOTE: No SGP lagged baseline projections available in SGP slot. No baseline referenced student growth projection targets will be produced.\n")
		sgp.projections.lagged.baseline <- FALSE;
	}
	target.type <- c("sgp.projections", "sgp.projections.baseline", "sgp.projections.lagged", "sgp.projections.lagged.baseline")[
				c(sgp.projections, sgp.projections.baseline, sgp.projections.lagged, sgp.projections.lagged.baseline)]

	### Calculate Targets

	if ((sgp.projections | sgp.projections.baseline | sgp.projections.lagged | sgp.projections.lagged.baseline) & !sgp.target.scale.scores.only) {

		target.args <- get.target.arguments(SGP::SGPstateData[[state]][["Growth"]][["System_Type"]], target.type, projection.unit.label)

		for (target.type.iter in target.args[['target.type']]) {
			for (target.level.iter in target.args[['target.level']]) {
				tmp.data <- getTargetSGP(sgp_object, content_areas, state, years, target.type.iter, target.level.iter, max.sgp.target.years.forward)
				if (any(duplicated(tmp.data))) {
					duplicated.projections.tf <- TRUE
					tmp.data <- getPreferredSGP(tmp.data, state, type="TARGET")
				} else duplicated.projections.tf <- FALSE
				variables.to.merge <- names(tmp.data) %w/o% key(slot.data)
				tmp.index <- slot.data[tmp.data[,intersect(key(slot.data), names(tmp.data)), with=FALSE], which=TRUE]
				slot.data[tmp.index, variables.to.merge := tmp.data[, variables.to.merge, with=FALSE], with=FALSE, nomatch=0]
			}
		}
		if (duplicated.projections.tf) {
			tmp.messages <- c(tmp.messages, paste(
				"\tNOTE: Multiple Projections exist for individual students. Unique SGP Targets will be created using SGP Progression Preference Table for ", state, ".\n", sep=""))
		}

		### SGP_TARGET_CONTENT_AREA calculation

		terminal.content_areas <- unique(slot.data[!is.na(target.args[['my.sgp.target']][1])][['CONTENT_AREA']])
		if (!is.null(SGP::SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]])) {
			terminal.content_areas <- intersect(terminal.content_areas, sapply(SGP::SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]], tail, 1))
		}

		if (!is.null(sgp.target.content_areas) && sgp.target.content_areas) {
			for (my.sgp.target.content_area.iter in target.args[['my.sgp.target.content_area']]) {
				slot.data[!is.na(get(target.args[['my.sgp.target']][1])), target.args[['my.sgp.target.content_area']] :=
					getTargetSGPContentArea(GRADE[1], CONTENT_AREA[1], state, max.sgp.target.years.forward, target.args[['my.sgp.target.content_area.iter']]),
					by=list(GRADE, CONTENT_AREA), with=FALSE]
			}
		}

		### CATCH_UP_KEEP_UP_STATUS Calculation

		if ("CATCH_UP_KEEP_UP" %in% target.args[['target.level']] & (sgp.projections.lagged | sgp.projections.lagged.baseline)) {

			catch.up.keep.up.levels <- getTargetAchievementLevels(state, "CATCH_UP_KEEP_UP")
			slot.data[,CATCH_UP_KEEP_UP_STATUS_INITIAL:=getTargetInitialStatus(ACHIEVEMENT_LEVEL_PRIOR, state, status.type="CATCH_UP_KEEP_UP")]

			for (i in seq_along(target.args[['my.sgp']])) {
				if (length(grep("BASELINE", target.args[['my.sgp']][i]))==0) my.label <- "CATCH_UP_KEEP_UP_STATUS" else my.label <- "CATCH_UP_KEEP_UP_STATUS_BASELINE"
				if (my.label %in% names(slot.data)) slot.data[,my.label := NULL, with=FALSE]
				slot.data[,my.label := rep(as.character(NA), dim(slot.data)[1]), with=FALSE]

				slot.data[CATCH_UP_KEEP_UP_STATUS_INITIAL == "Keeping Up" & get(target.args[['my.sgp']][i]) >= get(target.args[['my.sgp.target']][i]),
					my.label := "Keep Up: Yes", with=FALSE]
				slot.data[CATCH_UP_KEEP_UP_STATUS_INITIAL == "Keeping Up" & get(target.args[['my.sgp']][i]) < get(target.args[['my.sgp.target']][i]),
					my.label := "Keep Up: No", with=FALSE]
				slot.data[CATCH_UP_KEEP_UP_STATUS_INITIAL == "Catching Up" & get(target.args[['my.sgp']][i]) >= get(target.args[['my.sgp.target']][i]),
					my.label := "Catch Up: Yes", with=FALSE]
				slot.data[CATCH_UP_KEEP_UP_STATUS_INITIAL == "Catching Up" & get(target.args[['my.sgp']][i]) < get(target.args[['my.sgp.target']][i]),
					my.label := "Catch Up: No", with=FALSE]

				### CATCH_UP_KEEP_UP clean up based upon reality

				slot.data[CATCH_UP_KEEP_UP_STATUS_INITIAL == "Keeping Up" & get(my.label) == "Keep Up: Yes" &
					ACHIEVEMENT_LEVEL %in% catch.up.keep.up.levels[['NO']], my.label := "Keep Up: No", with=FALSE]
				slot.data[CATCH_UP_KEEP_UP_STATUS_INITIAL == "Catching Up" & get(my.label) == "Catch Up: No" &
					ACHIEVEMENT_LEVEL %in% catch.up.keep.up.levels[['YES']], my.label := "Catch Up: Yes", with=FALSE]
				slot.data[CATCH_UP_KEEP_UP_STATUS_INITIAL == "Catching Up" & get(my.label) == "Catch Up: Yes" &
					ACHIEVEMENT_LEVEL %in% catch.up.keep.up.levels[['NO']] &
					GRADE == max(type.convert(GRADE[!is.na(get(target.args[['my.sgp.target']]))], as.is=TRUE)) &
					CONTENT_AREA %in% terminal.content_areas, my.label := "Catch Up: No", with=FALSE]
				slot.data[CATCH_UP_KEEP_UP_STATUS_INITIAL == "Keeping Up" & get(my.label) == "Keep Up: No" &
					ACHIEVEMENT_LEVEL %in% catch.up.keep.up.levels[['YES']] &
					GRADE == max(type.convert(GRADE[!is.na(get(target.args[['my.sgp.target']]))], as.is=TRUE)) &
					CONTENT_AREA %in% terminal.content_areas, my.label := "Keep Up: Yes", with=FALSE]
				slot.data[,my.label := as.factor(get(my.label)), with=FALSE]
			}
		}


		### MOVE_UP_STAY_UP_STATUS Calculation

		if ("MOVE_UP_STAY_UP" %in% target.args[['target.level']] & (sgp.projections.lagged | sgp.projections.lagged.baseline)) {

			move.up.stay.up.levels <- getTargetAchievementLevels(state, "MOVE_UP_STAY_UP")
			slot.data[,MOVE_UP_STAY_UP_STATUS_INITIAL:=getTargetInitialStatus(ACHIEVEMENT_LEVEL_PRIOR, state, status.type="MOVE_UP_STAY_UP")]

			for (i in seq_along(target.args[['my.sgp']])) {
				if (length(grep("BASELINE", target.args[['my.sgp']][i]))==0) my.label <- "MOVE_UP_STAY_UP_STATUS" else my.label <- "MOVE_UP_STAY_UP_STATUS_BASELINE"
				if (my.label %in% names(slot.data)) slot.data[,my.label := NULL, with=FALSE]
				slot.data[,my.label := rep(as.character(NA), dim(slot.data)[1]), with=FALSE]

				slot.data[MOVE_UP_STAY_UP_STATUS_INITIAL == "Staying Up" & get(target.args[['my.sgp']][i]) >= get(target.args[['my.sgp.target.move.up.stay.up']][i]),
					my.label := "Stay Up: Yes", with=FALSE]
				slot.data[MOVE_UP_STAY_UP_STATUS_INITIAL == "Staying Up" & get(target.args[['my.sgp']][i]) < get(target.args[['my.sgp.target.move.up.stay.up']][i]),
					my.label := "Stay Up: No", with=FALSE]
				slot.data[MOVE_UP_STAY_UP_STATUS_INITIAL == "Moving Up" & get(target.args[['my.sgp']][i]) >= get(target.args[['my.sgp.target.move.up.stay.up']][i]),
					my.label := "Move Up: Yes", with=FALSE]
				slot.data[MOVE_UP_STAY_UP_STATUS_INITIAL == "Moving Up" & get(target.args[['my.sgp']][i]) < get(target.args[['my.sgp.target.move.up.stay.up']][i]),
					my.label := "Move Up: No", with=FALSE]

				### MOVE_UP_STAY_UP clean up based upon reality

				slot.data[MOVE_UP_STAY_UP_STATUS_INITIAL == "Staying Up" & get(my.label) == "Stay Up: Yes" &
					ACHIEVEMENT_LEVEL %in% move.up.stay.up.levels[['NO']], my.label := "Stay Up: No", with=FALSE]
				slot.data[MOVE_UP_STAY_UP_STATUS_INITIAL == "Moving Up" & get(my.label) == "Move Up: No" &
					ACHIEVEMENT_LEVEL %in% move.up.stay.up.levels[['YES']], my.label := "Move Up: Yes", with=FALSE]
				slot.data[MOVE_UP_STAY_UP_STATUS_INITIAL == "Moving Up" & get(my.label) == "Move Up: Yes" &
					ACHIEVEMENT_LEVEL %in% move.up.stay.up.levels[['NO']] &
					GRADE == max(type.convert(GRADE[!is.na(get(target.args[['my.sgp.target.move.up.stay.up']]))], as.is=TRUE)) &
					CONTENT_AREA %in% terminal.content_areas, my.label := "Move Up: No", with=FALSE]
				slot.data[MOVE_UP_STAY_UP_STATUS_INITIAL == "Staying Up" & get(my.label) == "Stay Up: No" &
					ACHIEVEMENT_LEVEL %in% move.up.stay.up.levels[['YES']] &
					GRADE == max(type.convert(GRADE[!is.na(get(target.args[['my.sgp.target.move.up.stay.up']]))], as.is=TRUE)) &
					CONTENT_AREA %in% terminal.content_areas, my.label := "Stay Up: Yes", with=FALSE]
				slot.data[,my.label := as.factor(get(my.label)), with=FALSE]
			}
		}

		for (i in intersect(names(slot.data), c("CATCH_UP_KEEP_UP_STATUS_INITIAL", "MOVE_UP_STAY_UP_STATUS_INITIAL"))) {
			slot.data[,i:=NULL,with=FALSE]
		}

	} ## END sgp.projections.lagged | sgp.projections.lagged.baseline


	###################################################################################################
	### Create SGP Scale Score targets (Cohort and Baseline referenced) if requested
	###################################################################################################

	if (sgp.target.scale.scores) {

		if (!exists("target.args")) target.args <- get.target.arguments(SGP::SGPstateData[[state]][["Growth"]][["System_Type"]], target.type, projection.unit.label)
		tmp.target.list <- list()
		for (target.type.iter in target.args[['sgp.target.scale.scores.types']]) {
			for (target.level.iter in target.args[['target.level']]) {
				tmp.target.list[[paste(target.type.iter, target.level.iter)]] <-
					data.table(getTargetSGP(sgp_object, content_areas, state, years, target.type.iter, target.level.iter, max.sgp.target.years.forward, return.lagged.status=FALSE),
						key=c(getKey(sgp_object), "SGP_PROJECTION_GROUP"))
			}
		}
		tmp.target.data <- data.table(Reduce(function(x, y) merge(x, y, all=TRUE), tmp.target.list[!sapply(tmp.target.list, function(x) dim(x)[1]==0)],
			accumulate=FALSE), key=c("VALID_CASE", "CONTENT_AREA", "YEAR", "ID"))

		for (projection_group.iter in unique(tmp.target.data[['SGP_PROJECTION_GROUP']])) {
			for (target.type.iter in target.args[['sgp.target.scale.scores.types']]) {
				tmp.target.level.names <-
					as.character(sapply(target.args[['target.level']], function(x) getTargetName(state, target.type.iter, x, max.sgp.target.years.forward, "SGP_TARGET", projection.unit.label, projection_group.iter)))
				if (any(!tmp.target.level.names %in% names(tmp.target.data))) {
					tmp.target.data[,tmp.target.level.names[!tmp.target.level.names %in% names(tmp.target.data)]:=as.integer(NA)]
				}

				sgp_object <- getTargetScaleScore(
					sgp_object,
					state,
					getTargetData(tmp.target.data, projection_group.iter, tmp.target.level.names),
					target.type.iter,
					tmp.target.level.names,
					getYearsContentAreasGrades(state, years=unique(tmp.target.data[SGP_PROJECTION_GROUP==projection_group.iter][['YEAR']]), content_areas=unique(tmp.target.data[SGP_PROJECTION_GROUP==projection_group.iter][['CONTENT_AREA']])),
					sgp.config=sgp.config,
					projection_group.identifier=projection_group.iter,
					sgp.projections.equated=if (length(grep("baseline", target.type.iter) > 0)) NULL else sgp.projections.equated,
					SGPt=SGPt,
					parallel.config=parallel.config)
			}
		}
	} ### END if (sgp.target.scale.scores)


	### Put slot.data into @Data slot

	setkeyv(slot.data, getKey(slot.data))
	sgp_object@Data <- slot.data

	message(c(tmp.messages, paste("Finished combineSGP", date(), "in", convertTime(timetaken(started.at)), "\n"), sep=""))

	return(sgp_object)
} ## END combineSGP Function
