`prepareSGP` <-
function(data,
	data_supplementary=NULL,
	state=NULL,
	var.names=NULL,
	create.additional.variables=TRUE,
	fix.duplicates=NULL) {

	VALID_CASE <- ID <- CONTENT_AREA <- YEAR <- ID <- GRADE <- SCALE_SCORE <- DUPLICATED_CASES <- tmp.dups.index.to.remove <- NULL
	SGPstateData <- SGP::SGPstateData ### Needed due to possible assignment of values to SGPstateData

	## Print start time

	started.at <- proc.time()
	messageSGP(paste("\nStarted prepareSGP", date()), "\n")
	messageSGP(match.call())

	## Get state (if possible)

		if (is.null(state)) {
			tmp.name <- toupper(gsub("_", " ", deparse(substitute(data))))
			state <- getStateAbbreviation(tmp.name, "prepareSGP")
		}


	### Utility functions

	## getNames

	getNames <- function(names.provided, var.names) {

		## Get the names of the provided variables and create data frame

		variable.names.provided <- data.frame(
			column.provided=seq_along(names.provided),
			names.provided=names.provided,
			stringsAsFactors=FALSE)

		## Required variables and default data frame

		required.names <- c("ID", "CONTENT_AREA", "YEAR", "GRADE", "SCALE_SCORE", "VALID_CASE")

		default.variable.names <- data.frame(
			names.provided=required.names,
			names.sgp=required.names,
			names.type=SGPstateData[["DEMO"]][["Variable_Name_Lookup"]]$names.sgp.type[match(required.names, SGPstateData[["DEMO"]][["Variable_Name_Lookup"]]$names.sgp)],
			names.info=SGPstateData[["DEMO"]][["Variable_Name_Lookup"]]$names.sgp.info[match(required.names, SGPstateData[["DEMO"]][["Variable_Name_Lookup"]]$names.sgp)],
			names.output=SGPstateData[["DEMO"]][["Variable_Name_Lookup"]]$names.sgp.output[match(required.names, SGPstateData[["DEMO"]][["Variable_Name_Lookup"]]$names.sgp)],
			stringsAsFactors=FALSE)

		## Check names.provided

		if (!is.null(var.names)) {
			if (!class(var.names) %in% c("list", "data.frame")) {
				stop("\tNOTE: Supplied argument to 'var.names' must be of class list or data.frame. Please supply argument of appropriate class.")
			}
			if (identical("data.frame", class(var.names))) {
				if (!identical(names(var.names), c("names.provided", "names.sgp", "names.type", "names.info", "names.output"))) {
					stop("\tNOTE: Supplied data.frame to arugment 'var.names' does not include all required columns: 'names.provided', 'names.sgp', 'names.type', 'names.info', 'names.output'")
				}
			}
			if (identical("list", class(var.names))) {
				var.names <- data.frame(
					names.provided=unlist(var.names),
					names.sgp=toupper(names(var.names)),
					names.type=SGPstateData[["DEMO"]][["Variable_Name_Lookup"]]$names.sgp.type[
						match(toupper(names(var.names)), SGPstateData[["DEMO"]][["Variable_Name_Lookup"]]$names.sgp)],
					names.info=SGPstateData[["DEMO"]][["Variable_Name_Lookup"]]$names.sgp.info[
						match(toupper(names(var.names)), SGPstateData[["DEMO"]][["Variable_Name_Lookup"]]$names.sgp)],
					names.output=SGPstateData[["DEMO"]][["Variable_Name_Lookup"]]$names.sgp.output[
						match(toupper(names(var.names)), SGPstateData[["DEMO"]][["Variable_Name_Lookup"]]$names.sgp)],
					stringsAsFactors=FALSE, row.names=NULL)
			}

			## Include default variable names (as needed)

			tmp.names <- default.variable.names[!(default.variable.names[["names.sgp"]] %in% var.names[["names.sgp"]]),]
			variable.names <- rbind(var.names, tmp.names)
		} else {
			if (is.null(state) || state=="DEMO" || is.null(SGPstateData[[state]][["Variable_Name_Lookup"]])) {
				variable.names <- default.variable.names
			} else {
				variable.names <- SGPstateData[[state]][["Variable_Name_Lookup"]]
			}
		}

		## Compile the provided variable names and the corresponding (capitalized) variable names used in SGP

		variable.names <- merge(variable.names.provided, variable.names, by.x="names.provided", by.y="names.provided", all.x=TRUE)
		variable.names <- subset(variable.names, select=c("names.provided", "names.sgp", "names.type", "names.info", "names.output", "column.provided"))
		variable.names <- merge(variable.names, SGPstateData[["DEMO"]][["Variable_Name_Lookup"]], by.x="names.provided", by.y="names.sgp", all.x=TRUE)
		variable.names[is.na(variable.names$names.type) & !is.na(variable.names$names.sgp.type),][,c("names.sgp", "names.type", "names.info", "names.output")] <-
			variable.names[is.na(variable.names$names.type) & !is.na(variable.names$names.sgp.type),][,c("names.provided", "names.sgp.type", "names.sgp.info", "names.sgp.output")]
		variable.names$names.sgp.type <- variable.names$names.sgp.info <- variable.names$names.sgp.output <- NULL

		## Check see if any of the required variables are missing

		if (!all(required.names %in% variable.names$names.sgp)) {
			stop(paste("\tNOTE: The {data} object is missing the following column name(s): ", paste(required.names[(required.names %in% variable.names$names.sgp)==FALSE], collapse=", "),
			". Please identify the variable using the {var.names} argument.", sep=""))
		}
		return(data.frame(variable.names[order(variable.names$column.provided),][,c("names.provided", "names.sgp", "names.type", "names.info", "names.output")], row.names=NULL, stringsAsFactors=FALSE))
	} ## END getNames


	###################################################################
	###
	### prepare object depending upon whether it is of class SGP
	###
	###################################################################

	if (is.SGP(data)) {

		if (!is.null(state) & is.null(var.names)) {
			if (!identical(state, "DEMO") & !identical(data@Names, SGPstateData[[state]][["Variable_Name_Lookup"]])) {
				data@Names <- SGPstateData[[state]][["Variable_Name_Lookup"]]
			}
			if (identical(state, "DEMO") & !identical(data@Names, SGPstateData[[state]][["Variable_Name_Lookup"]])) {
				data@Names <- getNames(SGPstateData[['DEMO']][['Variable_Name_Lookup']][['names.sgp']], var.names)
			}
		}

		if (!is.null(var.names)) {
			data@Names <- getNames(names(data@Data), var.names)
		}

		## run checkSGP

		data <- checkSGP(data, state=state)

		## Key data.table and check for duplicate cases

		setkeyv(data@Data, getKey(data))
		if (is.null(fix.duplicates) && any(duplicated(data@Data["VALID_CASE"]))) {
			message(paste("\tWARNING: @Data keyed by", paste(getKey(data), collapse=", "), "has duplicate cases. Subsequent joins/merges will likely be corrupt."))
			message("\tNOTE: Duplicate cases are available in current workspace as 'DUPLICATED_CASES' and saved as 'DUPLICATED_CASES.Rdata'.")
			assign("DUPLICATED_CASES",
				data@Data["VALID_CASE"][unique(data@Data["VALID_CASE"][duplicated(data@Data["VALID_CASE"]), c("VALID_CASE", "CONTENT_AREA", "YEAR", "ID"), with=FALSE])])
			save(DUPLICATED_CASES, file="DUPLICATED_CASES.Rdata")
		}

		## Check for knots and boundaries

		if (is.null(SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]])) {
			SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]] <- createKnotsBoundaries(data@Data)
			assign(paste(state, "Knots_Boundaries", sep="_"), SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]])
			save(list=paste(state, "Knots_Boundaries", sep="_"), file=paste(state, "Knots_Boundaries.Rdata", sep="_"))
			message(paste("\tNOTE: Knots and Boundaries do not exist for state provided.\n\tThey have been produced and are available using state=", state, " for subsequent analyses and saved to your working directory '", getwd(), "'.", sep=""))
		}

		## Create FIRST_OBESRVATION, LAST_OBSERVATION if YEAR_WITHIN exists

		if ("YEAR_WITHIN" %in% names(data@Data) & !all(c("FIRST", "LAST") %in% names(data@Data))) {
			data@Data <- getFirstAndLastInGroup(data@Data)
			setkeyv(data@Data, getKey(data))
		}

		if (!.hasSlot(data, "Data_Supplementary")) data@Data_Supplementary <- NULL

		if (!is.null(data_supplementary)) {
			if (!identical(class(data_supplementary), "list")) {
				stop("\tNOTE: Supplied supplementary data to 'data_supplementary' must be data table(s) embedded in a wrapper list")
			} else {
				data@Data_Supplementary <- c(data@Data_Supplementary, data_supplementary)
			}
		}

		data@Version <- getVersion(data)
		sgp_object <- data
	} else {
		variable.names <- getNames(names(data), var.names)

		##  Create keyed data.table and check for duplicate cases

		data <- as.data.table(data)
		setnames(data, which(!is.na(variable.names$names.sgp)), variable.names$names.sgp[!is.na(variable.names$names.sgp)])
		setkeyv(data, getKey(data))
		if (is.null(fix.duplicates) && any(duplicated(data["VALID_CASE"]))) {
			message(paste("\tWARNING: Data keyed by", paste(getKey(data), collapse=", "), "has duplicate cases. Subsequent joins/merges will be corrupted."))
			message("\tNOTE: Duplicate cases are available in current workspace as 'DUPLICATED_CASES' and saved as 'DUPLICATED_CASES.Rdata'.")
			assign("DUPLICATED_CASES", data["VALID_CASE"][unique(data["VALID_CASE"][duplicated(data["VALID_CASE"]), c("VALID_CASE", "CONTENT_AREA", "YEAR", "ID"), with=FALSE])])
			save(DUPLICATED_CASES, file="DUPLICATED_CASES.Rdata")
		}

		## Check for knots and boundaries

		if (is.null(SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]])) {
			SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]] <- createKnotsBoundaries(data)
			assign(paste(state, "Knots_Boundaries", sep="_"), SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]])
			save(list=paste(state, "Knots_Boundaries", sep="_"), file=paste(state, "Knots_Boundaries.Rdata", sep="_"))
			save(SGPstateData, file="SGPstateData.Rdata")
			message(paste("\tNOTE: Knots and Boundaries do not exist for state provided.\n\tThey have been produced and are available using state=", state, " for subsequent analyses and saved to your working directory '", getwd(), "'.", sep=""))
		}

		## Create FIRST_OBESRVATION, LAST_OBSERVATION if YEAR_WITHIN exists

		if ("YEAR_WITHIN" %in% names(data) & !all(c("FIRST", "LAST") %in% names(data))) {
			data <- getFirstAndLastInGroup(data)
			setkeyv(data, getKey(data))
		}

		## Test data_supplementary argument

		if (!is.null(data_supplementary) && !identical(class(data_supplementary), "list")) {
			stop("\tNOTE: Supplied supplementary data to 'data_supplementary' must be data table(s) embedded in a wrapper list")
		}

		##  Create the SGP object

		sgp_object <- new("SGP", Data=data, Data_Supplementary=data_supplementary, Names=variable.names, Version=getVersion(data))
		sgp_object <- checkSGP(sgp_object, state=state)

	} ## END else


	#################################################################
	###
	### Add additional variables, correct duplicates, ...
	###
	#################################################################

	## Create ACHIEVEMENT_LEVEL is it doesn't exist

	if (!"ACHIEVEMENT_LEVEL" %in% names(sgp_object@Data) &
			(!is.null(SGPstateData[[state]][["Achievement"]][["Cutscores"]]) | !is.null(SGPstateData[[state]][["Achievement"]][["Cutscore_Information"]]))) {
		sgp_object@Data <- getAchievementLevel(sgp_object@Data, state=state)
		setkeyv(sgp_object@Data, getKey(sgp_object))
		message(paste("\tNOTE: Added variable ACHIEVEMENT_LEVEL to @Data using", state, "cutscores embedded in SGPstateData."))
	}

	if (create.additional.variables) {

		### HIGH_NEED_STATUS

		if (!"HIGH_NEED_STATUS" %in% names(sgp_object@Data) & "SCHOOL_NUMBER" %in% names(sgp_object@Data)) {
			sgp_object <- getHighNeedStatus(sgp_object)
		}

		### STATE_ENROLLMENT_STATUS, DISTRICT_ENROLLMENT_STATUS, SCHOOL_ENROLLMENT_STATUS

		tmp.enrollment_status.variables <- c("STATE_ENROLLMENT_STATUS", "DISTRICT_ENROLLMENT_STATUS", "SCHOOL_ENROLLMENT_STATUS")
		tmp.enrollment_status.levels <- c("STATE", "DISTRICT", "SCHOOL")

		for (i in seq_along(tmp.enrollment_status.variables)) {
			if (!tmp.enrollment_status.variables[i] %in% names(sgp_object@Data) & paste(tmp.enrollment_status.levels[i], "NUMBER", sep="_") %in% names(sgp_object@Data)) {
				sgp_object@Data[[tmp.enrollment_status.variables[i]]] <-
					factor(1, levels=0:1, labels=paste("Enrolled", capwords(tmp.enrollment_status.levels[i]), c(": No", ": Yes"), sep=""))
				sgp_object@Data[[tmp.enrollment_status.variables[i]]][sgp_object@Data[['VALID_CASE']]!="VALID_CASE"] <- NA
				message(paste("\tNOTE: Added variable", tmp.enrollment_status.variables[i], "to @Data."))
			}
		}

	} ### end if (create.additional.variables)

	if (!is.null(fix.duplicates) && !is.null(DUPLICATED_CASES)) {
		if (identical(toupper(fix.duplicates), "KEEP.ALL")) {
			if (all(unique(DUPLICATED_CASES$YEAR) %in% (tmp.last.year <- tail(sort(unique(sgp_object@Data$YEAR)), 1)))) {
				sgp_object@Data <- createUniqueLongData(sgp_object@Data)
				message("\tNOTE: Duplicate cases in current year made UNIQUE. Modified IDs include suffix '_DUPS_***' in @Data.")
			} else {
				message("\tNOTE: Duplicate case modification is only available when duplicates reside in last year of data. Duplicate cases are NOT fixed.")
			}
		} ### End KEEP.ALL
	}

	##  Print finish time
	message(paste("Finished prepareSGP", date(), "in", convertTime(timetaken(started.at)), "\n"))

	return(sgp_object)
} ## END prepareSGP function
