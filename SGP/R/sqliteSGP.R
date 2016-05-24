`sqliteSGP` <-
function(sgp_object,
	state=NULL,
	years=NULL,
	content_areas=NULL,
	other.student.groups,
	text.output=TRUE,
	json.output=TRUE,
	null.output.string=NULL,
	projection.years.for.target=3,
	output.directory=file.path("Data", "SchoolView")) {

	started.at <- proc.time()
	message(paste("\tStarted sqliteSGP in outputSGP", date()))

	YEAR <- DISTRICT_NUMBER <- SCHOOL_NUMBER <- CONTENT_AREA <- DISTRICT_ENROLLMENT_STATUS <- GRADE <- ETHNICITY <- STUDENTGROUP <- SCHOOL_ENROLLMENT_STATUS <- EMH_LEVEL <- MEDIAN_SGP <- NULL
	INSTRUCTOR_NUMBER <- INSTRUCTOR_ENROLLMENT_STATUS <- TMP_ID <- NULL

        ## Create state (if NULL) from sgp_object (if possible)

	        if (is.null(state)) {
			tmp.name <- toupper(gsub("_", " ", deparse(substitute(sgp_object))))
			state <- getStateAbbreviation(tmp.name, "sqliteSGP")
	        }

	## Create group variable names

		if (!is.null(SGP::SGPstateData[[state]][["SGP_Configuration"]][["output.groups"]])) {
			output.groups <- SGP::SGPstateData[[state]][["SGP_Configuration"]][["output.groups"]]
		} else {
			output.groups <- c("DISTRICT", "SCHOOL")
		}

		group.number <- paste(output.groups, "NUMBER", sep="_")
		group.enroll.status <- paste(output.groups, "ENROLLMENT_STATUS", sep="_")
		group.enroll.status.label <- paste("Enrolled ", sapply(output.groups, capwords), ": Yes", sep="")


	## Create/Set database

		if (state %in% c(datasets::state.abb, "DEMO")) {
			tmp.state <- gsub(" ", "_", c(datasets::state.name, "Demonstration")[state==c(datasets::state.abb, "DEMO")])
		} else {
			tmp.state <- gsub(" ", "_", state)
		}

		sqlite.output.directory <- file.path(output.directory, "SQLITE")
		dir.create(sqlite.output.directory, recursive=TRUE, showWarnings=FALSE)
		if (file.exists(file.path(sqlite.output.directory, paste(tmp.state, "_Data_SQLITE.sqlite", sep="")))) file.remove(file.path(sqlite.output.directory, paste(tmp.state, "_Data_SQLITE.sqlite", sep="")))
		db <- dbConnect(SQLite(), dbname=file.path(sqlite.output.directory, paste(tmp.state, "_Data_SQLITE.sqlite", sep="")))
		if (text.output) {
			text.output.directory <- file.path(output.directory, "TEXT")
			dir.create(text.output.directory, recursive=TRUE, showWarnings=FALSE)
		}
		if (json.output) {
			json.output.directory <- file.path(output.directory, "JSON")
			dir.create(json.output.directory, recursive=TRUE, showWarnings=FALSE)
		}


	## Utility functions

		strtail <- function(s, n=1) {
			if (n < 0) substring(s, 1-n)
			else substring(s, nchar(s)-n+1)
		}

		strhead <- function(s,n=1) {
			if (n < 0) substr(s, 1, nchar(s)+n)
			else substr(s, 1, n)
		}

		sqlite.create.table <- function(table.name, field.types, primary.key) {
			tmp.sql <- paste("CREATE TABLE ", table.name, " (", paste(field.types, collapse=", "),
				", PRIMARY KEY (", paste(primary.key, collapse=", "), "))", sep="")
			return(tmp.sql)
		}

		"%w/o%" <- function(x, y) x[!x %in% y]

		convert.variables <- function(tmp.df, factor.variables=NULL) {
			if (length(grep("_", tmp.df$YEAR)) > 0) {
				tmp.df$YEAR <- sapply(strsplit(tmp.df$YEAR, "_"), '[', 2)
			}
			if (is.character(tmp.df$CONTENT_AREA)) {
				tmp.df$CONTENT_AREA <- as.factor(tmp.df$CONTENT_AREA)
			}
			tmp.factor.names <- c(factor.variables, names(tmp.df)[sapply(tmp.df, class)=="factor"] %w/o% c(group.number[2], group.number[1], "INSTRUCTOR_NUMBER"))
			for (i in tmp.factor.names) {
				tmp.df[[i]] <- unclass(as.factor(tmp.df[[i]]))
			}
			tmp.df[sapply(tmp.df, is.nan)] <- NA
			return(tmp.df)
		}

		get.grade <- function(grade) {
			if (SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Test_Season"]]=="Fall") grade-1 else grade
		}

		get.year <- function(year) {
			if (SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Test_Season"]]=="Fall") {
				yearIncrement(year, -1)
			} else {
				return(year)
			}
		}

		convert.names <- function(my.data) {
			names(my.data)[names(my.data)=="PERCENT_CATCHING_UP_KEEPING_UP"] <- "PERCENT_AT_ABOVE_TARGET"
			names(my.data)[names(my.data)==paste("MEDIAN_SGP_TARGET", projection.years.for.target, "YEAR", sep="_")] <- "MEDIAN_SGP_TARGET"
			if ("EMH_LEVEL" %in% names(my.data) && is.numeric(my.data[['EMH_LEVEL']])) {
				my.data[['EMH_LEVEL']] <- as.character(factor(my.data[['EMH_LEVEL']], levels=1:3, labels=c("E", "H", "M")))
			}
			if ("EMH_LEVEL" %in% names(my.data) && is.character(my.data[['EMH_LEVEL']])) {
				my.data[['EMH_LEVEL']] <- substr(my.data[['EMH_LEVEL']],1,1)
			}
			if ("GENDER" %in% names(my.data)) {
				my.data[['STUDENTGROUP']][my.data[['STUDENTGROUP']]=="Female"] <- "F"
				my.data[['STUDENTGROUP']][my.data[['STUDENTGROUP']]=="Male"] <- "M"
			}
			if ("INSTRUCTOR_NUMBER" %in% names(my.data)) names(my.data)[names(my.data)=="INSTRUCTOR_NUMBER"] <- "TEACHER_USID"
			if (group.number[1]!="DISTRICT") names(my.data)[names(my.data)==group.number[1]] <- "DISTRICT_NUMBER"
			if (group.number[2]!="SCHOOL") names(my.data)[names(my.data)==group.number[2]] <- "SCHOOL_NUMBER"
			return(my.data)
		}

	## Create relevant variables

		if (is.null(years)) years <- unique(sgp_object@Data$YEAR) %w/o% NA
		if (is.null(content_areas)) content_areas <- unique(sgp_object@Data$CONTENT_AREA) %w/o% NA
		if (!is.null(SGP::SGPstateData[[state]][["SGP_Configuration"]][["null.output.string"]])) {
			my.null.string <- SGP::SGPstateData[[state]][["SGP_Configuration"]][["null.output.string"]]
		} else {
			my.null.string <- "NULL"
		}

	## Create tmp.school.and.district.by.year table

		setkeyv(sgp_object@Data, c("YEAR", group.number[1], group.number[2]))
		tmp.school.and.district.by.year  <- as.data.frame(convert.variables(unique(sgp_object@Data)[, c("YEAR", group.number[1], group.number[2]), with=FALSE]))


	###
	### Construct tables
	###

	### Table 1. DISTRICT

		field.types <- c(
			"DISTRICT_NUMBER TEXT NOT NULL",
			"CONTENT_AREA TEXT NOT NULL",
			"YEAR INTEGER NOT NULL",
			"MEDIAN_SGP REAL",
			"MEDIAN_SGP_TARGET REAL",
			"PERCENT_AT_ABOVE_TARGET REAL",
			"PERCENT_AT_ABOVE_PROFICIENT REAL",
			"MEDIAN_SGP_COUNT INTEGER",
			"PERCENT_AT_ABOVE_PROFICIENT_COUNT INTEGER")

		tmp <- as.data.frame(convert.variables(subset(sgp_object@Summary[[group.number[1]]][[paste(group.number[1], "CONTENT_AREA__YEAR", group.enroll.status[1], sep="__")]],
			!is.na(get(group.number[1])) & CONTENT_AREA %in% content_areas & YEAR %in% years & get(group.enroll.status[1])==group.enroll.status.label[1] & !is.na(MEDIAN_SGP))))
		tmp <- convert.names(tmp)
		tmp <- tmp[, sapply(strsplit(field.types, " "), function(x) head(x,1))]

		dbGetQuery(db, sqlite.create.table("DISTRICT", field.types, c("YEAR", "DISTRICT_NUMBER", "CONTENT_AREA")))
		dbWriteTable(db, "DISTRICT", tmp, row.names=FALSE, append=TRUE)

		if (text.output) write.table(tmp, file=file.path(text.output.directory, "DISTRICT.dat"), row.names=FALSE, na=my.null.string, quote=FALSE, sep="|")
		if (json.output) cat(toJSON(tmp), file=file.path(json.output.directory, "DISTRICT.json"))


	### Table 2. DISTRICT_GRADE

		field.types <- c(
			"DISTRICT_NUMBER TEXT NOT NULL",
			"CONTENT_AREA TEXT NOT NULL",
			"YEAR INTEGER NOT NULL",
			"GRADE INTEGER NOT NULL",
			"MEDIAN_SGP REAL",
			"MEDIAN_SGP_TARGET REAL",
			"PERCENT_AT_ABOVE_TARGET REAL",
			"PERCENT_AT_ABOVE_PROFICIENT REAL",
			"MEDIAN_SGP_COUNT INTEGER",
			"PERCENT_AT_ABOVE_PROFICIENT_COUNT INTEGER")

		tmp <- as.data.frame(convert.variables(subset(sgp_object@Summary[[group.number[1]]][[paste(group.number[1], "CONTENT_AREA__YEAR__GRADE", group.enroll.status[1], sep="__")]],
			!is.na(get(group.number[1])) & CONTENT_AREA %in% content_areas & YEAR %in% years & !is.na(GRADE) & get(group.enroll.status[1])==group.enroll.status.label[1] & !is.na(MEDIAN_SGP))))
		tmp <- convert.names(tmp)
		tmp <- tmp[, sapply(strsplit(field.types, " "), function(x) head(x,1))]

		dbGetQuery(db, sqlite.create.table("DISTRICT_GRADE", field.types, c("YEAR", "DISTRICT_NUMBER", "CONTENT_AREA", "GRADE")))
		dbWriteTable(db, "DISTRICT_GRADE", tmp[, sapply(strsplit(field.types, " "), function(x) head(x,1))], row.names=FALSE, append=TRUE)

		if (text.output) write.table(tmp, file=file.path(text.output.directory, "DISTRICT_GRADE.dat"), row.names=FALSE, na=my.null.string, quote=FALSE, sep="|")
		if (json.output) cat(toJSON(tmp), file=file.path(json.output.directory, "DISTRICT_GRADE.json"))


	### Table 3. DISTRICT_ETHNICITY

		field.types <- c(
			"DISTRICT_NUMBER TEXT NOT NULL",
			"CONTENT_AREA TEXT NOT NULL",
			"YEAR INTEGER NOT NULL",
			"ETHNICITY INTEGER NOT NULL",
			"MEDIAN_SGP REAL",
			"MEDIAN_SGP_TARGET REAL",
			"PERCENT_AT_ABOVE_TARGET REAL",
			"PERCENT_AT_ABOVE_PROFICIENT REAL",
			"MEDIAN_SGP_COUNT INTEGER",
			"PERCENT_AT_ABOVE_PROFICIENT_COUNT INTEGER",
			"ENROLLMENT_PERCENTAGE REAL")

		tmp <- as.data.frame(convert.variables(subset(sgp_object@Summary[[group.number[1]]][[paste(group.number[1], "CONTENT_AREA__YEAR__ETHNICITY", group.enroll.status[1], sep="__")]],
			!is.na(get(group.number[1])) & CONTENT_AREA %in% content_areas & YEAR %in% years & !is.na(ETHNICITY) & get(group.enroll.status[1])==group.enroll.status.label[1] &
			!is.na(MEDIAN_SGP)), factor.variables="ETHNICITY"))
		tmp <- convert.names(tmp)
		tmp$ENROLLMENT_PERCENTAGE <- NA
		tmp <- tmp[, sapply(strsplit(field.types, " "), function(x) head(x,1))]

		dbGetQuery(db, sqlite.create.table("DISTRICT_ETHNICITY", field.types, c("YEAR", "DISTRICT_NUMBER", "CONTENT_AREA", "ETHNICITY")))
		dbWriteTable(db, "DISTRICT_ETHNICITY", tmp, row.names=FALSE, append=TRUE)

		if (text.output) write.table(tmp, file=file.path(text.output.directory, "DISTRICT_ETHNICITY.dat"), row.names=FALSE, na=my.null.string, quote=FALSE, sep="|")
		if (json.output) cat(toJSON(tmp), file=file.path(json.output.directory, "DISTRICT_ETHNICITY.json"))


	### Table 4. DISTRICT_GRADE_ETHNICITY

		field.types <- c(
			"DISTRICT_NUMBER TEXT NOT NULL",
			"CONTENT_AREA TEXT NOT NULL",
			"YEAR INTEGER NOT NULL",
			"GRADE TEXT NOT NULL",
			"ETHNICITY TEXT NOT NULL",
			"MEDIAN_SGP REAL",
			"MEDIAN_SGP_TARGET REAL",
			"PERCENT_AT_ABOVE_TARGET REAL",
			"PERCENT_AT_ABOVE_PROFICIENT REAL",
			"MEDIAN_SGP_COUNT INTEGER",
			"PERCENT_AT_ABOVE_PROFICIENT_COUNT INTEGER")

		tmp <- as.data.frame(convert.variables(subset(sgp_object@Summary[[group.number[1]]][[paste(group.number[1], "CONTENT_AREA__YEAR__GRADE__ETHNICITY", group.enroll.status[1], sep="__")]],
			!is.na(get(group.number[1])) & CONTENT_AREA %in% content_areas & YEAR %in% years & !is.na(GRADE) & !is.na(ETHNICITY) &
			get(group.enroll.status[1])==group.enroll.status.label[1] & !is.na(MEDIAN_SGP)), factor.variables="ETHNICITY"))
		tmp <- convert.names(tmp)
		tmp <- tmp[, sapply(strsplit(field.types, " "), function(x) head(x,1))]

		dbGetQuery(db, sqlite.create.table("DISTRICT_GRADE_ETHNICITY", field.types, c("YEAR", "DISTRICT_NUMBER", "CONTENT_AREA", "GRADE", "ETHNICITY")))
		dbWriteTable(db, "DISTRICT_GRADE_ETHNICITY", tmp, row.names=FALSE, append=TRUE)

		if (text.output) write.table(tmp, file=file.path(text.output.directory, "DISTRICT_GRADE_ETHNICITY.dat"), row.names=FALSE, na=my.null.string, quote=FALSE, sep="|")
		if (json.output) cat(toJSON(tmp), file=file.path(json.output.directory, "DISTRICT_GRADE_ETHNICITY.json"))


	### Table 5. DISTRICT_STUDENTGROUP

		field.types <- c(
			"DISTRICT_NUMBER TEXT NOT NULL",
			"CONTENT_AREA TEXT NOT NULL",
			"YEAR INTEGER NOT NULL",
			"STUDENTGROUP TEXT NOT NULL",
			"MEDIAN_SGP REAL",
			"MEDIAN_SGP_TARGET REAL",
			"PERCENT_AT_ABOVE_TARGET REAL",
			"PERCENT_AT_ABOVE_PROFICIENT REAL",
			"MEDIAN_SGP_COUNT INTEGER",
			"PERCENT_AT_ABOVE_PROFICIENT_COUNT INTEGER",
			"ENROLLMENT_PERCENTAGE REAL")

		tmp.list <- list()
		for (i in other.student.groups %w/o% grep("ETHNICITY", other.student.groups, value=TRUE)) {
			tmp.list[[i]] <- sgp_object@Summary[[group.number[1]]][[paste(group.number[1], "CONTENT_AREA__YEAR", i, group.enroll.status[1], sep="__")]]
		}

		for (i in seq_along(tmp.list)) {
			setnames(tmp.list[[i]], 4, "STUDENTGROUP")
		}

		tmp <- as.data.frame(convert.variables(subset(rbindlist(tmp.list, fill=TRUE),
			!is.na(get(group.number[1])) & CONTENT_AREA %in% content_areas & YEAR %in% years & !is.na(STUDENTGROUP) & get(group.enroll.status[1])==group.enroll.status.label[1] &
			!is.na(MEDIAN_SGP)), factor.variables="STUDENTGROUP"))

		tmp <- convert.names(tmp)
		tmp$ENROLLMENT_PERCENTAGE <- NA
		tmp <- data.table(tmp[, sapply(strsplit(field.types, " "), function(x) head(x,1))], key=c("YEAR", "DISTRICT_NUMBER", "CONTENT_AREA", "STUDENTGROUP"))
		tmp <- as.data.frame(data.table(tmp[!duplicated(tmp)]))

		dbGetQuery(db, sqlite.create.table("DISTRICT_STUDENTGROUP", field.types, c("YEAR", "DISTRICT_NUMBER", "CONTENT_AREA", "STUDENTGROUP")))
		dbWriteTable(db, "DISTRICT_STUDENTGROUP", tmp, row.names=FALSE, append=TRUE)

		if (text.output) write.table(tmp, file=file.path(text.output.directory, "DISTRICT_STUDENTGROUP.dat"), row.names=FALSE, na=my.null.string, quote=FALSE, sep="|")
		if (json.output) cat(toJSON(tmp), file=file.path(json.output.directory, "DISTRICT_STUDENTGROUP.json"))


	### Table 6. DISTRICT_GRADE_STUDENTGROUP

		field.types <- c(
			"DISTRICT_NUMBER TEXT NOT NULL",
			"CONTENT_AREA TEXT NOT NULL",
			"YEAR INTEGER NOT NULL",
			"GRADE TEXT NOT NULL",
			"STUDENTGROUP TEXT NOT NULL",
			"MEDIAN_SGP REAL",
			"MEDIAN_SGP_TARGET REAL",
			"PERCENT_AT_ABOVE_TARGET REAL",
			"PERCENT_AT_ABOVE_PROFICIENT REAL",
			"MEDIAN_SGP_COUNT INTEGER",
			"PERCENT_AT_ABOVE_PROFICIENT_COUNT INTEGER")

		tmp.list <- list()
		for (i in other.student.groups %w/o% grep("ETHNICITY", other.student.groups, value=TRUE)) {
			tmp.list[[i]] <- sgp_object@Summary[[group.number[1]]][[paste(group.number[1], "CONTENT_AREA__YEAR__GRADE", i, group.enroll.status[1], sep="__")]]
		}

		for (i in seq_along(tmp.list)) {
			setnames(tmp.list[[i]], 5, "STUDENTGROUP")
		}

		tmp <- as.data.frame(convert.variables(subset(rbindlist(tmp.list, fill=TRUE),
			!is.na(get(group.number[1])) & YEAR %in% years & CONTENT_AREA %in% content_areas & !is.na(STUDENTGROUP) & get(group.enroll.status[1])==group.enroll.status.label[1] &
			!is.na(MEDIAN_SGP)), factor.variables="STUDENTGROUP"))
		tmp <- convert.names(tmp)
                tmp <- data.table(tmp[, sapply(strsplit(field.types, " "), function(x) head(x,1))], key=c("YEAR", "DISTRICT_NUMBER", "CONTENT_AREA", "GRADE", "STUDENTGROUP"))
                tmp <- as.data.frame(tmp[!duplicated(tmp)])

		dbGetQuery(db, sqlite.create.table("DISTRICT_GRADE_STUDENTGROUP", field.types, c("YEAR", "DISTRICT_NUMBER", "CONTENT_AREA", "GRADE", "STUDENTGROUP")))
		dbWriteTable(db, "DISTRICT_GRADE_STUDENTGROUP", tmp, row.names=FALSE, append=TRUE)

		if (text.output) write.table(tmp, file=file.path(text.output.directory, "DISTRICT_GRADE_STUDENTGROUP.dat"), row.names=FALSE, na=my.null.string, quote=FALSE, sep="|")
		if (json.output) cat(toJSON(tmp), file=file.path(json.output.directory, "DISTRICT_GRADE_STUDENTGROUP.json"))


	## Table 7. SCHOOL

		field.types <- c(
			"DISTRICT_NUMBER TEXT NOT NULL",
			"SCHOOL_NUMBER TEXT NOT NULL",
			"EMH_LEVEL TEXT NOT NULL",
			"CONTENT_AREA TEXT NOT NULL",
			"YEAR INTEGER NOT NULL",
			"MEDIAN_SGP REAL",
			"MEDIAN_SGP_TARGET REAL",
			"PERCENT_AT_ABOVE_TARGET REAL",
			"PERCENT_AT_ABOVE_PROFICIENT REAL",
			"MEDIAN_SGP_COUNT INTEGER",
			"PERCENT_AT_ABOVE_PROFICIENT_COUNT INTEGER")

		tmp <- as.data.frame(convert.variables(subset(sgp_object@Summary[[group.number[2]]][[paste(group.number[2], "EMH_LEVEL__CONTENT_AREA__YEAR", group.enroll.status[2], sep="__")]],
			!is.na(get(group.enroll.status[2])) & !is.na(get(group.number[2])) & !is.na(EMH_LEVEL) & CONTENT_AREA %in% content_areas & YEAR %in% years & get(group.enroll.status[2])==group.enroll.status.label[2] &
			!is.na(MEDIAN_SGP))))
		tmp <- as.data.frame(merge(tmp, as.data.frame(tmp.school.and.district.by.year), all.x=TRUE))
		tmp <- convert.names(tmp)
		tmp <- tmp[, sapply(strsplit(field.types, " "), function(x) head(x,1))]

		dbGetQuery(db, sqlite.create.table("SCHOOL", field.types, c("YEAR", "DISTRICT_NUMBER", "SCHOOL_NUMBER", "EMH_LEVEL", "CONTENT_AREA")))
		dbWriteTable(db, "SCHOOL", tmp, row.names=FALSE, append=TRUE)

		if (text.output) write.table(tmp, file=file.path(text.output.directory, "SCHOOL.dat"), row.names=FALSE, na=my.null.string, quote=FALSE, sep="|")
		if (json.output) cat(toJSON(tmp), file=file.path(json.output.directory, "SCHOOL.json"))


	## Table 8. SCHOOL_GRADE

		field.types <- c(
			"DISTRICT_NUMBER TEXT NOT NULL",
			"SCHOOL_NUMBER TEXT NOT NULL",
			"EMH_LEVEL TEXT NOT NULL",
			"CONTENT_AREA TEXT NOT NULL",
			"YEAR INTEGER NOT NULL",
			"GRADE TEXT NOT NULL",
			"MEDIAN_SGP REAL",
			"MEDIAN_SGP_TARGET REAL",
			"PERCENT_AT_ABOVE_TARGET REAL",
			"PERCENT_AT_ABOVE_PROFICIENT REAL",
			"MEDIAN_SGP_COUNT INTEGER",
			"PERCENT_AT_ABOVE_PROFICIENT_COUNT INTEGER")

		tmp <- as.data.frame(convert.variables(subset(sgp_object@Summary[[group.number[2]]][[paste(group.number[2], "EMH_LEVEL__CONTENT_AREA__YEAR__GRADE", group.enroll.status[2], sep="__")]],
			!is.na(get(group.number[2])) & !is.na(EMH_LEVEL) & CONTENT_AREA %in% content_areas & YEAR %in% years & !is.na(GRADE) & get(group.enroll.status[2])==group.enroll.status.label[2] &
			!is.na(MEDIAN_SGP))))
		tmp <- data.frame(merge(tmp, as.data.frame(tmp.school.and.district.by.year), all.x=TRUE))
		tmp <- convert.names(tmp)
		tmp <- tmp[, sapply(strsplit(field.types, " "), function(x) head(x,1))]

		dbGetQuery(db, sqlite.create.table("SCHOOL_GRADE", field.types, c("YEAR", "DISTRICT_NUMBER", "SCHOOL_NUMBER", "EMH_LEVEL", "GRADE", "CONTENT_AREA")))
		dbWriteTable(db, "SCHOOL_GRADE", tmp, row.names=FALSE, append=TRUE)

		if (text.output) write.table(tmp, file=file.path(text.output.directory, "SCHOOL_GRADE.dat"), row.names=FALSE, na=my.null.string, quote=FALSE, sep="|")
		if (json.output) cat(toJSON(tmp), file=file.path(json.output.directory, "SCHOOL_GRADE.json"))


	## Table 9. SCHOOL_ETHNICITY

		field.types <- c(
			"DISTRICT_NUMBER TEXT NOT NULL",
			"SCHOOL_NUMBER TEXT NOT NULL",
			"EMH_LEVEL TEXT NOT NULL",
			"CONTENT_AREA TEXT NOT NULL",
			"YEAR INTEGER NOT NULL",
			"ETHNICITY TEXT NOT NULL",
			"MEDIAN_SGP REAL",
			"MEDIAN_SGP_TARGET REAL",
			"PERCENT_AT_ABOVE_TARGET REAL",
			"PERCENT_AT_ABOVE_PROFICIENT REAL",
			"MEDIAN_SGP_COUNT INTEGER",
			"PERCENT_AT_ABOVE_PROFICIENT_COUNT INTEGER",
			"ENROLLMENT_PERCENTAGE REAL")

		tmp <- as.data.frame(convert.variables(subset(sgp_object@Summary[[group.number[2]]][[paste(group.number[2], "EMH_LEVEL__CONTENT_AREA__YEAR__ETHNICITY", group.enroll.status[2], sep="__")]],
			!is.na(get(group.number[2])) & !is.na(EMH_LEVEL) & CONTENT_AREA %in% content_areas & YEAR %in% years & !is.na(ETHNICITY) & get(group.enroll.status[2])==group.enroll.status.label[2] &
			!is.na(MEDIAN_SGP)), factor.variables="ETHNICITY"))
		tmp <- data.frame(merge(tmp, as.data.frame(tmp.school.and.district.by.year), all.x=TRUE))
		tmp <- convert.names(tmp)
		tmp$ENROLLMENT_PERCENTAGE <- NA
		tmp <- tmp[, sapply(strsplit(field.types, " "), function(x) head(x,1))]

		dbGetQuery(db, sqlite.create.table("SCHOOL_ETHNICITY", field.types, c("YEAR", "DISTRICT_NUMBER", "SCHOOL_NUMBER", "EMH_LEVEL", "CONTENT_AREA", "ETHNICITY")))
		dbWriteTable(db, "SCHOOL_ETHNICITY", tmp, row.names=FALSE, append=TRUE)

		if (text.output) write.table(tmp, file=file.path(text.output.directory, "SCHOOL_ETHNICITY.dat"), row.names=FALSE, na=my.null.string, quote=FALSE, sep="|")
		if (json.output) cat(toJSON(tmp), file=file.path(json.output.directory, "SCHOOL_ETHNICITY.json"))


	## Table 10. SCHOOL_STUDENTGROUP

		field.types <- c(
			"DISTRICT_NUMBER TEXT NOT NULL",
			"SCHOOL_NUMBER TEXT NOT NULL",
			"EMH_LEVEL TEXT NOT NULL",
			"CONTENT_AREA TEXT NOT NULL",
			"YEAR INTEGER NOT NULL",
			"STUDENTGROUP TEXT NOT NULL",
			"MEDIAN_SGP REAL",
			"MEDIAN_SGP_TARGET REAL",
			"PERCENT_AT_ABOVE_TARGET REAL",
			"PERCENT_AT_ABOVE_PROFICIENT REAL",
			"MEDIAN_SGP_COUNT INTEGER",
			"PERCENT_AT_ABOVE_PROFICIENT_COUNT INTEGER",
			"ENROLLMENT_PERCENTAGE REAL")

		tmp.list <- list()
		for (i in other.student.groups %w/o% grep("ETHNICITY", other.student.groups, value=TRUE)) {
			tmp.list[[i]] <- sgp_object@Summary[[group.number[2]]][[paste(group.number[2], "EMH_LEVEL__CONTENT_AREA__YEAR", i, group.enroll.status[2], sep="__")]]
		}

		for (i in seq_along(tmp.list)) {
			setnames(tmp.list[[i]], 5, "STUDENTGROUP")
		}

		tmp <- as.data.frame(convert.variables(subset(rbindlist(tmp.list, fill=TRUE),
			!is.na(get(group.number[2])) & !is.na(EMH_LEVEL) & CONTENT_AREA %in% content_areas & YEAR %in% years & !is.na(STUDENTGROUP) & get(group.enroll.status[2])==group.enroll.status.label[2] &
			!is.na(MEDIAN_SGP)), factor.variables="STUDENTGROUP"))
		tmp <- as.data.frame(merge(tmp, as.data.frame(tmp.school.and.district.by.year), all.x=TRUE))
		tmp <- convert.names(tmp)
		tmp$ENROLLMENT_PERCENTAGE <- NA
                tmp <- data.table(tmp[, sapply(strsplit(field.types, " "), function(x) head(x,1))], key=c("YEAR", "DISTRICT_NUMBER", "SCHOOL_NUMBER", "EMH_LEVEL", "CONTENT_AREA", "STUDENTGROUP"))
                tmp <- as.data.frame(tmp[!duplicated(tmp)])

		dbGetQuery(db, sqlite.create.table("SCHOOL_STUDENTGROUP", field.types,
			c("YEAR", "DISTRICT_NUMBER", "SCHOOL_NUMBER", "EMH_LEVEL", "CONTENT_AREA", "STUDENTGROUP")))
		dbWriteTable(db, "SCHOOL_STUDENTGROUP", tmp, row.names=FALSE, append=TRUE)

		if (text.output) write.table(tmp, file=file.path(text.output.directory, "SCHOOL_STUDENTGROUP.dat"), row.names=FALSE, na=my.null.string, quote=FALSE, sep="|")
		if (json.output) cat(toJSON(tmp), file=file.path(json.output.directory, "SCHOOL_STUDENTGROUP.json"))


	## Table 11. SCHOOL_TEACHER

	if (any(c(paste(group.number[2], "INSTRUCTOR_NUMBER__EMH_LEVEL__CONTENT_AREA__YEAR", sep="__"),
		paste(group.number[2], "INSTRUCTOR_NUMBER__EMH_LEVEL__CONTENT_AREA__YEAR__INSTRUCTOR_ENROLLMENT_STATUS", sep="__")) %in% names(sgp_object@Summary[[group.number[2]]]))) {

		field.types <- c(
			"DISTRICT_NUMBER TEXT NOT NULL",
			"SCHOOL_NUMBER TEXT NOT NULL",
			"EMH_LEVEL TEXT NOT NULL",
			"TEACHER_USID TEXT NOT NULL",
			"CONTENT_AREA TEXT NOT NULL",
			"YEAR INTEGER NOT NULL",
			"MEDIAN_SGP REAL",
			"MEDIAN_SGP_TARGET REAL",
			"PERCENT_AT_ABOVE_TARGET REAL",
			"PERCENT_AT_ABOVE_PROFICIENT REAL",
			"MEDIAN_SGP_COUNT INTEGER",
			"PERCENT_AT_ABOVE_PROFICIENT_COUNT INTEGER")

		if (paste(group.number[2], "INSTRUCTOR_NUMBER__EMH_LEVEL__CONTENT_AREA__YEAR__INSTRUCTOR_ENROLLMENT_STATUS", sep="__") %in% names(sgp_object@Summary[[group.number[2]]])) {
			tmp.table.name <- paste(group.number[2], "INSTRUCTOR_NUMBER__EMH_LEVEL__CONTENT_AREA__YEAR__INSTRUCTOR_ENROLLMENT_STATUS", sep="__")
			tmp <- as.data.frame(convert.variables(subset(sgp_object@Summary[[group.number[2]]][[tmp.table.name]],
				!is.na(get(group.number[2])) & !is.na(INSTRUCTOR_NUMBER) & !is.na(EMH_LEVEL) & CONTENT_AREA %in% content_areas & YEAR %in% years &
				INSTRUCTOR_ENROLLMENT_STATUS=="Enrolled Instructor: Yes" & !is.na(MEDIAN_SGP))))
		} else {
			tmp.table.name <- paste(group.number[2], "INSTRUCTOR_NUMBER__EMH_LEVEL__CONTENT_AREA__YEAR", sep="__")
			tmp <- as.data.frame(convert.variables(subset(sgp_object@Summary[[group.number[2]]][[tmp.table.name]],
				!is.na(get(group.number[2])) & !is.na(INSTRUCTOR_NUMBER) & !is.na(EMH_LEVEL) & CONTENT_AREA %in% content_areas & YEAR %in% years &
				!is.na(MEDIAN_SGP))))
		}

		tmp <- data.frame(merge(tmp, as.data.frame(tmp.school.and.district.by.year), all.x=TRUE))
		tmp <- convert.names(tmp)
		tmp <- tmp[, sapply(strsplit(field.types, " "), function(x) head(x,1))]

		dbGetQuery(db, sqlite.create.table("SCHOOL_TEACHER", field.types, c("YEAR", "DISTRICT_NUMBER", "SCHOOL_NUMBER", "TEACHER_USID", "EMH_LEVEL", "CONTENT_AREA")))
		dbWriteTable(db, "SCHOOL_TEACHER", tmp, row.names=FALSE, append=TRUE)

		if (text.output) write.table(tmp, file=file.path(text.output.directory, "SCHOOL_TEACHER.dat"), row.names=FALSE, na=my.null.string, quote=FALSE, sep="|")
		if (json.output) cat(toJSON(tmp), file=file.path(json.output.directory, "SCHOOL_TEACHER.json"))

	} ### END SCHOOL_TEACHER table


	## Table 12. KEY_VALUE_LOOKUP (ADD CONTENT_AREA and YEAR)

		field.types <- c(
			"KEY_VALUE_ID INTEGER NOT NULL",
			"KEY_VALUE_KEY TEXT",
			"KEY_VALUE_CODE TEXT",
			"KEY_VALUE_TEXT TEXT")

		# CONTENT_AREA

		tmp <- subset(sgp_object@Summary[[group.number[1]]][[paste(group.number[1], "CONTENT_AREA__YEAR", group.enroll.status[1], sep="__")]],
			!is.na(get(group.number[1])) & CONTENT_AREA %in% content_areas & YEAR %in% years & get(group.enroll.status[1])==group.enroll.status.label[1])
		tmp.CONTENT_AREA <- data.frame(
			KEY_VALUE_KEY="CONTENT_AREA",
			KEY_VALUE_CODE=seq_along(unique(tmp$CONTENT_AREA)),
			KEY_VALUE_TEXT=sapply(sort(unique(tmp$CONTENT_AREA)), capwords))

		# YEAR

		tmp <- convert.variables(subset(sgp_object@Summary[[group.number[1]]][[paste(group.number[1], "CONTENT_AREA__YEAR", group.enroll.status[1], sep="__")]],
			!is.na(get(group.number[1])) & CONTENT_AREA %in% content_areas & YEAR %in% years & get(group.enroll.status[1])==group.enroll.status.label[1]))
		tmp.YEAR <- data.frame(
			KEY_VALUE_KEY="YEAR",
			KEY_VALUE_CODE=sort(unique(tmp$YEAR)),
			KEY_VALUE_TEXT=paste(as.numeric(sapply(sort(unique(tmp$YEAR)), get.year))-1, "-", sapply(sort(unique(tmp$YEAR)), get.year), sep=""))

		# GRADE

		tmp <- subset(as.data.frame(sgp_object@Summary[[group.number[1]]][[paste(group.number[1], "CONTENT_AREA__YEAR__GRADE", group.enroll.status[1], sep="__")]]),
			!is.na(get(group.number[1])) & CONTENT_AREA %in% content_areas & YEAR %in% years & !is.na(GRADE) & get(group.enroll.status[1])==group.enroll.status.label[1])
		tmp.GRADE <- data.frame(
			KEY_VALUE_KEY="GRADE",
			KEY_VALUE_CODE=sort(unique(as.integer(tmp$GRADE))),
			KEY_VALUE_TEXT=paste("Grade", get.grade(sort(unique(as.integer(tmp$GRADE))))))

		# EMH_LEVEL

		tmp <- subset(sgp_object@Summary[[group.number[2]]][[paste(group.number[2], "EMH_LEVEL__CONTENT_AREA__YEAR", group.enroll.status[2], sep="__")]],
			!is.na(get(group.number[2])) & !is.na(EMH_LEVEL) & CONTENT_AREA %in% content_areas & YEAR %in% years & get(group.enroll.status[2])==group.enroll.status.label[2])
		if (!is.factor(tmp$EMH_LEVEL)) tmp[['EMH_LEVEL']] <- as.factor(tmp[['EMH_LEVEL']])
		tmp.EMH <- data.frame(
			KEY_VALUE_KEY="EMH_LEVEL",
			KEY_VALUE_CODE=strhead(levels(as.factor(tmp$EMH_LEVEL))[sort(unique(as.integer(as.factor(tmp$EMH_LEVEL))))], 1), ## TEMP fix until EMH_LEVEL is fixed up
			KEY_VALUE_TEXT= levels(as.factor(tmp$EMH_LEVEL))[sort(unique(as.integer(as.factor(tmp$EMH_LEVEL))))])

		# ETHNICITY

		tmp <- subset(as.data.frame(sgp_object@Summary[[group.number[1]]][[paste(group.number[1], "CONTENT_AREA__YEAR__ETHNICITY", group.enroll.status[1], sep="__")]]),
			!is.na(get(group.number[1])) & CONTENT_AREA %in% content_areas & YEAR %in% years & !is.na(ETHNICITY) & get(group.enroll.status[1])==group.enroll.status.label[1])
		tmp.ETHNICITY <- data.frame(
			KEY_VALUE_KEY="ETHNICITY",
			KEY_VALUE_CODE=sort(unique(as.integer(as.factor(tmp$ETHNICITY)))),
			KEY_VALUE_TEXT=levels(as.factor(tmp$ETHNICITY))[sort(unique(as.integer(as.factor(tmp$ETHNICITY))))])

		# STUDENTGROUP

		tmp.list <- list()
		for (i in other.student.groups %w/o% grep("ETHNICITY", other.student.groups, value=TRUE)) {
			tmp.list[[i]] <- sgp_object@Summary[[group.number[1]]][[paste(group.number[1], "CONTENT_AREA__YEAR", i, group.enroll.status[1], sep="__")]]
		}

		for (i in seq_along(tmp.list)) {
			setnames(tmp.list[[i]], 4, "STUDENTGROUP")
		}

		tmp <- data.table(convert.names(convert.variables(subset(rbindlist(tmp.list, fill=TRUE),
			!is.na(get(group.number[1])) & !is.na(STUDENTGROUP) & get(group.enroll.status[1])==group.enroll.status.label[1]))),
			key=c("YEAR", "DISTRICT_NUMBER", "CONTENT_AREA", "STUDENTGROUP"))
                tmp <- as.data.frame(data.table(tmp[!duplicated(tmp)]))

		tmp.STUDENTGROUP <- data.frame(
			KEY_VALUE_KEY="STUDENT_GROUP", ### NOTE: Must have underscore. It's an older version of the table
			KEY_VALUE_CODE=sort(unique(as.integer(as.factor(tmp$STUDENTGROUP)))),
			KEY_VALUE_TEXT=levels(as.factor(tmp$STUDENTGROUP))[sort(unique(as.integer(as.factor(tmp$STUDENTGROUP))))])


		tmp <- rbind(tmp.CONTENT_AREA, tmp.YEAR, tmp.GRADE, tmp.EMH, tmp.ETHNICITY, tmp.STUDENTGROUP)
		tmp <- data.frame(KEY_VALUE_ID=1:dim(tmp)[1], tmp)

		dbGetQuery(db, sqlite.create.table("KEY_VALUE_LOOKUP", field.types, "KEY_VALUE_ID"))
		dbWriteTable(db, "KEY_VALUE_LOOKUP", tmp, row.names=FALSE, append=TRUE)

		if (text.output) write.table(tmp, file=file.path(text.output.directory, "KEY_VALUE_LOOKUP.dat"), row.names=FALSE, na=my.null.string, quote=FALSE, sep="|")
		if (json.output) cat(toJSON(tmp), file=file.path(json.output.directory, "KEY_VALUE_LOOKUP.json"))

###
### Disconnect database
###

	dbDisconnect(db)

	message(paste("\tFinished sqliteSGP in outputSGP", date(), "in", convertTime(timetaken(started.at)), "\n"))

} ### END sqliteSGP
