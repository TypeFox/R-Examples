`rliSGP` <-
function(sgp_object,
	additional.data=NULL,
	state=NULL,
	content_areas=c("MATHEMATICS", "READING", "EARLY_LITERACY"),
	testing.window=NULL, ### FALL, WINTER, SPRING
	eow.or.update="UPDATE", ### UPDATE or EOW
	update.save.shell.only=FALSE,
	configuration.year=NULL,
	sgp.percentiles.baseline=TRUE,
	sgp.projections.baseline=TRUE,
	sgp.projections.lagged.baseline=FALSE,
	sgp.target.scale.scores=TRUE,
	update.ids=NULL,
	SGPt=TRUE,
	save.intermediate.results=FALSE,
	coefficient.matrices=NULL,
	goodness.of.fit.print=FALSE,
	return.updated.shell=FALSE,
	fix.duplicates="KEEP.ALL",
	eow.calculate.sgps=FALSE,
	parallel.config=NULL) {

	YEAR <- GRADE <- ID <- NEW_ID <- .EACHI <- DATE <- NULL

	started.at <- proc.time()
	message(paste("\nStarted rliSGP", date()), "\n")

        if (is.null(state)) {
                tmp.name <- toupper(gsub("_", " ", deparse(substitute(sgp_object))))
                state <- getStateAbbreviation(tmp.name, "abcSGP")
        }

	if (!state %in% c("RLI", "RLI_UK")) stop("\tNOTE: 'rliSGP' only works with states RLI or RLI_UK currently")


	### Utility functions

	convertToBaseline <- function(baseline_matrices) {
		tmp.list <- list()
		if (is.null(baseline_matrices)) {
			return(NULL)
		} else {
			for (i in names(baseline_matrices)) {
				for (j in seq_along(baseline_matrices[[i]])) {
					baseline_matrices[[i]][[j]]@Time <- list(rep("BASELINE", length(unlist(baseline_matrices[[i]][[j]]@Time))))
				}
				names(baseline_matrices[[i]]) <- sub("[.][1234]_", "_", names(baseline_matrices[[i]]))
			}

			tmp.content_areas <- unique(sapply(strsplit(names(baseline_matrices), "[.]"), '[', 1))
			for (i in tmp.content_areas) {
				tmp.list[[paste(i, "BASELINE", sep=".")]] <- unlist(baseline_matrices[grep(i, names(baseline_matrices))], recursive=FALSE)
			}
			return(tmp.list)
		}
	}

	updateIDS <- function(my.data, id.lookup) {
		setnames(id.lookup, 1:2, c("ID", "NEW_ID"))
		id.lookup[,ID:=as.character(ID)]; id.lookup[,NEW_ID:=as.character(NEW_ID)]
		setkey(id.lookup, ID)
		if (is.SGP(my.data)) {
			tmp.dt <- copy(my.data@Data)
			setkey(tmp.dt, ID)
			tmp.dt[id.lookup, ID:=NEW_ID, by=.EACHI]
			sgp_object@Data <- tmp.dt
			setkeyv(sgp_object@Data, c("VALID_CASE", "CONTENT_AREA", "YEAR", "ID"))
			return(sgp_object)
		}
		if (is.data.frame(my.data)) {
			setkey(my.data, ID)
			my.data[id.lookup, ID:=NEW_ID, by=.EACHI]
			return(my.data)
		}
	}

	getRLIConfig <- function(content_areas, configuration.year, testing.window, SGPt) {
		tmp.list <- list()
		for (i in content_areas) {
			tmp.list[[i]] <- SGP::SGPstateData$RLI$SGP_Configuration$sgp.config.function$value(configuration.year, i, testing.window)
		}
		return(unlist(tmp.list, recursive=FALSE))
	}


	### Tests for arguments

	if (!is.null(additional.data) && !is.data.table(additional.data)) additional.data <- as.data.table(additional.data)

	if ("DATE" %in% names(additional.data)) additional.data[,DATE:=as.Date(DATE)]

	if (!is.null(update.ids) && !is.data.table(update.ids)) update.ids <- as.data.table(update.ids)

	if (state=="RLI_UK") content_areas <- "READING"

	if (long.data.supplied <- is.data.frame(sgp_object)) {
		tmp.last.year <- tail(sort(unique(sgp_object[['YEAR']])), 1)
		additional.data <- sgp_object[YEAR==tmp.last.year]
		sgp_object <- new("SGP", Data=suppressMessages(prepareSGP(sgp_object[YEAR!=tmp.last.year], state=state)@Data), Version=getVersion(sgp_object))
		gc(FALSE)
	}

	if (length(find.package("RLImatrices", quiet=TRUE))==0) stop("Package RLImatrices required from GitHub.")
	if (is.null(coefficient.matrices)) {
		eval(parse(text="require(RLImatrices)"))
		SGPstateData[[state]][["Baseline_splineMatrix"]][["Coefficient_Matrices"]] <-
		eval(parse(text=paste(paste(state, "SGPt_Baseline_Matrices", sep="_"), "$", paste(state, "SGPt_Baseline_Matrices", max(tail(sort(unique(additional.data[['YEAR']])), 1), "2013_2014.1"), sep="_"), sep="")))
	} else {
		SGPstateData[[state]][["Baseline_splineMatrix"]][["Coefficient_Matrices"]] <- coefficient.matrices
	}

	if (!is.null(testing.window) && (length(testing.window) != 1 || !testing.window %in% c("FALL", "WINTER", "SPRING"))) {
		stop("\tPlease supply either 'FALL', 'WINTER', or 'SPRING' for the testing.window argument.")
	} else {
		testing.window <- c("FALL", "WINTER", "SPRING")[as.numeric(tail(unlist(strsplit(tail(sort(unique(additional.data[['YEAR']])), 1), '[.]')), 1))]
	}

	if (is.null(configuration.year)) configuration.year <- head(unlist(strsplit(tail(sort(unique(additional.data[['YEAR']])), 1), '[.]')), 1)

	### Create variables

	if (is.null(SGPt)) update.shell.name <- paste(state, "SGP_UPDATE_SHELL", sep="_") else update.shell.name <- paste(state, "SGPt_UPDATE_SHELL", sep="_")
	if (testing.window=="FALL") num.windows.to.keep <- 5 else num.windows.to.keep <- 6


	### Update IDS if requested

	if (!is.null(update.ids)) {
		sgp_object <- updateIDS(sgp_object, update.ids)
		additional.data <- updateIDS(additional.data, update.ids)
	}


	########################################################################
	###
	### WITHIN_WINDOW UPDATE scripts
	###
	########################################################################

	if (eow.or.update=="UPDATE") {

		sgp_object <- updateSGP(
			what_sgp_object=sgp_object,
			with_sgp_data_LONG=additional.data,
			state=state,
			steps=c("prepareSGP", "analyzeSGP", "combineSGP", "outputSGP"),
			save.intermediate.results=save.intermediate.results,
			sgp.percentiles=FALSE,
			sgp.projections=FALSE,
			sgp.projections.lagged=FALSE,
			sgp.percentiles.baseline=sgp.percentiles.baseline,
			sgp.projections.baseline=sgp.projections.baseline,
			sgp.projections.lagged.baseline=sgp.projections.lagged.baseline,
			sgp.target.scale.scores=sgp.target.scale.scores,
			sgp.target.scale.scores.only=TRUE,
			outputSGP.output.type="RLI",
			goodness.of.fit.print=goodness.of.fit.print,
			update.old.data.with.new=FALSE,
			SGPt=SGPt,
			fix.duplicates=fix.duplicates,
			parallel.config=parallel.config,
			sgp.config=getRLIConfig(content_areas, configuration.year, testing.window, SGPt))

		if (!is.null(update.ids)) {
			assign(update.shell.name, sgp_object)
			save(list=update.shell.name, paste(update.shell.name, "Rdata", sep="."))
		}

		if (update.save.shell.only) {
			assign(update.shell.name, prepareSGP(subset(sgp_object@Data, YEAR %in% tail(head(sort(unique(sgp_object@Data$YEAR)), -1), num.windows.to.keep)),
				state=state, create.additional.variables=FALSE))
			save(list=update.shell.name, file=paste(update.shell.name, "Rdata", sep="."))
		}
	} ### END UPDATE scripts


	###############################################################################
	###
	### END_OF_WINDOW UPDATE scripts
	###
	###############################################################################

	if (eow.or.update=="EOW") {

		if (update.save.shell.only) {
			tmp.data <- rbindlist(list(sgp_object@Data, additional.data), fill=TRUE)
			assign(update.shell.name, prepareSGP(subset(tmp.data, YEAR %in% tail(sort(unique(tmp.data$YEAR)), num.windows.to.keep)), state=state, create.additional.variables=FALSE))
			save(list=update.shell.name, file=paste(update.shell.name, "Rdata", sep="."))
		} else {
			if (eow.calculate.sgps) my.steps <- c("prepareSGP", "analyzeSGP", "combineSGP", "outputSGP") else steps <- c("prepareSGP", "analyzeSGP")
			sgp_object <- updateSGP(
				what_sgp_object=sgp_object,
				with_sgp_data_LONG=additional.data,
				state=state,
				steps=steps,
				save.intermediate.results=save.intermediate.results,
				sgp.percentiles=TRUE,
				sgp.projections=FALSE,
				sgp.projections.lagged=FALSE,
				sgp.percentiles.baseline=sgp.percentiles.baseline & eow.calculate.sgps,
				sgp.projections.baseline=sgp.projections.baseline & eow.calculate.sgps,
				sgp.projections.lagged.baseline=sgp.projections.lagged.baseline & eow.calculate.sgps,
				sgp.target.scale.scores=sgp.target.scale.scores & eow.calculate.sgps,
				sgp.target.scale.scores.only=TRUE,
				outputSGP.output.type="RLI",
				update.old.data.with.new=TRUE,
				goodness.of.fit.print=goodness.of.fit.print,
				SGPt=SGPt,
				sgp.percentiles.calculate.sgps=eow.calculate.sgps,
				parallel.config=parallel.config,
				sgp.config=getRLIConfig(content_areas, configuration.year, testing.window, SGPt))

			### Create and save new UPDATE_SHELL

			if (!long.data.supplied) {
				assign(update.shell.name, prepareSGP(sgp_object@Data[YEAR %in% tail(sort(unique(sgp_object@Data$YEAR)), num.windows.to.keep)],
					state=state, create.additional.variables=FALSE))
				save(list=update.shell.name, file=paste(update.shell.name, "Rdata", sep="."))
			}


			### Convert and save coefficient matrices for inclusion in RLImatrices package

			if (testing.window=="FALL") {
				matrix.window <- paste(configuration.year, 3, sep=".")
			} else {
				matrix.window <- paste(yearIncrement(configuration.year, 1), c(3, 1, 2)[match(testing.window, c("FALL", "WINTER", "SPRING"))], sep=".")
			}
			new.matrices <-convertToBaseline(sgp_object@SGP$Coefficient_Matrices[grep(configuration.year, names(sgp_object@SGP$Coefficient_Matrices))])
			old.matrices <- SGPstateData[[state]][["Baseline_splineMatrix"]][["Coefficient_Matrices"]]
			year.to.replace <- head(sort(unique(sapply(lapply(sapply(names(old.matrices[['READING.BASELINE']]), strsplit, '[.]'), '[', 2:3), paste, collapse="."))), 1)
			for (content_area.iter in c("EARLY_LITERACY.BASELINE", "READING.BASELINE", "MATHEMATICS.BASELINE")) {
				old.matrices[[content_area.iter]][grep(year.to.replace, names(old.matrices[[content_area.iter]]))] <- NULL
				old.matrices[[content_area.iter]] <- c(old.matrices[[content_area.iter]], new.matrices[[content_area.iter]])
			}
			eval(parse(text=paste(paste(state, "SGPt_Baseline_Matrices$", sep="_"), paste(state, "SGPt_Baseline_Matrices", matrix.window, sep="_"), " <- old.matrices", sep="")))
			save(list=paste(state, "SGPt_Baseline_Matrices", sep="_"), file=paste(paste(state, "SGPt_Baseline_Matrices", sep="_"), "rda", sep="."), compress="xz")
			message(paste("\tNOTE: ", paste(state, "SGPt_Baseline_Matrices", sep="_"), " saved to working directory contains matrices for use in ", matrix.window, ".", sep=""))
			message(paste("\t\tAdd", paste(paste(state, "SGPt_Baseline_Matrices", sep="_"), "rda", sep="."), "to the RLImatrices GitHub repo 'data' directory,"))
			message("\t\tupdate version number/date, tag repo and commit tagged version to GitHub.\n")
		}
	} ### END END_OF_WINDOW scripts


	### Return SGP object if requested

	if (return.updated.shell) return(sgp_object)

	message(paste("Finished rliSGP", date(), "in", convertTime(timetaken(started.at)), "\n"))
} ### END rliSGP
