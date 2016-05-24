`updateSGP` <-
function(what_sgp_object=NULL,
	with_sgp_data_LONG=NULL,
	with_sgp_data_INSTRUCTOR_NUMBER=NULL,
	state=NULL,
	steps=c("prepareSGP", "analyzeSGP", "combineSGP", "summarizeSGP", "visualizeSGP", "outputSGP"),
	years=NULL,
	content_areas=NULL,
	grades=NULL,
	sgp.percentiles=TRUE,
	sgp.projections=TRUE,
	sgp.projections.lagged=TRUE,
	sgp.percentiles.baseline=TRUE,
	sgp.projections.baseline=TRUE,
	sgp.projections.lagged.baseline=TRUE,
	simulate.sgps=TRUE,
	save.old.summaries=TRUE,
	save.intermediate.results=TRUE,
	calculate.simex=NULL,
	calculate.simex.baseline=NULL,
	sgp.use.my.coefficient.matrices=NULL,
	sgp.target.scale.scores=FALSE,
	sgp.target.scale.scores.only=FALSE,
	overwrite.existing.data=FALSE,
	update.old.data.with.new=TRUE,
	sgPlot.demo.report=TRUE,
	plot.types=c("bubblePlot", "studentGrowthPlot", "growthAchievementPlot"),
	outputSGP.output.type=c("LONG_Data", "LONG_FINAL_YEAR_Data", "WIDE_Data", "INSTRUCTOR_Data"),
	sgp.config=NULL,
	goodness.of.fit.print=TRUE,
	parallel.config=NULL,
	sgp.sqlite=NULL,
	SGPt=NULL,
	sgp.percentiles.equated=NULL,
	sgp.percentiles.equating.method=NULL,
	sgp.percentiles.calculate.sgps=TRUE,
	fix.duplicates=NULL,
	...) {

	SGPstateData <- SGP::SGPstateData ### Needed due to possible assignment of values to SGPstateData

	started.at <- proc.time()
	messageSGP(paste("\nStarted updateSGP", date()), "\n")
	messageSGP(match.call())


	### Create state (if NULL) from sgp_object (if possible)

	if (is.null(state)) {
		tmp.name <- toupper(gsub("_", " ", deparse(substitute(what_sgp_object))))
		state <- getStateAbbreviation(tmp.name, "updateSGP")
	}

	if (!is.null(calculate.simex) | !is.null(calculate.simex.baseline)) {
		if (is.null(SGPstateData[[state]][["Assessment_Program_Information"]][["CSEM"]])) {
			message("\tNOTE: CSEMs are required in 'SGPstateData' (either as a data.frame of CSEMs or as a variable name of CSEMsin @Data) to produce SIMEX corrected SGPs. SIMEX corrected SGPs will NOT be calculated.")
			calculate.simex <- calculate.simex.baseline <- NULL
		}
	}

	if (identical(calculate.simex, TRUE)) {
		##  Enforce that simex.use.my.coefficient.matrices must be TRUE for updating COHORT SIMEX SGPs (ONLY WHEN USING PRE-EXISTING COEFFICIENT MATRICES)
		if (is.character(csem.variable <- SGPstateData[[state]][["Assessment_Program_Information"]][["CSEM"]])) {
			calculate.simex <- list(csem.data.vnames=csem.variable, lambda=seq(0,2,0.5), simulation.iterations=75, simex.sample.size=5000, extrapolation="linear", save.matrices=TRUE, simex.use.my.coefficient.matrices = sgp.use.my.coefficient.matrices)
		} else 	calculate.simex <- list(state=state, lambda=seq(0,2,0.5), simulation.iterations=75, simex.sample.size=5000, extrapolation="linear", save.matrices=TRUE, simex.use.my.coefficient.matrices = sgp.use.my.coefficient.matrices)
	}

	if (identical(calculate.simex.baseline, TRUE)) {
		##  Enforce that simex.use.my.coefficient.matrices must be TRUE for updating BASELINE SIMEX SGPs
		if (is.character(csem.variable <- SGPstateData[[state]][["Assessment_Program_Information"]][["CSEM"]])) {
			calculate.simex.baseline <- list(csem.data.vnames=csem.variable, lambda=seq(0,2,0.5), simulation.iterations=75, simex.sample.size=5000, extrapolation="linear", save.matrices=TRUE, simex.use.my.coefficient.matrices = TRUE)
		} else 	calculate.simex.baseline <- list(state=state, lambda=seq(0,2,0.5), simulation.iterations=75, simex.sample.size=5000, extrapolation="linear", save.matrices=TRUE, simex.use.my.coefficient.matrices = TRUE)
	}

	### Utility functions

	"%w/o%" <- function(x,y) x[!x %in% y]

	### Argument checks

	if (is.null(what_sgp_object)) {
		stop("\tNOTE: Argument 'what_sgp_object' must be supplied to updateSGP (at a minimum). See man page for 'updateSGP' for details.")
	}

	if (is.null(with_sgp_data_LONG)) {
		sgp.use.my.coefficient.matrices <- TRUE
	}
	if (identical(sgp.use.my.coefficient.matrices, FALSE)) {
		sgp.use.my.coefficient.matrices <- NULL
	}

	matrix.names <- names(what_sgp_object@SGP[['Coefficient_Matrices']])


	##############################################################################
	### DOESN'T supply 'with_sgp_data_LONG'
	##############################################################################

	if (is.null(with_sgp_data_LONG)) {

		tmp.years <- tmp.content_areas.years <- list()

		if (is.null(content_areas)) {
			content_areas <- sort(unique(sapply(strsplit(matrix.names, "[.]"), '[', 1)))
		}

		for (i in content_areas) {
			tmp.content_area.matrix.names <- grep(i, matrix.names, value=TRUE)
			tmp.years[[i]] <- sort(unique(sapply(strsplit(tmp.content_area.matrix.names, "[.]"), '[', 2))) %w/o% "BASELINE"
		}

		if (!is.null(years)) {
			for (i in content_areas) {
				tmp.years[[i]] <- intersect(tmp.years[[i]], years)
			}
		}

		for (i in names(tmp.years)) {
			tmp.content_areas.years[[i]] <- paste(i, tmp.years[[i]], sep=".")
		}
		tmp.content_areas.years <- as.character(unlist(tmp.content_areas.years))

		### NULL out existing results to be re-calculated

		what_sgp_object@SGP[['Goodness_of_Fit']][grep(paste(tmp.content_areas.years, collapse="|"), names(what_sgp_object@SGP[['Goodness_of_Fit']]))] <- NULL
		what_sgp_object@SGP[['SGPercentiles']][grep(paste(tmp.content_areas.years, collapse="|"), names(what_sgp_object@SGP[['SGPercentiles']]))] <- NULL
		what_sgp_object@SGP[['SGProjections']][grep(paste(tmp.content_areas.years, collapse="|"), names(what_sgp_object@SGP[['SGProjections']]))] <- NULL
		what_sgp_object@SGP[['Simulated_SGPs']][grep(paste(tmp.content_areas.years, collapse="|"), names(what_sgp_object@SGP[['Simulated_SGPs']]))] <- NULL

		### Add in INSTRUCTOR_NUMBER data is supplied

		if (!is.null(with_sgp_data_INSTRUCTOR_NUMBER)) {
			what_sgp_object@Data_Supplementary[['INSTRUCTOR_NUMBER']] <-
				data.table(rbindlist(list(what_sgp_object@Data_Supplementary[['INSTRUCTOR_NUMBER']], with_sgp_data_INSTRUCTOR_NUMBER), fill=TRUE),
					key=c("ID", "CONTENT_AREA", "YEAR"))
		}

		### Update results

		sgp_object <- abcSGP(
					sgp_object=what_sgp_object,
					state=state,
					steps=steps,
					years=years,
					content_areas=content_areas,
					grades=grades,
					sgp.percentiles=sgp.percentiles,
					sgp.projections=sgp.projections,
					sgp.projections.lagged=sgp.projections.lagged,
					sgp.percentiles.baseline=sgp.percentiles.baseline,
					sgp.projections.baseline=sgp.projections.baseline,
					sgp.projections.lagged.baseline=sgp.projections.lagged.baseline,
					save.intermediate.results=save.intermediate.results,
					save.old.summaries=save.old.summaries,
					sgPlot.demo.report=sgPlot.demo.report,
					sgp.use.my.coefficient.matrices=sgp.use.my.coefficient.matrices,
					calculate.simex = calculate.simex,
					calculate.simex.baseline=calculate.simex.baseline,
					simulate.sgps = simulate.sgps,
					sgp.target.scale.scores=sgp.target.scale.scores,
					sgp.target.scale.scores.only=sgp.target.scale.scores.only,
					sgp.config=sgp.config,
					plot.types=plot.types,
					goodness.of.fit.print=goodness.of.fit.print,
					outputSGP.output.type=outputSGP.output.type,
					SGPt=SGPt,
					sgp.percentiles.equated=sgp.percentiles.equated,
					sgp.percentiles.equating.method=sgp.percentiles.equating.method,
					sgp.percentiles.calculate.sgps=sgp.percentiles.calculate.sgps,
					parallel.config=parallel.config,
					...
					)

		### Print finish and return SGP object

		message(paste("Finished updateSGP", date(), "in", convertTime(timetaken(started.at)), "\n"))
		return(sgp_object)
	} ### END is.null(with_sgp_data_LONG)


	#############################################################################################
	### DOES supply 'with_sgp_data_LONG'
	#############################################################################################

	if (!is.null(with_sgp_data_LONG)) {

		HIGH_NEED_STATUS <- YEAR <- ID <- VALID_CASE <- CONTENT_AREA <- FIRST_OBSERVATION <- LAST_OBSERVATION <- NULL
		tmp_sgp_object <- prepareSGP(with_sgp_data_LONG, state=state, create.additional.variables=FALSE, fix.duplicates=NULL)
		if (!is.null(sgp.config)) years <- unique(sapply(lapply(sgp.config, '[[', 'sgp.panel.years'), tail, 1))
		if (is.null(years)) update.years <- sort(unique(tmp_sgp_object@Data$YEAR)) else update.years <- years
		if (is.null(content_areas)) update.content_areas <- sort(unique(tmp_sgp_object@Data["VALID_CASE"]$CONTENT_AREA)) else update.content_areas <- content_areas
		if (is.null(grades)) update.grades <- sort(unique(tmp_sgp_object@Data["VALID_CASE"]$GRADE)) else update.grades <- grades

		if (overwrite.existing.data) {
				what_sgp_object@Data <- rbindlist(list(what_sgp_object@Data[which(YEAR!=update.years)], tmp_sgp_object@Data), fill=TRUE)
				what_sgp_object@SGP[['Goodness_of_Fit']][grep(update.years, names(what_sgp_object@SGP[['Goodness_of_Fit']]))] <- NULL
				what_sgp_object@SGP[['Linkages']][grep(update.years, names(what_sgp_object@SGP[['Linkages']]))] <- NULL
				what_sgp_object@SGP[['SGPercentiles']][grep(update.years, names(what_sgp_object@SGP[['SGPercentiles']]))] <- NULL
				what_sgp_object@SGP[['SGProjections']][grep(update.years, names(what_sgp_object@SGP[['SGProjections']]))] <- NULL
				what_sgp_object@SGP[['Simulated_SGPs']][grep(update.years, names(what_sgp_object@SGP[['Simulated_SGPs']]))] <- NULL
				if (is.null(sgp.use.my.coefficient.matrices)) {
					what_sgp_object@SGP[['Coefficient_Matrices']][grep(update.years, names(what_sgp_object@SGP[['Coefficient_Matrices']]))] <- NULL
				}

			if ("HIGH_NEED_STATUS" %in% names(what_sgp_object@Data)) {
				what_sgp_object@Data[, HIGH_NEED_STATUS := NULL]
				what_sgp_object <- suppressMessages(prepareSGP(what_sgp_object, state=state, fix.duplicates=fix.duplicates))
			}

			what_sgp_object <- abcSGP(
						what_sgp_object,
						steps=steps,
						years=update.years,
						content_areas=update.content_areas,
						grades=update.grades,
						state=state,
						sgp.percentiles=sgp.percentiles,
						sgp.projections=sgp.projections,
						sgp.projections.lagged=sgp.projections.lagged,
						sgp.percentiles.baseline=sgp.percentiles.baseline,
						sgp.projections.baseline=sgp.projections.baseline,
						sgp.projections.lagged.baseline=sgp.projections.lagged.baseline,
						save.intermediate.results=save.intermediate.results,
						save.old.summaries=save.old.summaries,
						sgPlot.demo.report=sgPlot.demo.report,
						sgp.use.my.coefficient.matrices=sgp.use.my.coefficient.matrices,
						calculate.simex = calculate.simex,
						calculate.simex.baseline=calculate.simex.baseline,
						simulate.sgps = simulate.sgps,
						sgp.target.scale.scores=sgp.target.scale.scores,
						sgp.target.scale.scores.only=sgp.target.scale.scores.only,
						sgp.config=sgp.config,
						plot.types=plot.types,
						goodness.of.fit.print=goodness.of.fit.print,
						outputSGP.output.type=outputSGP.output.type,
						SGPt=SGPt,
						sgp.percentiles.equated=sgp.percentiles.equated,
						sgp.percentiles.equating.method=sgp.percentiles.equating.method,
						sgp.percentiles.calculate.sgps=sgp.percentiles.calculate.sgps,
						parallel.config=parallel.config,
						...)

			### Print finish and return SGP object

			message(paste("Finished updateSGP", date(), "in", convertTime(timetaken(started.at)), "\n"))
			return(what_sgp_object)

		} else {
			if (!is.null(sgp.use.my.coefficient.matrices)) {
				# Extract score histories.  Don't use CONTENT_AREA due to potential use of EOCT course progressions.
				tmp.long.data <- rbindlist(list(data.table(what_sgp_object@Data, key=c("VALID_CASE", "ID"))[
					unique(data.table(tmp_sgp_object@Data, key=c("VALID_CASE", "ID"))[,list(VALID_CASE, ID)]), nomatch=0], tmp_sgp_object@Data), fill=TRUE)
				if ("YEAR_WITHIN" %in% names(tmp.long.data)) {
					tmp.long.data$FIRST_OBSERVATION <- NULL
					tmp.long.data$LAST_OBSERVATION <- NULL
				}
				tmp.sgp_object.update <- prepareSGP(tmp.long.data, state=state, create.additional.variables=FALSE, fix.duplicates=fix.duplicates)
				tmp.sgp_object.update@SGP$Coefficient_Matrices <- what_sgp_object@SGP$Coefficient_Matrices

				if (is.null(SGPstateData[[state]][["SGP_Configuration"]])) {
					SGPstateData[[state]][["SGP_Configuration"]] <- list(return.prior.scale.score.standardized=FALSE)
				} else SGPstateData[[state]][["SGP_Configuration"]][["return.prior.scale.score.standardized"]] <- FALSE

				tmp.sgp_object.update <- analyzeSGP(
							tmp.sgp_object.update,
							years=update.years,
							content_areas=update.content_areas,
							grades=update.grades,
							state=state,
							sgp.percentiles=sgp.percentiles,
							sgp.projections=sgp.projections,
							sgp.projections.lagged=sgp.projections.lagged,
							sgp.percentiles.baseline=sgp.percentiles.baseline,
							sgp.projections.baseline=sgp.projections.baseline,
							sgp.projections.lagged.baseline=sgp.projections.lagged.baseline,
							sgp.use.my.coefficient.matrices=sgp.use.my.coefficient.matrices,
							calculate.simex=calculate.simex,
							calculate.simex.baseline=calculate.simex.baseline,
							simulate.sgps=simulate.sgps,
							sgp.config=sgp.config,
							goodness.of.fit.print=FALSE,
							sgp.sqlite=sgp.sqlite,
							SGPt=SGPt,
							sgp.percentiles.equated=sgp.percentiles.equated,
							sgp.percentiles.equating.method=sgp.percentiles.equating.method,
							sgp.percentiles.calculate.sgps=sgp.percentiles.calculate.sgps,
							parallel.config=parallel.config,
							...)

				if ("combineSGP" %in% steps) {
					tmp.sgp_object.update <- suppressMessages(combineSGP(tmp.sgp_object.update,
						state=state,
						sgp.percentiles=sgp.percentiles,
						sgp.projections=sgp.projections,
						sgp.projections.lagged=sgp.projections.lagged,
						sgp.percentiles.baseline=sgp.percentiles.baseline,
						sgp.projections.baseline=sgp.projections.baseline,
						sgp.projections.lagged.baseline=sgp.projections.lagged.baseline,
						sgp.target.scale.scores=sgp.target.scale.scores,
						sgp.target.scale.scores.only=sgp.target.scale.scores.only,
						SGPt=SGPt,
						parallel.config=parallel.config))
				}

				### Output of INTERMEDIATE results including full student history

				dir.create(file.path("Data", "Updated_Data"), recursive=TRUE, showWarnings=FALSE)
				tmp.file.name <- paste(gsub(" ", "_", toupper(getStateAbbreviation(state, type="name"))), "SGP_Update", paste(update.years, collapse=","), sep="_")
				assign(tmp.file.name, tmp.sgp_object.update)
				save(list=tmp.file.name, file=file.path("Data", "Updated_Data", paste(tmp.file.name, "Rdata", sep=".")))
				outputSGP(tmp.sgp_object.update, state=state, output.type=union(outputSGP.output.type, intersect(outputSGP.output.type, "LONG_FINAL_YEAR_Data")),
					outputSGP.directory=file.path("Data", "Updated_Data"))

				### Merge update with original SGP object

				what_sgp_object@Data <- data.table(rbindlist(list(what_sgp_object@Data, tmp_sgp_object@Data), fill=TRUE), key=getKey(what_sgp_object))
				if ("HIGH_NEED_STATUS" %in% names(what_sgp_object@Data)) {
					what_sgp_object@Data[, HIGH_NEED_STATUS := NULL]
					what_sgp_object <- suppressMessages(prepareSGP(what_sgp_object, state=state, fix.duplicates=fix.duplicates))
				}
				what_sgp_object@SGP <- mergeSGP(what_sgp_object@SGP, tmp.sgp_object.update@SGP)

				### Filter out any exact duplicates in projections and percentiles (duplicates because score histories not subsetted based on CONTENT_AREA).
				if (sgp.projections | sgp.projections.lagged | sgp.projections.baseline | sgp.projections.lagged.baseline) {
					for (ca in grep(update.years, names(what_sgp_object@SGP[["SGProjections"]]), value=TRUE)) {
						tmp_proj <- data.table(what_sgp_object@SGP[["SGProjections"]][[ca]])
						setkey(tmp_proj)
						what_sgp_object@SGP[["SGProjections"]][[ca]]<- tmp_proj[!duplicated(tmp_proj)]
					}
				}

				if (sgp.percentiles | sgp.percentiles.baseline) {
					for (ca in grep(update.years, names(what_sgp_object@SGP[["SGPercentiles"]]), value=TRUE)) {
						tmp_sgp <- data.table(what_sgp_object@SGP[["SGPercentiles"]][[ca]])
						setkeyv(tmp_sgp, names(tmp_sgp)[grep("ID|SGP|SGP_NORM_GROUP", names(tmp_sgp))])
						what_sgp_object@SGP[["SGPercentiles"]][[ca]]<- tmp_sgp[!duplicated(tmp_sgp)]
					}
				}

				if ("combineSGP" %in% steps) {
					what_sgp_object <- combineSGP(
						what_sgp_object,
						years=update.years,
						state=state,
						sgp.percentiles=sgp.percentiles,
						sgp.projections=sgp.projections,
						sgp.projections.lagged=sgp.projections.lagged,
						sgp.percentiles.baseline=sgp.percentiles.baseline,
						sgp.projections.baseline=sgp.projections.baseline,
						sgp.projections.lagged.baseline=sgp.projections.lagged.baseline,
						sgp.target.scale.scores=sgp.target.scale.scores,
						sgp.target.scale.scores.only=sgp.target.scale.scores.only,
						SGPt=SGPt,
						parallel.config=parallel.config)
				}


				### Add in INSTRUCTOR_NUMBER data for summarizeSGP if supplied

				if (!is.null(with_sgp_data_INSTRUCTOR_NUMBER)) {
					what_sgp_object@Data_Supplementary[['INSTRUCTOR_NUMBER']] <-
						data.table(rbindlist(list(what_sgp_object@Data_Supplementary[['INSTRUCTOR_NUMBER']], with_sgp_data_INSTRUCTOR_NUMBER), fill=TRUE),
							key=c("ID", "CONTENT_AREA", "YEAR"))
				}

				if ("summarizeSGP" %in% steps) what_sgp_object <- summarizeSGP(what_sgp_object, state=state, parallel.config=parallel.config)
				if ("visualizeSGP" %in% steps) visualizeSGP(what_sgp_object, state=state, plot.types=plot.types, sgPlot.demo.report=sgPlot.demo.report)
				if ("outputSGP" %in% steps) outputSGP(what_sgp_object, state=state, output.type=outputSGP.output.type)


				### Re-establish FIRST_ & LAST_OBSERVATION variables
				if ("YEAR_WITHIN" %in% names(what_sgp_object@Data)) {
					what_sgp_object@Data[, LAST_OBSERVATION := NULL]
					what_sgp_object@Data[, FIRST_OBSERVATION := NULL]
					what_sgp_object <- suppressMessages(prepareSGP(what_sgp_object, state=state, create.additional.variables=FALSE, fix.duplicates=fix.duplicates))
				}

				### Print finish and return SGP object

				message(paste("Finished updateSGP", date(), "in", convertTime(timetaken(started.at)), "\n"))
				return(what_sgp_object)
			} else {
				if (update.old.data.with.new) {
					what_sgp_object@Data <- data.table(rbindlist(list(what_sgp_object@Data, tmp_sgp_object@Data), fill=TRUE), key=getKey(what_sgp_object))
				} else {
					what_sgp_object@Data <- data.table(rbindlist(list(what_sgp_object@Data[which(ID %in% tmp_sgp_object@Data$ID)], tmp_sgp_object@Data), fill=TRUE),
						key=getKey(what_sgp_object))
				}

				### Re-establish FIRST_ & LAST_OBSERVATION variables
				if ("YEAR_WITHIN" %in% names(what_sgp_object@Data)) {
					what_sgp_object@Data[, LAST_OBSERVATION := NULL]
					what_sgp_object@Data[, FIRST_OBSERVATION := NULL]
				}

				### prepareSGP
				if ("HIGH_NEED_STATUS" %in% names(what_sgp_object@Data)) {
					what_sgp_object@Data[, HIGH_NEED_STATUS := NULL]
					what_sgp_object <- prepareSGP(what_sgp_object, state=state, fix.duplicates=fix.duplicates)
				} else {
					what_sgp_object <- prepareSGP(what_sgp_object, state=state, create.additional.variables=FALSE, fix.duplicates=fix.duplicates)
				}

				### Add in INSTRUCTOR_NUMBER data for summarizeSGP if supplied

				if (!is.null(with_sgp_data_INSTRUCTOR_NUMBER)) {
					what_sgp_object@Data_Supplementary[['INSTRUCTOR_NUMBER']] <-
						data.table(rbindlist(list(what_sgp_object@Data_Supplementary[['INSTRUCTOR_NUMBER']], with_sgp_data_INSTRUCTOR_NUMBER), fill=TRUE),
							key=c("ID", "CONTENT_AREA", "YEAR"))
				}

				### abcSGP

				what_sgp_object <- abcSGP(
							what_sgp_object,
							steps=(steps %w/o% "prepareSGP"),
							years=update.years,
							content_areas=update.content_areas,
							grades=update.grades,
							state=state,
							sgp.percentiles=sgp.percentiles,
							sgp.projections=sgp.projections,
							sgp.projections.lagged=sgp.projections.lagged,
							sgp.percentiles.baseline=sgp.percentiles.baseline,
							sgp.projections.baseline=sgp.projections.baseline,
							sgp.projections.lagged.baseline=sgp.projections.lagged.baseline,
							save.intermediate.results=save.intermediate.results,
							save.old.summaries=save.old.summaries,
							sgPlot.demo.report=sgPlot.demo.report,
							sgp.use.my.coefficient.matrices=sgp.use.my.coefficient.matrices,
							calculate.simex = calculate.simex,
							calculate.simex.baseline=calculate.simex.baseline,
							simulate.sgps = simulate.sgps,
							sgp.target.scale.scores=sgp.target.scale.scores,
							sgp.target.scale.scores.only=sgp.target.scale.scores.only,
							sgp.config=sgp.config,
							plot.types=plot.types,
							goodness.of.fit.print=goodness.of.fit.print,
							outputSGP.output.type=outputSGP.output.type,
							SGPt=SGPt,
							sgp.percentiles.equated=sgp.percentiles.equated,
							sgp.percentiles.equating.method=sgp.percentiles.equating.method,
							sgp.percentiles.calculate.sgps=sgp.percentiles.calculate.sgps,
							parallel.config=parallel.config,
							...)

				### Print finish and return SGP object

				messageSGP(paste("Finished updateSGP", date(), "in", convertTime(timetaken(started.at)), "\n"))
				return(what_sgp_object)
			} ### END if else (!is.null(sgp.use.my.coefficient.matrices))
		} ### END if (overwrite.existing.data)
	} ### END !is.null(with_sgp_data_LONG)
} ## END updateSGP Function
