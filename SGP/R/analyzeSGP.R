`analyzeSGP` <-
function(sgp_object,
         state=NULL,
         years=NULL,
         content_areas=NULL,
         grades=NULL,
         sgp.percentiles=TRUE,
         sgp.projections=TRUE,
         sgp.projections.lagged=TRUE,
         sgp.percentiles.baseline=TRUE,
         sgp.projections.baseline=TRUE,
         sgp.projections.lagged.baseline=TRUE,
         sgp.percentiles.baseline.max.order=3,
         sgp.projections.baseline.max.order=3,
         sgp.projections.lagged.baseline.max.order=3,
         sgp.projections.max.forward.progression.years=3,
         sgp.projections.max.forward.progression.grade=NULL,
         sgp.minimum.default.panel.years=NULL,
         sgp.use.my.coefficient.matrices=NULL,
         sgp.use.my.sgp_object.baseline.coefficient.matrices=NULL,
         simulate.sgps=TRUE,
         calculate.simex=NULL,
         calculate.simex.baseline=NULL,
         goodness.of.fit.print=TRUE,
         sgp.config=NULL,
         sgp.config.drop.nonsequential.grade.progression.variables=TRUE,
         sgp.baseline.panel.years=NULL,
         sgp.baseline.config=NULL,
         trim.sgp.config=TRUE,
         parallel.config=NULL,
         verbose.output=FALSE,
         print.other.gp=NULL,
         sgp.projections.projection.unit="YEAR",
         get.cohort.data.info=FALSE,
         sgp.sqlite=NULL,
         sgp.percentiles.equated=NULL,
         sgp.percentiles.equating.method=NULL,
         sgp.percentiles.calculate.sgps=TRUE,
         SGPt=NULL,
         ...) {

	started.at <- proc.time()
	messageSGP(paste("\nStarted analyzeSGP", date()), "\n")
	messageSGP(match.call())

	VALID_CASE <- CONTENT_AREA <- YEAR <- GRADE <- ID <- YEAR_WITHIN <- SCALE_SCORE <- NULL
	SGPstateData <- SGP::SGPstateData ### Needed due to possible assignment of values to SGPstateData

	###
	### Create relevant analyzeSGP variables
	###

	### Create state (if NULL) from sgp_object (if possible)

	if (is.null(state)) {
		tmp.name <- toupper(gsub("_", " ", deparse(substitute(sgp_object))))
		state <- getStateAbbreviation(tmp.name, "analyzeSGP")
	}


	###
	### Tests associated with provided arguments
	###

	if (!(sgp.percentiles | sgp.percentiles.baseline)) {
		simulate.sgps <- FALSE
	}

	if (simulate.sgps) {
		if (is.null(SGPstateData[[state]][["Assessment_Program_Information"]][["CSEM"]])) {
			messageSGP("\tNOTE: CSEMs are required in 'SGPstateData' (either as a data.frame of CSEMs or as a variable name of CSEMsin @Data) to simulate SGPs for confidence interval calculations. SGP standard errors will not be calculated.")
			calculate.confidence.intervals <- csem.variable <- NULL
		} else {
            if (is.data.frame(SGPstateData[[state]][["Assessment_Program_Information"]][["CSEM"]])) {
                calculate.confidence.intervals <- state
                csem.variable <- NULL
            }
            if (is.character(SGPstateData[[state]][["Assessment_Program_Information"]][["CSEM"]])) {
                calculate.confidence.intervals <- csem.variable <- SGPstateData[[state]][["Assessment_Program_Information"]][["CSEM"]]
            }
		}
	} else {
		calculate.confidence.intervals <- csem.variable <- NULL
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][["sgp.config.drop.nonsequential.grade.progression.variables"]])) {
		sgp.config.drop.nonsequential.grade.progression.variables <- SGPstateData[[state]][["SGP_Configuration"]][["sgp.config.drop.nonsequential.grade.progression.variables"]]
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][["sgp.loss.hoss.adjustment"]])) {
		sgp.loss.hoss.adjustment <- SGPstateData[[state]][["SGP_Configuration"]][["sgp.loss.hoss.adjustment"]]
	} else {
		sgp.loss.hoss.adjustment <- NULL
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][["return.norm.group.scale.scores"]])) {
		return.norm.group.scale.scores <- SGPstateData[[state]][["SGP_Configuration"]][["return.norm.group.scale.scores"]]
	} else {
		return.norm.group.scale.scores <- NULL
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][["return.norm.group.dates"]])) {
		return.norm.group.dates <- SGPstateData[[state]][["SGP_Configuration"]][["return.norm.group.dates"]]
	} else {
		return.norm.group.dates <- NULL
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][["return.projection.group.scale.scores"]])) {
		return.projection.group.scale.scores <- SGPstateData[[state]][["SGP_Configuration"]][["return.projection.group.scale.scores"]]
	} else {
		return.projection.group.scale.scores <- NULL
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][["return.projection.group.scale.scores"]])) {
		return.projection.group.scale.scores <- SGPstateData[[state]][["SGP_Configuration"]][["return.projection.group.scale.scores"]]
	} else {
		return.projection.group.scale.scores <- NULL
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][["return.projection.group.dates"]])) {
		return.projection.group.dates <- SGPstateData[[state]][["SGP_Configuration"]][["return.projection.group.dates"]]
	} else {
		return.projection.group.dates <- NULL
	}

	if (!is.null(SGPstateData[[state]][["Growth"]][["Cutscores"]][["Cuts"]])) {
		percentile.trajectory.values <- unique(c(SGPstateData[[state]][["Growth"]][["Cutscores"]][["Cuts"]], 50))
	} else {
		percentile.trajectory.values <- c(35, 50, 65)
	}

	if (!is.null(SGPstateData[[state]][["Student_Report_Information"]][["Projection_Fan_Limits"]])) {
		percentile.trajectory.values <- sort(c(SGPstateData[[state]][["Student_Report_Information"]][["Projection_Fan_Limits"]], percentile.trajectory.values))
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][["sgp.projections.baseline.max.order"]])) {
		sgp.projections.baseline.max.order <- SGPstateData[[state]][["SGP_Configuration"]][["sgp.projections.baseline.max.order"]]
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][["sgp.projections.lagged.baseline.max.order"]])) {
		sgp.projections.lagged.baseline.max.order <- SGPstateData[[state]][["SGP_Configuration"]][["sgp.projections.lagged.baseline.max.order"]]
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][["return.prior.scale.score.standardized"]])) {
		return.prior.scale.score.standardized <- SGPstateData[[state]][["SGP_Configuration"]][["return.prior.scale.score.standardized"]]
	} else {
		return.prior.scale.score.standardized <- TRUE
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][["max.n.for.coefficient.matrices"]])) {
		max.n.for.coefficient.matrices <- SGPstateData[[state]][["SGP_Configuration"]][["max.n.for.coefficient.matrices"]]
	} else {
		max.n.for.coefficient.matrices <- NULL
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][["sgp.cohort.size"]]) & is.null(sgp.use.my.coefficient.matrices)) {
		tmp.cohort.size <- SGPstateData[[state]][["SGP_Configuration"]][["sgp.cohort.size"]]
	} else tmp.cohort.size <- NULL

	if (!is.null(sgp.config) && sgp.config.drop.nonsequential.grade.progression.variables) {
		sgp.config.drop.nonsequential.grade.progression.variables <- FALSE
	}

	# lower.level.parallel.config <- TRUE
	if (all(c("PERCENTILES", "TAUS") %in% names(parallel.config[['WORKERS']]))) {
		# lower.level.parallel.config <- FALSE
		# if (as.numeric(parallel.config[['WORKERS']][['PERCENTILES']])!=1) stop("Both TAUS and PERCENTILES can not be executed in Parallel at the same time.")
		if (.Platform$OS.type != "unix" | "SNOW_TEST" %in% names(parallel.config)) stop("Both TAUS and PERCENTILES can not be executed in Parallel at the same time in Windows OS or using SNOW type backends.")
		messageSGP("\n\tCAUTION:  Running higher- and lower-level processes in parallel at the same time.  Make sure you have enough CPU cores and memory to support this.\n")
	}
	if (all(c("PERCENTILES", "SIMEX") %in% names(parallel.config[['WORKERS']]))) {
		if (.Platform$OS.type != "unix" | "SNOW_TEST" %in% names(parallel.config)) stop("Both SIMEX and PERCENTILES can not be executed in Parallel at the same time in Windows OS or using SNOW type backends.")
		messageSGP("\n\tCAUTION:  Running higher- and lower-level processes in parallel at the same time.  Make sure you have enough CPU cores and memory to support this.\n")
	}

	if (all(c("BASELINE_PERCENTILES", "TAUS") %in% names(parallel.config[['WORKERS']]))) {
		if (.Platform$OS.type != "unix" | "SNOW_TEST" %in% names(parallel.config)) stop("Both TAUS and BASELINE_PERCENTILES can not be executed in Parallel at the same time in Windows OS or using SNOW type backends.")
		messageSGP("\n\tCAUTION:  Running higher- and lower-level processes in parallel at the same time.  Make sure you have enough CPU cores and memory to support this.\n")
	}
	if (all(c("BASELINE_PERCENTILES", "SIMEX") %in% names(parallel.config[['WORKERS']]))) {
		if (.Platform$OS.type != "unix" | "SNOW_TEST" %in% names(parallel.config)) stop("Both SIMEX and BASELINE_PERCENTILES can not be executed in Parallel at the same time in Windows OS or using SNOW type backends.")
		messageSGP("\n\tCAUTION:  Running higher- and lower-level processes in parallel at the same time.  Make sure you have enough CPU cores and memory to support this.\n")
	}

	if (any(c("SIMEX", "TAUS") %in% names(parallel.config[['WORKERS']]))) {
		lower.level.parallel.config <- parallel.config
		if (length(grep("SUMMARY|GA_PLOTS|SG_PLOTS", names(parallel.config[['WORKERS']]), value=TRUE, invert=TRUE)) <= 2) parallel.config <- NULL # NULL out parallel.config when passed from abcSGP, etc
	} else lower.level.parallel.config <- NULL


	if (!is.null(calculate.simex) | !is.null(calculate.simex.baseline)) {
		if (is.null(SGPstateData[[state]][["Assessment_Program_Information"]][["CSEM"]])) {
			messageSGP("\tNOTE: CSEMs are required in 'SGPstateData' (either as a data.frame of CSEMs or as a variable name of CSEMsin @Data) to produce SIMEX corrected SGPs. SIMEX corrected SGPs will NOT be calculated.")
			calculate.simex <- calculate.simex.baseline <- NULL
		}
	}

	if (identical(calculate.simex, TRUE)) {
		if (is.character(csem.variable <- SGPstateData[[state]][["Assessment_Program_Information"]][["CSEM"]])) {
			calculate.simex <- list(csem.data.vnames=csem.variable, lambda=seq(0,2,0.5), simulation.iterations=75, simex.sample.size=5000, extrapolation="linear", save.matrices=TRUE)
		} else 	{
			calculate.simex <- list(state=state, lambda=seq(0,2,0.5), simulation.iterations=75, simex.sample.size=5000, extrapolation="linear", save.matrices=TRUE)
			csem.variable <- NULL
		}
	}

	if (identical(calculate.simex.baseline, TRUE)) {
		if (is.character(csem.variable <- SGPstateData[[state]][["Assessment_Program_Information"]][["CSEM"]])) {
			calculate.simex.baseline <- list(csem.data.vnames=csem.variable, lambda=seq(0,2,0.5), simulation.iterations=75, simex.sample.size=5000, extrapolation="linear", save.matrices=TRUE)
		} else 	{
			calculate.simex.baseline <- list(state=state, lambda=seq(0,2,0.5), simulation.iterations=75, simex.sample.size=5000, extrapolation="linear", save.matrices=TRUE)
			csem.variable <- NULL
		}
	}

	if (is.null(sgp.minimum.default.panel.years) & !is.null(SGPstateData[[state]][["SGP_Configuration"]][['sgp.minimum.default.panel.years']])) {
		sgp.minimum.default.panel.years <- SGPstateData[[state]][["SGP_Configuration"]][['sgp.minimum.default.panel.years']]
	}

	if (is.null(sgp.minimum.default.panel.years) & is.null(SGPstateData[[state]][["SGP_Configuration"]][['sgp.minimum.default.panel.years']])) {
		if (length(unique(sgp_object@Data$YEAR))==2) {
			sgp.minimum.default.panel.years <- 2
			messageSGP("\tNOTE: Only two years of data present. Minimum default of 3 years of panel data for SGP analyses changed to 2. Please confirm this is consistent with analyses you wish to perform.")
		} else {
			sgp.minimum.default.panel.years <- 3
		}
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][['sgp.projections.max.forward.progression.grade']])) {
		sgp.projections.max.forward.progression.grade <- SGPstateData[[state]][["SGP_Configuration"]][['sgp.projections.max.forward.progression.grade']]
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][['print.other.gp']])) {
		print.other.gp <- SGPstateData[[state]][["SGP_Configuration"]][['print.other.gp']]
	}
	if (is.null(print.other.gp)) print.other.gp <- FALSE

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][['sgp.projections.max.forward.progression.years']])) {
		if (SGPstateData[[state]][["SGP_Configuration"]][['sgp.projections.max.forward.progression.years']]==FALSE) {
			sgp.projections.max.forward.progression.years <- NULL
		} else {
			sgp.projections.max.forward.progression.years <- SGPstateData[[state]][["SGP_Configuration"]][['sgp.projections.max.forward.progression.years']]
		}
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][['sgp.projections.projection.unit']])) {
		sgp.projections.projection.unit <- SGPstateData[[state]][["SGP_Configuration"]][['sgp.projections.projection.unit']]
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][['sgp.projections.projection.unit.label']])) {
		sgp.projections.projection.unit.label <- SGPstateData[[state]][["SGP_Configuration"]][['sgp.projections.projection.unit.label']]
	} else {
		sgp.projections.projection.unit.label <- sgp.projections.projection.unit
	}

	if (is.character(goodness.of.fit.print)) {
		if (goodness.of.fit.print =="GROB") {
			goodness.of.fit.print <- FALSE
			goodness.of.fit.print.arg <- state
		} else goodness.of.fit.print <- as.logical(goodness.of.fit.print)
	} else {
		if (!goodness.of.fit.print){
			goodness.of.fit.print.arg <- FALSE
		} else {
            if (identical(SGPstateData[[state]][["SGP_Configuration"]][["goodness.of.fit.achievement.level.prior"]], FALSE)) {   ### For RLI and RLI_UK
			         goodness.of.fit.print.arg <- TRUE
            } else {
			        goodness.of.fit.print.arg <- state
            }
		}
	}

	if (!is.null(SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][["Year"]])) {
		if (SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][["Year"]]!={tmp.last.year <- tail(sort(unique(sgp_object@Data$YEAR)), 1)}) {
			sgp.percentiles.equated <- FALSE
            if ("SCALE_SCORE_EQUATED" %in% names(sgp_object@Data)) sgp_object@Data[["SCALE_SCORE_EQUATED"]][sgp_object@Data[["YEAR"]] >= tmp.last.year] <- sgp_object@Data[["SCALE_SCORE"]][sgp_object@Data[["YEAR"]] >= tmp.last.year]
		} else {
			if (!identical(sgp.percentiles.equated, FALSE)) sgp.percentiles.equated <- TRUE
		}
	} else {
		if (identical(sgp.percentiles.equated, TRUE)) {
			messageSGP("\t\tNOTE: 'sgp.percentiles.equated' has been set to TRUE but no meta-data exists in 'SGPstateData' associated with that assessment transition. Equated/linked SGP analyses require meta-data embedded in 'SGPstateData' to correctly work. Contact package administrators on how such data can be added to the package.")
		}
		sgp.percentiles.equated <- FALSE
	}

	if (!is.null(SGPt)) {
		if (identical(SGPt, TRUE)) SGPt <- "DATE"
		if (!all(SGPt %in% names(sgp_object@Data))) {
			tmp.messages <- c(tmp.messages, "\t\tNOTE: Variables", paste(SGPt, collapse=", "), "are not all contained in the supplied 'sgp_object@Data'. 'SGPt' is set to NULL.\n")
			SGPt <- NULL
		}
        SGPt.max.time <- SGPstateData[[state]][['SGP_Configuration']][['SGPt.max.time']]
	} else {
        SGPt.max.time <- NULL
    }

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][['sgp.use.my.sgp_object.baseline.coefficient.matrices']]) && is.null(sgp.use.my.sgp_object.baseline.coefficient.matrices)) {
		sgp.use.my.sgp_object.baseline.coefficient.matrices <- SGPstateData[[state]][["SGP_Configuration"]][['sgp.use.my.sgp_object.baseline.coefficient.matrices']]
	}

	if (identical(sgp.use.my.sgp_object.baseline.coefficient.matrices, FALSE)) sgp.use.my.sgp_object.baseline.coefficient.matrices <- NULL

	if (!is.null(sgp.use.my.sgp_object.baseline.coefficient.matrices) && length(grep("BASELINE", names(sgp_object@SGP$Coefficient_Matrices)))==0) {
		sgp.use.my.sgp_object.baseline.coefficient.matrices <- NULL
	}

	if (!is.null(SGPstateData[[state]][["SGP_Configuration"]][['lagged.percentile.trajectory.values']])) {
		lagged.percentile.trajectory.values <- SGPstateData[[state]][["SGP_Configuration"]][['lagged.percentile.trajectory.values']]
	} else {
        lagged.percentile.trajectory.values <- NULL
    }


	###
	### Utility functions
	###

	## Function to export/print goodness of fit results as pdf files to directory Goodness_of_Fit

	gof.print <- function(sgp_object) {
		if (length(sgp_object@SGP[["Goodness_of_Fit"]]) > 0) {
			for (i in names(sgp_object@SGP[["Goodness_of_Fit"]])) {
				dir.create(paste("Goodness_of_Fit/", i, sep=""), recursive=TRUE, showWarnings=FALSE)
					for (output.format in c("PDF", "PNG")) {
						for (j in names(sgp_object@SGP[["Goodness_of_Fit"]][[i]])) {
							tmp.path <- file.path("Goodness_of_Fit", i, j)
							if (!identical(.Platform$OS.type, "unix") & nchar(tmp.path) > 250) {
								tmp.content_area <- unlist(strsplit(j, "[.]"))[1]
								tmp.path <- gsub(tmp.content_area, substr(tmp.content_area, 1, 1), tmp.path)
							}
							if (output.format=="PDF") {
								pdf(file=paste(tmp.path, ".pdf", sep=""), width=8.5, height=11)
							}
							if (output.format=="PNG") {
								Cairo(file=paste(tmp.path, ".png", sep=""),
								      width=8.5, height=11, units="in", dpi=144, pointsize=10.5, bg="transparent")
							}
							grid.draw(sgp_object@SGP[["Goodness_of_Fit"]][[i]][[j]])
							dev.off()
						}
					}
				}
		} else {
			messageSGP("\tNOTE: No Goodness of Fit tables available to print. No tables will be produced.")
		}
	}

	## Function to merge coefficient matrices from coefficient matrix productions

	merge.coefficient.matrices <- function(list.of.matrices, simex=FALSE) {
		tmp.list <- list()
		tmp.coefficient.matrices <- unlist(list.of.matrices, recursive=FALSE)
		if (simex) {
			for (tmp.names in unique(names(tmp.coefficient.matrices))) {
				tmp1 <- unlist(tmp.coefficient.matrices[grep(tmp.names, names(tmp.coefficient.matrices))], recursive=FALSE)
				names(tmp1) <- sapply(strsplit(names(tmp1), "[.]"), function(x) x[4])
				tmp.list[[tmp.names]] <- tmp1
			}
		} else {
			for (tmp.names in unique(names(tmp.coefficient.matrices))) {
				tmp1 <- unlist(tmp.coefficient.matrices[grep(tmp.names, names(tmp.coefficient.matrices))], recursive=FALSE)
				names(tmp1) <- sapply(strsplit(names(tmp1), "[.]"), function(x) x[3])
				tmp.list[[tmp.names]] <- tmp1
			}
		}
	tmp.list
	}

	get.simulate.sgps.arg <- function(calculate.confidence.intervals, state, sgp.iter) {
		if (is.null(calculate.confidence.intervals) || calculate.confidence.intervals==state) {
			return(calculate.confidence.intervals)
		} else {
			return(gsub("[.]+$", "", paste(calculate.confidence.intervals, tail(sgp.iter[['sgp.panel.years']], 1), tail(sgp.iter[['sgp.content.areas']], 1),  tail(sgp.iter[['sgp.panel.years.within']], 1), sep=".")))
		}
	}

	get.calculate.simex.arg <- function(calculate.simex, sgp.iter) {
		if (is.null(calculate.simex)) return(NULL) # If not NULL, must be a list
		if (is.null(calculate.simex$csem.data.vnames)) return(calculate.simex)
		calculate.simex$csem.data.vnames <- gsub("[.]+$", "", paste(calculate.simex$csem.data.vnames, sgp.iter[['sgp.panel.years']], sgp.iter[['sgp.content.areas']],  sgp.iter[['sgp.panel.years.within']], sep="."))
		return(calculate.simex)
	}

	selectCoefficientMatrices <- function(tmp_sgp_object, coefficient.matrix.type=NULL) {

		if (is.null(coefficient.matrix.type)) {
			return(tmp_sgp_object[['Coefficient_Matrices']][
				setdiff(names(tmp_sgp_object[['Coefficient_Matrices']]), grep('BASELINE|EQUATED', names(tmp_sgp_object[['Coefficient_Matrices']]), value=TRUE))])
		}

		if (coefficient.matrix.type=="BASELINE") {
			return(tmp_sgp_object[["Coefficient_Matrices"]][grep("BASELINE", names(tmp_sgp_object[["Coefficient_Matrices"]]))])
		}

		if (coefficient.matrix.type=="EQUATED") {
			return(tmp_sgp_object[["Coefficient_Matrices"]][grep("EQUATED", names(tmp_sgp_object[["Coefficient_Matrices"]]))])
		}
	}

	#######################################################################################################################
	##   Set up the temporary sgp list object.  Fill with necessary old results if they exist first.
	##   Create subset of @Data containing essential data elements for analyses
	#######################################################################################################################

	if (sgp.percentiles.equated) {
		year.for.equate <- tail(sort(unique(sgp_object@Data$YEAR)), 1)

        if (is.null(SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][['Baseline_Projections_in_Transition_Year']]) &
            (sgp.percentiles.baseline | sgp.projections.baseline | sgp.projections.lagged.baseline)) {
                messageSGP("\tNOTE: Baseline SGP related analyses are not possible across an assessment transition with equating. Arguments related to baseline analyses are set to FALSE.")
                sgp.percentiles.baseline <- sgp.projections.baseline <- sgp.projections.lagged.baseline <- FALSE
        }

        if (is.null(sgp.use.my.coefficient.matrices)) {
            Scale_Score_Linkages <- list()
            dir.create(file.path("Data", paste("Linkages", year.for.equate, sep="_"), "Figures"), recursive=TRUE, showWarnings=FALSE)
		    content_areas.for.equate <- unique(sgp_object@Data[YEAR==year.for.equate]$CONTENT_AREA)

            if (is.null(sgp.percentiles.equating.method)) {
                messageSGP("\tNOTE: Analyses involving equating will be performed using each of: 'identity', 'mean', 'linear', and 'equipercentile' methods. See documentation associated with the 'sgp.percentiles.equating.method' argument in 'analyzeSGP'.")
                sgp.percentiles.equating.method <- c("identity", "mean", "linear", "equipercentile")
            }

            if (!identical(years, year.for.equate)) {
                messageSGP(paste("\tNOTE: Analyses involving equating only occur in most recent year. 'years' argument changed to ", year.for.equate, ".", sep=""))
                years <- year.for.equate
            }

            if (!all(paste(content_areas.for.equate, year.for.equate, sep=".") %in% names(SGPstateData[[state]][['Achievement']][['Knots_Boundaries']]))) {
                tmp.knots.boundaries <- createKnotsBoundaries(sgp_object@Data[YEAR==year.for.equate])
                names(tmp.knots.boundaries) <- paste(names(tmp.knots.boundaries), year.for.equate, sep=".")
                SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]] <- c(SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]], tmp.knots.boundaries)
                assign(paste(state, "Knots_Boundaries", sep="_"), SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]])
                save(list=paste(state, "Knots_Boundaries", sep="_"), file=paste(state, "Knots_Boundaries.Rdata", sep="_"))
                messageSGP(paste("\tNOTE: Knots and Boundaries do not exist for ", year.for.equate, " in state provided.\n\tThey have been produced, embedded in SGPstateData, and are available using state=", state, " for subsequent analyses and saved to your working directory '", getwd(), "'.", sep=""))
            }

            data.for.equate <- copy(sgp_object@Data)
            sgp_object@SGP[['Linkages']] <- Linkages <- equateSGP(data.for.equate, state, year.for.equate, sgp.percentiles.equating.method)
            setkey(data.for.equate, VALID_CASE, CONTENT_AREA, YEAR, GRADE, SCALE_SCORE)
            for (conversion.type.iter in c("OLD_TO_NEW", "NEW_TO_OLD")) {
                for (sgp.percentiles.equating.method.iter in sgp.percentiles.equating.method) {
                    data.for.equate <- convertScaleScore(data.for.equate, year.for.equate, sgp_object@SGP[['Linkages']],
			            conversion.type=conversion.type.iter, sgp.percentiles.equating.method.iter, state)
                    Scale_Score_Linkages[[conversion.type.iter]][[toupper(sgp.percentiles.equating.method.iter)]] <-
                        unique(data.for.equate)[!is.na(SCALE_SCORE) & VALID_CASE=="VALID_CASE", intersect(names(data.for.equate),
                            c("CONTENT_AREA", "YEAR", "GRADE", "SCALE_SCORE", "SCALE_SCORE_ACTUAL", paste("SCALE_SCORE_EQUATED", toupper(sgp.percentiles.equating.method.iter), conversion.type.iter, sep="_"))), with=FALSE]
                    write.table(Scale_Score_Linkages[[conversion.type.iter]][[toupper(sgp.percentiles.equating.method.iter)]],
                        file=paste(paste("Data/", paste("Linkages", year.for.equate, sep="_"), "/", sep=""), paste(gsub(" ", "_", getStateAbbreviation(state, type="LONG")), "Scale_Score_Linkages", toupper(sgp.percentiles.equating.method.iter), conversion.type.iter, sep="_"), ".txt", sep=""), row.names=FALSE, na="", quote=FALSE, sep="|")
                    linkagePlot(Scale_Score_Linkages[[conversion.type.iter]][[toupper(sgp.percentiles.equating.method.iter)]], conversion.type.iter, sgp.percentiles.equating.method.iter, year.for.equate, state)
                }
            }
            save(Linkages, file=paste(paste("Data/", paste("Linkages", year.for.equate, sep="_"), "/", sep=""), "Linkages.Rdata", sep=""))
            assign(paste(gsub(" ", "_", getStateAbbreviation(state, type="LONG")), "Scale_Score_Linkages", sep="_"), Scale_Score_Linkages)
            save(list=paste(gsub(" ", "_", getStateAbbreviation(state, type="LONG")), "Scale_Score_Linkages", sep="_"),
                file=paste(paste("Data/", paste("Linkages", year.for.equate, sep="_"), "/", sep=""), paste(gsub(" ", "_", getStateAbbreviation(state, type="LONG")), "Scale_Score_Linkages", sep="_"), ".Rdata", sep=""))
                setkey(data.for.equate, VALID_CASE, CONTENT_AREA, YEAR, ID)
            data.for.equate <- data.for.equate[,c(names(sgp_object@Data), 'SCALE_SCORE_EQUATED_EQUIPERCENTILE_OLD_TO_NEW'), with=FALSE]
            setnames(data.for.equate, 'SCALE_SCORE_EQUATED_EQUIPERCENTILE_OLD_TO_NEW', 'SCALE_SCORE_EQUATED')
            sgp_object@Data <- data.for.equate
        } ### END if (is.null(sgp.use.my.coefficient.matrices))

		equate.variable <- "SCALE_SCORE_EQUATED"
		equate.label <- coefficient.matrix.type <- "EQUATED"
		sgp.percentiles.equated <- TRUE
		sgp.projections.equated <- list(State=state, Year=year.for.equate, Linkages=sgp_object@SGP[['Linkages']])
		tmp_sgp_object <- list(Coefficient_Matrices=sgp_object@SGP[["Coefficient_Matrices"]], Knots_Boundaries=sgp_object@SGP[["Knots_Boundaries"]], Linkages=sgp_object@SGP[['Linkages']])
	} else {
		sgp.percentiles.equated <- FALSE
		equate.variable <- equate.label <- year.for.equate <- sgp.projections.equated <- coefficient.matrix.type <- NULL
		tmp_sgp_object <- list(Coefficient_Matrices=sgp_object@SGP[["Coefficient_Matrices"]], Knots_Boundaries=sgp_object@SGP[["Knots_Boundaries"]])
	} ### END if (sgp.percentiles.equated)

	variables.to.get <- c("VALID_CASE", "YEAR", "CONTENT_AREA", "GRADE", "ID", "SCALE_SCORE", "ACHIEVEMENT_LEVEL", "YEAR_WITHIN", "FIRST_OBSERVATION", "LAST_OBSERVATION",
				"STATE", csem.variable, equate.variable, SGPt)

	if (is.null(SGPt) && is.null(sgp.sqlite) && as.numeric(strsplit(format(object.size(sgp_object@Data), units="GB"), " Gb")[[1]]) > 1) sgp.sqlite <- TRUE else sgp.sqlite <- FALSE
	if (toupper(sgp.sqlite)=="KEEP") {keep.sqlite <- TRUE; sgp.sqlite <- TRUE} else keep.sqlite <- FALSE

	if (sgp.sqlite) {
		del.dir <- dir.create("Data/tmp_data", recursive=TRUE, showWarnings=FALSE)
		tmp_sgp_data_for_analysis <- dbConnect(SQLite(), dbname = "Data/tmp_data/TMP_SGP_Data.sqlite")
		dbWriteTable(tmp_sgp_data_for_analysis, name = "sgp_data", overwrite = TRUE,
			value=sgp_object@Data[,intersect(names(sgp_object@Data), variables.to.get), with=FALSE]["VALID_CASE"], row.names=0)
		sgp.data.names <- dbListFields(tmp_sgp_data_for_analysis, "sgp_data")
		dbDisconnect(tmp_sgp_data_for_analysis)
	} else {
		tmp_sgp_data_for_analysis <- sgp_object@Data[,intersect(names(sgp_object@Data), variables.to.get), with=FALSE]["VALID_CASE"]
		sgp.data.names <- names(tmp_sgp_data_for_analysis)
		if ("YEAR_WITHIN" %in% sgp.data.names) {
			setkey(tmp_sgp_data_for_analysis, VALID_CASE, CONTENT_AREA, YEAR, GRADE, YEAR_WITHIN)
		} else {
			setkey(tmp_sgp_data_for_analysis, VALID_CASE, CONTENT_AREA, YEAR, GRADE)
		}
	}


	#######################################################################################################################
	##   Baseline SGP - compute matrices first if they are not in SGPstateData or merge them into sgp_object if in SGPstateData
	#######################################################################################################################

	if (sgp.percentiles.baseline | sgp.projections.baseline | sgp.projections.lagged.baseline) {

		if (is.null(sgp.use.my.sgp_object.baseline.coefficient.matrices) && is.null(SGPstateData[[state]][["Baseline_splineMatrix"]])) {
			if (is.null(sgp.baseline.config)) {
				sgp.baseline.config <- getSGPBaselineConfig(sgp_object, content_areas, grades, sgp.baseline.panel.years, sgp.percentiles.baseline.max.order)
			} else {
				sgp.baseline.config <- checkConfig(sgp.baseline.config, "Baseline")
			}

			messageSGP("\n\tStarted Baseline Coefficient Matrix Calculation:\n")

			if (!is.null(parallel.config)) { ### PARALLEL BASELINE COEFFICIENT MATRIX CONSTRUCTION
				par.start <- startParallel(parallel.config, 'BASELINE_MATRICES')

				###  FOREACH flavor
				if (toupper(parallel.config[["BACKEND"]]) == "FOREACH") {
					tmp <- foreach(sgp.iter=iter(sgp.baseline.config), .packages="SGP", .errorhandling = "pass", .inorder=FALSE,
						.options.multicore=par.start$foreach.options, .options.mpi=par.start$foreach.options, .options.redis=par.start$foreach.options) %dopar% {
						return(baselineSGP(
							sgp_object,
							state=state,
							sgp.baseline.config=list(sgp.iter), ## NOTE: list of sgp.iter must be passed for proper iteration
							return.matrices.only=TRUE,
							calculate.baseline.sgps=FALSE))
					}
					tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
					if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
						tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
							sgp.percentiles.baseline.=getErrorReports(tmp, tmp.tf, rev(par.sgp.config[['sgp.percentiles.baseline']])))
					}
				} else {
					if (par.start$par.type=="SNOW") {
						tmp <- clusterApplyLB(par.start$internal.cl, sgp.baseline.config, function(sgp.iter) baselineSGP(
								sgp_object,
								state=state,
								sgp.baseline.config=list(sgp.iter), ## NOTE: list of sgp.iter must be passed for proper iteration
								return.matrices.only=TRUE,
								calculate.baseline.sgps=FALSE))

						tmp_sgp_object <- mergeSGP(tmp_sgp_object, list(Coefficient_Matrices=merge.coefficient.matrices(tmp)))
					} # END if (SNOW)

					if (par.start$par.type=="MULTICORE") {
						tmp <- mclapply(sgp.baseline.config, function(sgp.iter) baselineSGP(
									sgp_object,
									state=state,
									sgp.baseline.config=list(sgp.iter), ## NOTE: list of sgp.iter must be passed for proper iteration
									return.matrices.only=TRUE,
									calculate.baseline.sgps=FALSE),
								mc.cores=par.start$workers, mc.preschedule=FALSE)
						tmp_sgp_object <- mergeSGP(tmp_sgp_object, list(Coefficient_Matrices=merge.coefficient.matrices(tmp)))
					} # END if (MULTICORE)
					stopParallel(parallel.config, par.start)
				} #  END  if parallel
			} else {
				##  SEQUENTIAL BASELINE COEFFICIENT MATRIX CONSTRUCTION
				##  Or, run TAUS in parallel in studentGrowthPercentiles using lower.level.parallel.config
				##  Useful if many more cores/workers available than configs to iterate over.
				tmp <- list()
				for (sgp.iter in seq_along(sgp.baseline.config)) {
					tmp[[sgp.iter]] <- baselineSGP(
						sgp_object,
						state=state,
						sgp.baseline.config=sgp.baseline.config[sgp.iter], ## NOTE: must pass list, [...], not vector, [[...]].
						return.matrices.only=TRUE,
						calculate.baseline.sgps=FALSE,
						parallel.config=lower.level.parallel.config)
				}
				tmp_sgp_object <- mergeSGP(tmp_sgp_object, list(Coefficient_Matrices=merge.coefficient.matrices(tmp)))
			}
			rm(tmp)

			assign(paste(state, "_Baseline_Matrices", sep=""), list())
			for (tmp.matrix.label in grep("BASELINE", names(tmp_sgp_object$Coefficient_Matrices), value=TRUE)) {
				eval(parse(text=paste(state, "_Baseline_Matrices[['", tmp.matrix.label, "']] <- tmp_sgp_object[['Coefficient_Matrices']][['", tmp.matrix.label, "']]", sep="")))
			}
			save(list=paste(state, "_Baseline_Matrices", sep=""), file=paste(state, "_Baseline_Matrices.Rdata", sep=""))
			messageSGP("\n\tFinished Calculating Baseline Coefficient Matrices\n")
		} else {
			if (is.null(sgp.use.my.sgp_object.baseline.coefficient.matrices)) tmp_sgp_object <- mergeSGP(tmp_sgp_object, SGPstateData[[state]][["Baseline_splineMatrix"]])
		}
	} # END Get/Compute baseline coefficient matrices


	#######################################################################################################################
	##   SIMEX Baseline SGP - compute matrices first if they are not in sgp_object (NOTE: Not stored in SGPstateData)
	#######################################################################################################################

	if (sgp.percentiles.baseline & !is.null(calculate.simex.baseline)) {

		if (!is.null(sgp.config)) {
			tmp.subjects <- unique(sapply(sgp.config, function(x) tail(x[["sgp.content.areas"]],1)))
		} else {
			if (!is.null(content_areas)) tmp.subjects <- content_areas else tmp.subjects <- unique(sgp_object@Data["VALID_CASE"][["CONTENT_AREA"]])
		}

		###  Calculate BASELINE SIMEX matrices if they are not present
		if (!all(find.matrices <- paste(tmp.subjects, ".BASELINE.SIMEX", sep="") %in% names(tmp_sgp_object[["Coefficient_Matrices"]]))) {

			if (is.null(sgp.baseline.config)) {
				sgp.baseline.config <- getSGPBaselineConfig(sgp_object, content_areas = tmp.subjects, grades, sgp.baseline.panel.years, sgp.percentiles.baseline.max.order, calculate.simex.baseline)
			} else {
				sgp.baseline.config <- checkConfig(sgp.baseline.config, "Baseline")
			}

			sgp.baseline.config <- sgp.baseline.config[which(sapply(sgp.baseline.config, function(x) tail(x[["sgp.baseline.content.areas"]],1)) %in% tmp.subjects[!find.matrices])]

			messageSGP("\n\tStarted SIMEX Baseline Coefficient Matrix Calculation:\n")

			##  Enforce that simex.use.my.coefficient.matrices must be FALSE for BASELINE SIMEX matrix production
			calculate.simex.baseline$simex.use.my.coefficient.matrices <- NULL

			if (!is.null(parallel.config)) { ### PARALLEL BASELINE COEFFICIENT MATRIX CONSTRUCTION

				par.start <- startParallel(parallel.config, 'BASELINE_MATRICES')

				##  FOREACH flavor
				if (toupper(parallel.config[["BACKEND"]]) == "FOREACH") {
					tmp <- foreach(sgp.iter=iter(sgp.baseline.config), .packages="SGP", .errorhandling = "pass", .inorder=FALSE,
						.options.multicore=par.start$foreach.options, .options.mpi=par.start$foreach.options, .options.redis=par.start$foreach.options) %dopar% {
						return(baselineSGP(
							sgp_object,
							state=state,
							sgp.baseline.config=list(sgp.iter), ## NOTE: list of sgp.iter must be passed for proper iteration
							return.matrices.only=TRUE,
							calculate.baseline.sgps=FALSE,
							calculate.simex.baseline=calculate.simex.baseline,
							parallel.config=parallel.config))
					}
					tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
					if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
						tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
							sgp.percentiles.baseline.=getErrorReports(tmp, tmp.tf, rev(par.sgp.config[['sgp.percentiles.baseline']])))
					}
				} else {  ## SNOW and MULTICORE flavors
					if (par.start$par.type=="SNOW") {
						tmp <- clusterApplyLB(par.start$internal.cl, sgp.baseline.config, function(sgp.iter) baselineSGP(
							sgp_object,
							state=state,
							sgp.baseline.config=list(sgp.iter), ## NOTE: list of sgp.iter must be passed for proper iteration
							return.matrices.only=TRUE,
							calculate.baseline.sgps=FALSE,
							calculate.simex.baseline=calculate.simex.baseline,
							parallel.config=parallel.config))

						tmp_sgp_object <- mergeSGP(tmp_sgp_object, list(Coefficient_Matrices=merge.coefficient.matrices(tmp, simex=TRUE)))
					} # END if (SNOW)

					if (par.start$par.type=="MULTICORE") {
						tmp <- mclapply(sgp.baseline.config, function(sgp.iter) baselineSGP(
							sgp_object,
							state=state,
							sgp.baseline.config=list(sgp.iter), ## NOTE: list of sgp.iter must be passed for proper iteration
							return.matrices.only=TRUE,
							calculate.baseline.sgps=FALSE,
							calculate.simex.baseline=calculate.simex.baseline,
							parallel.config=parallel.config),
							mc.cores=par.start$workers, mc.preschedule=FALSE)

						tmp_sgp_object <- mergeSGP(tmp_sgp_object, list(Coefficient_Matrices=merge.coefficient.matrices(tmp, simex=TRUE)))
					} # END if (MULTICORE)
					stopParallel(parallel.config, par.start)
				} #  END FOREACH, SNOW and MULTICORE
			} else {
				## SEQUENTIAL BASELINE COEFFICIENT MATRIX CONSTRUCTION
				##  Or, run SIMEX simulation iterations in parallel in studentGrowthPercentiles using lower.level.parallel.config
				##  Useful if many more cores/workers available than configs to iterate over.
				tmp <- list()
				for (sgp.iter in seq_along(sgp.baseline.config)) {
					tmp[[sgp.iter]] <- baselineSGP(
						sgp_object,
						state=state,
						sgp.baseline.config=sgp.baseline.config[sgp.iter], ## NOTE: must pass list, [...], not vector, [[...]].
						return.matrices.only=TRUE,
						calculate.baseline.sgps=FALSE,
						calculate.simex.baseline=calculate.simex.baseline,
						parallel.config=lower.level.parallel.config)
				}

				tmp_sgp_object <- mergeSGP(tmp_sgp_object, list(Coefficient_Matrices=merge.coefficient.matrices(tmp, simex=TRUE)))
			}

			###  Save SIMEX BASELINE matrices
			assign(paste(state, "_SIMEX_Baseline_Matrices", sep=""), list())
			for (tmp.matrix.label in grep("BASELINE.SIMEX", names(tmp_sgp_object$Coefficient_Matrices), value=TRUE)) {
				eval(parse(text=paste(state, "_SIMEX_Baseline_Matrices[['", tmp.matrix.label, "']] <- tmp_sgp_object[['Coefficient_Matrices']][['", tmp.matrix.label, "']]", sep="")))
			}
			save(list=paste(state, "_SIMEX_Baseline_Matrices", sep=""), file=paste(state, "_SIMEX_Baseline_Matrices.Rdata", sep=""), compress="xz")

			messageSGP("\n\tFinished Calculating SIMEX Baseline Coefficient Matrices\n")
		} # END Compute SIMEX baseline coefficient matrices

		##  Enforce that simex.use.my.coefficient.matrices must be TRUE and save.matrices is FALSE for BASELINE SIMEX calculations below
		calculate.simex.baseline$simex.use.my.coefficient.matrices <- TRUE
		calculate.simex.baseline$save.matrices <- FALSE

	} # END check for SIMEX baseline matrices presence


	### Create par.sgp.config (for both parallel and sequential implementations) and par.sgp.config.projections for projections

	setkeyv(sgp_object@Data, getKey(sgp_object))
	par.sgp.config <- getSGPConfig(sgp_object, state, tmp_sgp_object, content_areas, years, grades, sgp.config, trim.sgp.config, sgp.percentiles, sgp.projections, sgp.projections.lagged,
		sgp.percentiles.baseline, sgp.projections.baseline, sgp.projections.lagged.baseline, sgp.config.drop.nonsequential.grade.progression.variables,
		sgp.minimum.default.panel.years, sgp.projections.max.forward.progression.years, sgp.use.my.coefficient.matrices, calculate.simex, calculate.simex.baseline, year.for.equate,
        sgp.percentiles.equated, SGPt)

	if (sgp.projections & length(par.sgp.config[['sgp.projections']])==0) {
		messageSGP("\tNOTE: No configurations are present for cohort referenced projections. No cohort referenced projections will be calculated.")
		sgp.projections <- FALSE
	}

	if (sgp.projections.lagged & length(par.sgp.config[['sgp.projections.lagged']])==0) {
		messageSGP("\tNOTE: No configurations are present for cohort referenced lagged projections. No lagged cohort referenced projections will be calculated.")
		sgp.projections.lagged <- FALSE
	}

	if (sgp.projections.baseline & length(par.sgp.config[['sgp.projections.baseline']])==0) {
		messageSGP("\tNOTE: No configurations are present for baseline projections. No baseline projections will be calculated.")
		sgp.projections.baseline <- sgp.projections.lagged.baseline <- FALSE
	}

	if (sgp.projections.lagged.baseline & length(par.sgp.config[['sgp.projections.lagged.baseline']])==0) {
		messageSGP("\tNOTE: No configurations are present for lagged baseline projections. No lagged baseline projections will be calculated.")
		sgp.projections.baseline <- sgp.projections.lagged.baseline <- FALSE
	}

	### Produce cohort data information

	if (get.cohort.data.info) {
		cohort_data_info <- getCohortDataInfo(tmp_sgp_data_for_analysis, par.sgp.config[['sgp.percentiles']])
		save(cohort_data_info, file="cohort_data_info.Rdata")
	}


	#######################################################################################################################
	#######################################################################################################################
	##   Percentiles, Baseline Percentiles, Projections, Lagged Projections -  PARALLEL FLAVORS FIRST
	#######################################################################################################################
	#######################################################################################################################

	if (!is.null(parallel.config)) {

	##################################
	###  PERCENTILES
	##################################

		if (sgp.percentiles) {
			par.start <- startParallel(parallel.config, 'PERCENTILES')

			###  FOREACH flavor
			if (toupper(parallel.config[["BACKEND"]]) == "FOREACH") {
				tmp <- foreach(sgp.iter=iter(rev(par.sgp.config[['sgp.percentiles']])), .packages="SGP", .errorhandling = "pass", .inorder=FALSE,
					.options.multicore=par.start$foreach.options, .options.mpi=par.start$foreach.options, .options.redis=par.start$foreach.options) %dopar% {
					return(studentGrowthPercentiles(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.percentiles", sgp.iter, csem.variable, equate.variable, SGPt=SGPt),
							Coefficient_Matrices=sgp.iter[["sgp.matrices"]],
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.percentiles")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						use.my.coefficient.matrices=sgp.use.my.coefficient.matrices,
						growth.levels=state,
						calculate.confidence.intervals=get.simulate.sgps.arg(calculate.confidence.intervals, state, sgp.iter),
						panel.data.vnames=getPanelDataVnames("sgp.percentiles", sgp.iter, sgp.data.names),
						additional.vnames.to.return=getPanelDataVnames("sgp.percentiles.to.return", sgp.iter, sgp.data.names),
						grade.progression=sgp.iter[["sgp.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.content.areas"]],
						year.progression=sgp.iter[["sgp.panel.years"]],
						max.order.for.percentile=SGPstateData[[state]][["SGP_Configuration"]][["max.order.for.percentile"]],
						percentile.cuts=SGPstateData[[state]][["SGP_Configuration"]][["percentile.cuts"]],
						drop.nonsequential.grade.progression.variables=FALSE, # taken care of with config
						exact.grade.progression.sequence=sgp.iter[["sgp.exact.grade.progression"]],
						sgp.loss.hoss.adjustment=sgp.loss.hoss.adjustment,
						sgp.cohort.size=tmp.cohort.size,
						return.norm.group.scale.scores=return.norm.group.scale.scores,
						return.norm.group.dates=return.norm.group.dates,
						return.prior.scale.score.standardized=return.prior.scale.score.standardized,
						goodness.of.fit=goodness.of.fit.print.arg,
						goodness.of.fit.minimum.n=SGPstateData[[state]][["SGP_Configuration"]][["goodness.of.fit.minimum.n"]],
						verbose.output=verbose.output,
						print.other.gp=print.other.gp,
						calculate.simex=get.calculate.simex.arg(sgp.iter[["sgp.calculate.simex"]], sgp.iter),
						max.n.for.coefficient.matrices=max.n.for.coefficient.matrices,
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.percentiles"),
                        SGPt.max.time=SGPt.max.time,
                        calculate.sgps=sgp.percentiles.calculate.sgps,
						parallel.config=par.start$Lower_Level_Parallel,
						...))
					}
				tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
				if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
					tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
						sgp.percentiles.=getErrorReports(tmp, tmp.tf, rev(par.sgp.config[['sgp.percentiles']])))
				}
			} else { # END FOREACH
				###    SNOW flavor
				if (par.start$par.type == 'SNOW') {
					tmp <- clusterApplyLB(par.start$internal.cl, rev(par.sgp.config[['sgp.percentiles']]), function(sgp.iter) studentGrowthPercentiles(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.percentiles", sgp.iter, csem.variable, equate.variable, SGPt=SGPt),
							Coefficient_Matrices=sgp.iter[["sgp.matrices"]],
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.percentiles")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						use.my.coefficient.matrices=sgp.use.my.coefficient.matrices,
						growth.levels=state,
						calculate.confidence.intervals=get.simulate.sgps.arg(calculate.confidence.intervals, state, sgp.iter),
						panel.data.vnames=getPanelDataVnames("sgp.percentiles", sgp.iter, sgp.data.names),
						additional.vnames.to.return=getPanelDataVnames("sgp.percentiles.to.return", sgp.iter, sgp.data.names),
						grade.progression=sgp.iter[["sgp.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.content.areas"]],
						year.progression=sgp.iter[["sgp.panel.years"]],
						max.order.for.percentile=SGPstateData[[state]][["SGP_Configuration"]][["max.order.for.percentile"]],
						percentile.cuts=SGPstateData[[state]][["SGP_Configuration"]][["percentile.cuts"]],
						drop.nonsequential.grade.progression.variables=FALSE, # taken care of with config
						exact.grade.progression.sequence=sgp.iter[["sgp.exact.grade.progression"]],
						sgp.loss.hoss.adjustment=sgp.loss.hoss.adjustment,
						sgp.cohort.size=tmp.cohort.size,
						return.norm.group.scale.scores=return.norm.group.scale.scores,
						return.norm.group.dates=return.norm.group.dates,
						return.prior.scale.score.standardized=return.prior.scale.score.standardized,
						goodness.of.fit=goodness.of.fit.print.arg,
						goodness.of.fit.minimum.n=SGPstateData[[state]][["SGP_Configuration"]][["goodness.of.fit.minimum.n"]],
						verbose.output=verbose.output,
						print.other.gp=print.other.gp,
						calculate.simex=get.calculate.simex.arg(sgp.iter[["sgp.calculate.simex"]], sgp.iter),
						max.n.for.coefficient.matrices=max.n.for.coefficient.matrices,
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.percentiles"),
                        SGPt.max.time=SGPt.max.time,
                        calculate.sgps=sgp.percentiles.calculate.sgps,
						parallel.config=par.start$Lower_Level_Parallel,
						...))

					tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
					if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
						tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
							sgp.percentiles.=getErrorReports(tmp, tmp.tf, rev(par.sgp.config[['sgp.percentiles']])))
					}
				} # END SNOW

				###  MULTICORE flavor
				if (par.start$par.type == 'MULTICORE') {
					tmp <- mclapply(rev(par.sgp.config[['sgp.percentiles']]), function(sgp.iter) studentGrowthPercentiles(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.percentiles", sgp.iter, csem.variable, equate.variable, SGPt=SGPt),
							Coefficient_Matrices=sgp.iter[["sgp.matrices"]],
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.percentiles")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						use.my.coefficient.matrices=sgp.use.my.coefficient.matrices,
						growth.levels=state,
						calculate.confidence.intervals=get.simulate.sgps.arg(calculate.confidence.intervals, state, sgp.iter),
						panel.data.vnames=getPanelDataVnames("sgp.percentiles", sgp.iter, sgp.data.names),
						additional.vnames.to.return=getPanelDataVnames("sgp.percentiles.to.return", sgp.iter, sgp.data.names),
						grade.progression=sgp.iter[["sgp.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.content.areas"]],
						year.progression=sgp.iter[["sgp.panel.years"]],
						max.order.for.percentile=SGPstateData[[state]][["SGP_Configuration"]][["max.order.for.percentile"]],
						percentile.cuts=SGPstateData[[state]][["SGP_Configuration"]][["percentile.cuts"]],
						drop.nonsequential.grade.progression.variables=FALSE, # taken care of with config
						exact.grade.progression.sequence=sgp.iter[["sgp.exact.grade.progression"]],
						sgp.loss.hoss.adjustment=sgp.loss.hoss.adjustment,
						sgp.cohort.size=tmp.cohort.size,
						return.norm.group.scale.scores=return.norm.group.scale.scores,
						return.norm.group.dates=return.norm.group.dates,
						return.prior.scale.score.standardized=return.prior.scale.score.standardized,
						goodness.of.fit=goodness.of.fit.print.arg,
						goodness.of.fit.minimum.n=SGPstateData[[state]][["SGP_Configuration"]][["goodness.of.fit.minimum.n"]],
						verbose.output=verbose.output,
						print.other.gp=print.other.gp,
						calculate.simex=get.calculate.simex.arg(sgp.iter[["sgp.calculate.simex"]], sgp.iter),
						max.n.for.coefficient.matrices=max.n.for.coefficient.matrices,
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.percentiles"),
                        SGPt.max.time=SGPt.max.time,
                        calculate.sgps=sgp.percentiles.calculate.sgps,
						parallel.config=par.start$Lower_Level_Parallel,
						...), mc.cores=par.start$workers, mc.preschedule=FALSE)

					tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
					if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
						tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
							sgp.percentiles.=getErrorReports(tmp, tmp.tf, rev(par.sgp.config[['sgp.percentiles']])))
					}
				} # End MULTICORE
			} # #END not FOREACH
			stopParallel(parallel.config, par.start)
		} #END if (sgp.percentiles)


	##################################
	###  PERCENTILES EQUATED
	##################################

		if (sgp.percentiles.equated) {
			par.start <- startParallel(parallel.config, 'PERCENTILES')

			###  FOREACH flavor
			if (toupper(parallel.config[["BACKEND"]]) == "FOREACH") {
				tmp <- foreach(sgp.iter=iter(rev(par.sgp.config[['sgp.percentiles.equated']])), .packages="SGP", .errorhandling = "pass", .inorder=FALSE,
					.options.multicore=par.start$foreach.options, .options.mpi=par.start$foreach.options, .options.redis=par.start$foreach.options) %dopar% {
					return(studentGrowthPercentiles(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.percentiles", sgp.iter, csem.variable, equate.variable, SGPt=SGPt),
							Coefficient_Matrices=sgp.iter[["sgp.equated.matrices"]],
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.percentiles")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1),
								my.subject=tail(sgp.iter[["sgp.content.areas"]], 1), my.extra.label=equate.label),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						use.my.coefficient.matrices=sgp.use.my.coefficient.matrices,
						growth.levels=state,
						calculate.confidence.intervals=get.simulate.sgps.arg(calculate.confidence.intervals, state, sgp.iter),
						panel.data.vnames=getPanelDataVnames("sgp.percentiles", sgp.iter, sgp.data.names, equate.variable),
						additional.vnames.to.return=getPanelDataVnames("sgp.percentiles.to.return", sgp.iter, sgp.data.names),
						grade.progression=sgp.iter[["sgp.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.content.areas"]],
						year.progression=sgp.iter[["sgp.panel.years"]],
						max.order.for.percentile=SGPstateData[[state]][["SGP_Configuration"]][["max.order.for.percentile"]],
						percentile.cuts=SGPstateData[[state]][["SGP_Configuration"]][["percentile.cuts"]],
						drop.nonsequential.grade.progression.variables=FALSE, # taken care of with config
						exact.grade.progression.sequence=sgp.iter[["sgp.exact.grade.progression"]],
						sgp.loss.hoss.adjustment=sgp.loss.hoss.adjustment,
						sgp.cohort.size=tmp.cohort.size,
						return.norm.group.scale.scores=return.norm.group.scale.scores,
						return.norm.group.dates=return.norm.group.dates,
						return.prior.scale.score.standardized=return.prior.scale.score.standardized,
						goodness.of.fit=goodness.of.fit.print.arg,
						goodness.of.fit.minimum.n=SGPstateData[[state]][["SGP_Configuration"]][["goodness.of.fit.minimum.n"]],
						verbose.output=verbose.output,
						print.other.gp=print.other.gp,
						max.n.for.coefficient.matrices=max.n.for.coefficient.matrices,
						sgp.percentiles.equated=sgp.projections.equated,
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.percentiles.equated"),
                        SGPt.max.time=SGPt.max.time,
						...))
					}
				tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
				if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
					tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
						sgp.percentiles.equated.=getErrorReports(tmp, tmp.tf, rev(par.sgp.config[['sgp.percentiles.equated']])))
				}
			} else { # END FOREACH
				###    SNOW flavor
				if (par.start$par.type == 'SNOW') {
					tmp <- clusterApplyLB(par.start$internal.cl, rev(par.sgp.config[['sgp.percentiles.equated']]), function(sgp.iter) studentGrowthPercentiles(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.percentiles", sgp.iter, csem.variable, equate.variable, SGPt=SGPt),
							Coefficient_Matrices=sgp.iter[["sgp.equated.matrices"]],
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.percentiles")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1),
								my.subject=tail(sgp.iter[["sgp.content.areas"]], 1), my.extra.label=equate.label),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						use.my.coefficient.matrices=sgp.use.my.coefficient.matrices,
						growth.levels=state,
						calculate.confidence.intervals=get.simulate.sgps.arg(calculate.confidence.intervals, state, sgp.iter),
						panel.data.vnames=getPanelDataVnames("sgp.percentiles", sgp.iter, sgp.data.names, equate.variable),
						additional.vnames.to.return=getPanelDataVnames("sgp.percentiles.to.return", sgp.iter, sgp.data.names),
						grade.progression=sgp.iter[["sgp.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.content.areas"]],
						year.progression=sgp.iter[["sgp.panel.years"]],
						max.order.for.percentile=SGPstateData[[state]][["SGP_Configuration"]][["max.order.for.percentile"]],
						percentile.cuts=SGPstateData[[state]][["SGP_Configuration"]][["percentile.cuts"]],
						drop.nonsequential.grade.progression.variables=FALSE, # taken care of with config
						exact.grade.progression.sequence=sgp.iter[["sgp.exact.grade.progression"]],
						sgp.loss.hoss.adjustment=sgp.loss.hoss.adjustment,
						sgp.cohort.size=tmp.cohort.size,
						return.norm.group.scale.scores=return.norm.group.scale.scores,
						return.norm.group.dates=return.norm.group.dates,
						return.prior.scale.score.standardized=return.prior.scale.score.standardized,
						goodness.of.fit=goodness.of.fit.print.arg,
						goodness.of.fit.minimum.n=SGPstateData[[state]][["SGP_Configuration"]][["goodness.of.fit.minimum.n"]],
						verbose.output=verbose.output,
						print.other.gp=print.other.gp,
						max.n.for.coefficient.matrices=max.n.for.coefficient.matrices,
						sgp.percentiles.equated=sgp.projections.equated,
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.percentiles.equated"),
                        SGPt.max.time=SGPt.max.time,
						...))

					tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
					if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
						tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
							sgp.percentiles.equated.=getErrorReports(tmp, tmp.tf, rev(par.sgp.config[['sgp.percentiles.equated']])))
					}
				} # END SNOW

				###  MULTICORE flavor
				if (par.start$par.type == 'MULTICORE') {
					tmp <- mclapply(rev(par.sgp.config[['sgp.percentiles.equated']]), function(sgp.iter) studentGrowthPercentiles(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.percentiles", sgp.iter, csem.variable, equate.variable, SGPt=SGPt),
							Coefficient_Matrices=sgp.iter[["sgp.equated.matrices"]],
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.percentiles")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1),
								my.subject=tail(sgp.iter[["sgp.content.areas"]], 1), my.extra.label=equate.label),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						use.my.coefficient.matrices=sgp.use.my.coefficient.matrices,
						growth.levels=state,
						calculate.confidence.intervals=get.simulate.sgps.arg(calculate.confidence.intervals, state, sgp.iter),
						panel.data.vnames=getPanelDataVnames("sgp.percentiles", sgp.iter, sgp.data.names, equate.variable),
						additional.vnames.to.return=getPanelDataVnames("sgp.percentiles.to.return", sgp.iter, sgp.data.names),
						grade.progression=sgp.iter[["sgp.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.content.areas"]],
						year.progression=sgp.iter[["sgp.panel.years"]],
						max.order.for.percentile=SGPstateData[[state]][["SGP_Configuration"]][["max.order.for.percentile"]],
						percentile.cuts=SGPstateData[[state]][["SGP_Configuration"]][["percentile.cuts"]],
						drop.nonsequential.grade.progression.variables=FALSE, # taken care of with config
						exact.grade.progression.sequence=sgp.iter[["sgp.exact.grade.progression"]],
						sgp.loss.hoss.adjustment=sgp.loss.hoss.adjustment,
						sgp.cohort.size=tmp.cohort.size,
						return.norm.group.scale.scores=return.norm.group.scale.scores,
						return.norm.group.dates=return.norm.group.dates,
						return.prior.scale.score.standardized=return.prior.scale.score.standardized,
						goodness.of.fit=goodness.of.fit.print.arg,
						goodness.of.fit.minimum.n=SGPstateData[[state]][["SGP_Configuration"]][["goodness.of.fit.minimum.n"]],
						verbose.output=verbose.output,
						print.other.gp=print.other.gp,
						max.n.for.coefficient.matrices=max.n.for.coefficient.matrices,
						sgp.percentiles.equated=sgp.projections.equated,
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.percentiles.equated"),
                        SGPt.max.time=SGPt.max.time,
						parallel.config=par.start$Lower_Level_Parallel,
						...), mc.cores=par.start$workers, mc.preschedule=FALSE)

					tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
					if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
						tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
							sgp.percentiles.equated.=getErrorReports(tmp, tmp.tf, rev(par.sgp.config[['sgp.percentiles.equated']])))
					}
				} # End MULTICORE
			} # #END not FOREACH
			stopParallel(parallel.config, par.start)
		} #END if (sgp.percentiles)


	####################################
	###  BASELINE PERCENTILES
	####################################

		if (sgp.percentiles.baseline & length(par.sgp.config[["sgp.percentiles.baseline"]]) > 0) {

			par.start <- startParallel(parallel.config, 'BASELINE_PERCENTILES')

			###  FOREACH flavor
			if (toupper(parallel.config[["BACKEND"]]) == "FOREACH") {
				tmp <- foreach(sgp.iter=iter(rev(par.sgp.config[['sgp.percentiles.baseline']])), .packages="SGP", .errorhandling = "pass", .inorder=FALSE,
					.options.multicore=par.start$foreach.options, .options.mpi=par.start$foreach.options, .options.redis=par.start$foreach.options) %dopar% {
					return(studentGrowthPercentiles(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.percentiles", sgp.iter, csem.variable, SGPt=SGPt),
							Coefficient_Matrices=sgp.iter[["sgp.baseline.matrices"]],
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.percentiles")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1),
							my.subject=tail(sgp.iter[["sgp.content.areas"]], 1), my.extra.label="BASELINE"),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						use.my.coefficient.matrices=list(my.year="BASELINE", my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						growth.levels=state,
						calculate.confidence.intervals=get.simulate.sgps.arg(calculate.confidence.intervals, state, sgp.iter),
						panel.data.vnames=getPanelDataVnames("sgp.percentiles.baseline", sgp.iter, sgp.data.names),
						additional.vnames.to.return=getPanelDataVnames("sgp.percentiles.to.return", sgp.iter, sgp.data.names),
						grade.progression=sgp.iter[["sgp.baseline.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.baseline.content.areas"]],
						year_lags.progression=sgp.iter[["sgp.baseline.panel.years.lags"]],
						num.prior=min(sgp.iter[["sgp.baseline.max.order"]], sgp.percentiles.baseline.max.order),
						percentile.cuts=SGPstateData[[state]][["SGP_Configuration"]][["percentile.cuts"]],
						drop.nonsequential.grade.progression.variables=FALSE, # taken care of with config
						exact.grade.progression.sequence=sgp.iter[["sgp.exact.grade.progression"]],
						sgp.loss.hoss.adjustment=sgp.loss.hoss.adjustment,
						return.norm.group.scale.scores=return.norm.group.scale.scores,
						return.norm.group.dates=return.norm.group.dates,
						return.prior.scale.score.standardized=return.prior.scale.score.standardized,
						goodness.of.fit=goodness.of.fit.print.arg,
						goodness.of.fit.minimum.n=SGPstateData[[state]][["SGP_Configuration"]][["goodness.of.fit.minimum.n"]],
						verbose.output=verbose.output,
						print.other.gp=print.other.gp,
						calculate.simex=get.calculate.simex.arg(sgp.iter[["sgp.calculate.simex.baseline"]], sgp.iter),
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.percentiles.baseline"),
                        SGPt.max.time=SGPt.max.time,
						...))
				}
				tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
				if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
					tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
						sgp.percentiles.baseline.=getErrorReports(tmp, tmp.tf, rev(par.sgp.config[['sgp.percentiles.baseline']])))
				}
			} else { # END FOREACH
				###    SNOW flavor
				if (par.start$par.type == 'SNOW') {
					tmp <- clusterApplyLB(par.start$internal.cl, rev(par.sgp.config[['sgp.percentiles.baseline']]), function(sgp.iter) studentGrowthPercentiles(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.percentiles", sgp.iter, csem.variable, SGPt=SGPt),
							Coefficient_Matrices=sgp.iter[["sgp.baseline.matrices"]],
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.percentiles")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1),
							my.subject=tail(sgp.iter[["sgp.content.areas"]], 1), my.extra.label="BASELINE"),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						use.my.coefficient.matrices=list(my.year="BASELINE", my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						growth.levels=state,
						calculate.confidence.intervals=get.simulate.sgps.arg(calculate.confidence.intervals, state, sgp.iter),
						panel.data.vnames=getPanelDataVnames("sgp.percentiles.baseline", sgp.iter, sgp.data.names),
						additional.vnames.to.return=getPanelDataVnames("sgp.percentiles.to.return", sgp.iter, sgp.data.names),
						grade.progression=sgp.iter[["sgp.baseline.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.baseline.content.areas"]],
						year_lags.progression=sgp.iter[["sgp.baseline.panel.years.lags"]],
						num.prior=min(sgp.iter[["sgp.baseline.max.order"]], sgp.percentiles.baseline.max.order),
						percentile.cuts=SGPstateData[[state]][["SGP_Configuration"]][["percentile.cuts"]],
						drop.nonsequential.grade.progression.variables=FALSE, # taken care of with config
						exact.grade.progression.sequence=sgp.iter[["sgp.exact.grade.progression"]],
						sgp.loss.hoss.adjustment=sgp.loss.hoss.adjustment,
						return.norm.group.scale.scores=return.norm.group.scale.scores,
						return.norm.group.dates=return.norm.group.dates,
						return.prior.scale.score.standardized=return.prior.scale.score.standardized,
						goodness.of.fit=goodness.of.fit.print.arg,
						goodness.of.fit.minimum.n=SGPstateData[[state]][["SGP_Configuration"]][["goodness.of.fit.minimum.n"]],
						verbose.output=verbose.output,
						print.other.gp=print.other.gp,
						calculate.simex=get.calculate.simex.arg(sgp.iter[["sgp.calculate.simex.baseline"]], sgp.iter),
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.percentiles.baseline"),
                        SGPt.max.time=SGPt.max.time,
						...))

					tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
					if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
						tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
							sgp.percentiles.baseline.=getErrorReports(tmp, tmp.tf, rev(par.sgp.config[['sgp.percentiles.baseline']])))
					}
				} # END SNOW

				###  MULTICORE flavor
				if (par.start$par.type == 'MULTICORE') {
					tmp <- mclapply(rev(par.sgp.config[['sgp.percentiles.baseline']]), function(sgp.iter) studentGrowthPercentiles(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.percentiles", sgp.iter, csem.variable, SGPt=SGPt),
							Coefficient_Matrices=sgp.iter[["sgp.baseline.matrices"]],
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.percentiles")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1),
							my.subject=tail(sgp.iter[["sgp.content.areas"]], 1), my.extra.label="BASELINE"),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						use.my.coefficient.matrices=list(my.year="BASELINE", my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						growth.levels=state,
						calculate.confidence.intervals=get.simulate.sgps.arg(calculate.confidence.intervals, state, sgp.iter),
						panel.data.vnames=getPanelDataVnames("sgp.percentiles.baseline", sgp.iter, sgp.data.names),
						additional.vnames.to.return=getPanelDataVnames("sgp.percentiles.to.return", sgp.iter, sgp.data.names),
						grade.progression=sgp.iter[["sgp.baseline.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.baseline.content.areas"]],
						year_lags.progression=sgp.iter[["sgp.baseline.panel.years.lags"]],
						num.prior=min(sgp.iter[["sgp.baseline.max.order"]], sgp.percentiles.baseline.max.order),
						percentile.cuts=SGPstateData[[state]][["SGP_Configuration"]][["percentile.cuts"]],
						drop.nonsequential.grade.progression.variables=FALSE, # taken care of with config
						exact.grade.progression.sequence=sgp.iter[["sgp.exact.grade.progression"]],
						sgp.loss.hoss.adjustment=sgp.loss.hoss.adjustment,
						return.norm.group.scale.scores=return.norm.group.scale.scores,
						return.norm.group.dates=return.norm.group.dates,
						return.prior.scale.score.standardized=return.prior.scale.score.standardized,
						goodness.of.fit=goodness.of.fit.print.arg,
						goodness.of.fit.minimum.n=SGPstateData[[state]][["SGP_Configuration"]][["goodness.of.fit.minimum.n"]],
						verbose.output=verbose.output,
						print.other.gp=print.other.gp,
						calculate.simex=get.calculate.simex.arg(sgp.iter[["sgp.calculate.simex.baseline"]], sgp.iter),
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.percentiles.baseline"),
                        SGPt.max.time=SGPt.max.time,
						parallel.config=par.start$Lower_Level_Parallel,
						...), mc.cores=par.start$workers, mc.preschedule=FALSE)

					tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
					if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
						tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
							sgp.percentiles.baseline.=getErrorReports(tmp, tmp.tf, rev(par.sgp.config[['sgp.percentiles.baseline']])))
					}
				} # End MULTICORE
			} # END parallel flavors
			stopParallel(parallel.config, par.start)
		} ## END if sgp.percentiles.baseline


	#######################################################
	###  PROJECTIONS (COHORT referenced)
	#######################################################

		if (sgp.projections) {

			par.start <- startParallel(parallel.config, 'PROJECTIONS')

			###  FOREACH flavor
			if (toupper(parallel.config[["BACKEND"]]) == "FOREACH") {
				tmp <- foreach(sgp.iter=iter(par.sgp.config[['sgp.projections']]), .packages="SGP", .errorhandling = "pass", .inorder=FALSE,
					.options.multicore=par.start$foreach.options, .options.mpi=par.start$foreach.options, .options.redis=par.start$foreach.options) %dopar% {
					return(studentGrowthProjections(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.projections", sgp.iter, sgp.scale.score.equated=equate.variable, SGPt=SGPt),
							Coefficient_Matrices=selectCoefficientMatrices(tmp_sgp_object, coefficient.matrix.type),
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.projections")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.projection.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.projection.content.areas"]], 1)),
						use.my.coefficient.matrices=list(my.year=tail(sgp.iter[["sgp.projection.panel.years"]], 1),
							my.subject=tail(sgp.iter[["sgp.projection.content.areas"]], 1), my.extra.label=equate.label),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.projection.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.projection.content.areas"]], 1)),
						performance.level.cutscores=state,
						max.forward.progression.grade=sgp.projections.max.forward.progression.grade,
						max.forward.progression.years=sgp.iter[['sgp.projections.max.forward.progression.years']],
						max.order.for.progression=getMaxOrderForProgression(tail(sgp.iter[["sgp.projection.panel.years"]], 1), tail(sgp.iter[["sgp.projection.content.areas"]], 1),
							state, sgp.projections.equated),
						percentile.trajectory.values=c(1, percentile.trajectory.values, 99),
						panel.data.vnames=getPanelDataVnames("sgp.projections", sgp.iter, sgp.data.names, equate.variable),
						grade.progression=sgp.iter[["sgp.projection.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.projection.content.areas"]],
						year_lags.progression=sgp.iter[["sgp.projection.panel.years.lags"]],
						grade.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						content_area.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						year_lags.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						calculate.sgps=!(tail(sgp.iter[["sgp.projection.panel.years"]], 1) %in%
							SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[tail(sgp.iter[["sgp.projection.content.areas"]], 1)]] &
							is.null(sgp.projections.equated)),
						sgp.exact.grade.progression=sgp.iter[["sgp.exact.grade.progression"]],
						projcuts.digits=SGPstateData[[state]][["SGP_Configuration"]][["projcuts.digits"]],
						projection.unit=sgp.projections.projection.unit,
						projection.unit.label=sgp.projections.projection.unit.label,
						return.projection.group.identifier=sgp.iter[["sgp.projection.sequence"]],
						return.projection.group.scale.scores=return.projection.group.scale.scores,
						return.projection.group.dates=return.projection.group.dates,
						sgp.projections.equated=sgp.projections.equated,
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.projections"),
						...))
				}
				tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
				if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
					tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
						sgp.projections.=getErrorReports(tmp, tmp.tf, par.sgp.config[['sgp.projections']]))
				}
			} else {# END FOREACH
				###   SNOW flavor
				if (par.start$par.type == 'SNOW') {
					tmp <- clusterApplyLB(par.start$internal.cl, par.sgp.config[['sgp.projections']], function(sgp.iter) studentGrowthProjections(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.projections", sgp.iter, sgp.scale.score.equated=equate.variable, SGPt=SGPt),
							Coefficient_Matrices=selectCoefficientMatrices(tmp_sgp_object, coefficient.matrix.type),
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.projections")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.projection.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.projection.content.areas"]], 1)),
						use.my.coefficient.matrices=list(my.year=tail(sgp.iter[["sgp.projection.panel.years"]], 1),
							my.subject=tail(sgp.iter[["sgp.projection.content.areas"]], 1), my.extra.label=equate.label),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.projection.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.projection.content.areas"]], 1)),
						performance.level.cutscores=state,
						max.forward.progression.grade=sgp.projections.max.forward.progression.grade,
						max.forward.progression.years=sgp.iter[['sgp.projections.max.forward.progression.years']],
						max.order.for.progression=getMaxOrderForProgression(tail(sgp.iter[["sgp.projection.panel.years"]], 1), tail(sgp.iter[["sgp.projection.content.areas"]], 1),
							state, sgp.projections.equated),
						percentile.trajectory.values=c(1, percentile.trajectory.values, 99),
						panel.data.vnames=getPanelDataVnames("sgp.projections", sgp.iter, sgp.data.names, equate.variable),
						grade.progression=sgp.iter[["sgp.projection.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.projection.content.areas"]],
						year_lags.progression=sgp.iter[["sgp.projection.panel.years.lags"]],
						grade.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						content_area.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						year_lags.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						calculate.sgps=!(tail(sgp.iter[["sgp.projection.panel.years"]], 1) %in%
							SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[tail(sgp.iter[["sgp.projection.content.areas"]], 1)]] &
							is.null(sgp.projections.equated)),
						sgp.exact.grade.progression=sgp.iter[["sgp.exact.grade.progression"]],
						projcuts.digits=SGPstateData[[state]][["SGP_Configuration"]][["projcuts.digits"]],
						projection.unit=sgp.projections.projection.unit,
						projection.unit.label=sgp.projections.projection.unit.label,
						return.projection.group.identifier=sgp.iter[["sgp.projection.sequence"]],
						return.projection.group.scale.scores=return.projection.group.scale.scores,
						return.projection.group.dates=return.projection.group.dates,
						sgp.projections.equated=sgp.projections.equated,
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.projections"),
						...))

					tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
					if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
						tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
							sgp.projections.=getErrorReports(tmp, tmp.tf, par.sgp.config[['sgp.projections']]))
					}
				} # END SNOW

				###  MULTICORE flavor
				if (par.start$par.type == 'MULTICORE') {
					tmp <- mclapply(par.sgp.config[['sgp.projections']], function(sgp.iter) studentGrowthProjections(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.projections", sgp.iter, sgp.scale.score.equated=equate.variable, SGPt=SGPt),
							Coefficient_Matrices=selectCoefficientMatrices(tmp_sgp_object, coefficient.matrix.type),
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.projections")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.projection.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.projection.content.areas"]], 1)),
						use.my.coefficient.matrices=list(my.year=tail(sgp.iter[["sgp.projection.panel.years"]], 1),
							my.subject=tail(sgp.iter[["sgp.projection.content.areas"]], 1), my.extra.label=equate.label),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.projection.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.projection.content.areas"]], 1)),
						performance.level.cutscores=state,
						max.forward.progression.grade=sgp.projections.max.forward.progression.grade,
						max.forward.progression.years=sgp.iter[['sgp.projections.max.forward.progression.years']],
						max.order.for.progression=getMaxOrderForProgression(tail(sgp.iter[["sgp.projection.panel.years"]], 1), tail(sgp.iter[["sgp.projection.content.areas"]], 1),
							state, sgp.projections.equated),
						percentile.trajectory.values=c(1, percentile.trajectory.values, 99),
						panel.data.vnames=getPanelDataVnames("sgp.projections", sgp.iter, sgp.data.names, equate.variable),
						grade.progression=sgp.iter[["sgp.projection.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.projection.content.areas"]],
						year_lags.progression=sgp.iter[["sgp.projection.panel.years.lags"]],
						grade.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						content_area.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						year_lags.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						calculate.sgps=!(tail(sgp.iter[["sgp.projection.panel.years"]], 1) %in%
							SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[tail(sgp.iter[["sgp.projection.content.areas"]], 1)]] &
							is.null(sgp.projections.equated)),
						sgp.exact.grade.progression=sgp.iter[["sgp.exact.grade.progression"]],
						projcuts.digits=SGPstateData[[state]][["SGP_Configuration"]][["projcuts.digits"]],
						projection.unit=sgp.projections.projection.unit,
						projection.unit.label=sgp.projections.projection.unit.label,
						return.projection.group.identifier=sgp.iter[["sgp.projection.sequence"]],
						return.projection.group.scale.scores=return.projection.group.scale.scores,
						return.projection.group.dates=return.projection.group.dates,
						sgp.projections.equated=sgp.projections.equated,
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.projections"),
						...), mc.cores=par.start$workers, mc.preschedule=FALSE)

					tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
					if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
						tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
							sgp.projections.=getErrorReports(tmp, tmp.tf, par.sgp.config[['sgp.projections']]))
					}
				} # End MULTICORE
			} # END parallel flavors
			stopParallel(parallel.config, par.start)
		} ## END if sgp.projections


	#######################################################
	###  PROJECTIONS (BASELINE referenced)
	#######################################################

		if (sgp.projections.baseline) {
			par.start <- startParallel(parallel.config, 'PROJECTIONS')

			###  FOREACH flavor
			if (toupper(parallel.config[["BACKEND"]]) == "FOREACH") {
				tmp <- foreach(sgp.iter=iter(par.sgp.config[['sgp.projections.baseline']]), .packages="SGP", .errorhandling = "pass", .inorder=FALSE,
					.options.multicore=par.start$foreach.options, .options.mpi=par.start$foreach.options, .options.redis=par.start$foreach.options) %dopar% {
					return(studentGrowthProjections(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.projections.baseline", sgp.iter, SGPt=SGPt),
							Coefficient_Matrices=selectCoefficientMatrices(tmp_sgp_object, "BASELINE"),
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.projections.baseline")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.projection.baseline.panel.years"]], 1),
							my.subject=tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1), my.extra.label="BASELINE"),
						use.my.coefficient.matrices=list(my.year="BASELINE", my.subject=tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1)),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.projection.baseline.panel.years"]], 1),
							my.subject=tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1)),
						performance.level.cutscores=state,
						max.forward.progression.grade=sgp.projections.max.forward.progression.grade,
						max.forward.progression.years=sgp.iter[['sgp.projections.max.forward.progression.years']],
						max.order.for.progression=sgp.projections.baseline.max.order,
						percentile.trajectory.values=c(1, percentile.trajectory.values, 99),
						panel.data.vnames=getPanelDataVnames("sgp.projections.baseline", sgp.iter, sgp.data.names),
						grade.progression=sgp.iter[["sgp.projection.baseline.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.projection.baseline.content.areas"]],
						year_lags.progression=sgp.iter[["sgp.projection.baseline.panel.years.lags"]],
						grade.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						content_area.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						year_lags.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						calculate.sgps=!(tail(sgp.iter[["sgp.projection.baseline.panel.years"]], 1) %in%
							SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1)]] &
							 is.null(sgp.projections.equated)),
						sgp.exact.grade.progression=sgp.iter[["sgp.exact.grade.progression"]],
						projcuts.digits=SGPstateData[[state]][["SGP_Configuration"]][["projcuts.digits"]],
						projection.unit=sgp.projections.projection.unit,
						projection.unit.label=sgp.projections.projection.unit.label,
						return.projection.group.identifier=sgp.iter[["sgp.projection.sequence"]],
						return.projection.group.scale.scores=return.projection.group.scale.scores,
						return.projection.group.dates=return.projection.group.dates,
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.projections.baseline"),
						...))
				}
				tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
				if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
					tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
						sgp.projections.baseline.=getErrorReports(tmp, tmp.tf, par.sgp.config[['sgp.projections.baseline']]))
				}
			} else {# END FOREACH
				###   SNOW flavor
				if (par.start$par.type == 'SNOW') {
					tmp <- clusterApplyLB(par.start$internal.cl, par.sgp.config[['sgp.projections.baseline']], function(sgp.iter) studentGrowthProjections(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.projections.baseline", sgp.iter, SGPt=SGPt),
							Coefficient_Matrices=selectCoefficientMatrices(tmp_sgp_object, "BASELINE"),
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.projections.baseline")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.projection.baseline.panel.years"]], 1),
							my.subject=tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1), my.extra.label="BASELINE"),
						use.my.coefficient.matrices=list(my.year="BASELINE", my.subject=tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1)),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.projection.baseline.panel.years"]], 1),
							my.subject=tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1)),
						performance.level.cutscores=state,
						max.forward.progression.grade=sgp.projections.max.forward.progression.grade,
						max.forward.progression.years=sgp.iter[['sgp.projections.max.forward.progression.years']],
						max.order.for.progression=sgp.projections.baseline.max.order,
						percentile.trajectory.values=c(1, percentile.trajectory.values, 99),
						panel.data.vnames=getPanelDataVnames("sgp.projections.baseline", sgp.iter, sgp.data.names),
						grade.progression=sgp.iter[["sgp.projection.baseline.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.projection.baseline.content.areas"]],
						year_lags.progression=sgp.iter[["sgp.projection.baseline.panel.years.lags"]],
						grade.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						content_area.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						year_lags.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						calculate.sgps=!(tail(sgp.iter[["sgp.projection.baseline.panel.years"]], 1) %in%
							SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1)]] &
							is.null(sgp.projections.equated)),
						sgp.exact.grade.progression=sgp.iter[["sgp.exact.grade.progression"]],
						projcuts.digits=SGPstateData[[state]][["SGP_Configuration"]][["projcuts.digits"]],
						projection.unit=sgp.projections.projection.unit,
						projection.unit.label=sgp.projections.projection.unit.label,
						return.projection.group.identifier=sgp.iter[["sgp.projection.sequence"]],
						return.projection.group.scale.scores=return.projection.group.scale.scores,
						return.projection.group.dates=return.projection.group.dates,
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.projections.baseline"),
						...))

					tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
					if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
						tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
							sgp.projections.baseline.=getErrorReports(tmp, tmp.tf, par.sgp.config[['sgp.projections.baseline']]))
					}
				} # END SNOW

				###  MULTICORE flavor
				if (par.start$par.type == 'MULTICORE') {
					tmp <- mclapply(par.sgp.config[['sgp.projections.baseline']], function(sgp.iter) studentGrowthProjections(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.projections.baseline", sgp.iter, SGPt=SGPt),
							Coefficient_Matrices=selectCoefficientMatrices(tmp_sgp_object, "BASELINE"),
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.projections.baseline")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.projection.baseline.panel.years"]], 1),
							my.subject=tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1), my.extra.label="BASELINE"),
						use.my.coefficient.matrices=list(my.year="BASELINE", my.subject=tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1)),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.projection.baseline.panel.years"]], 1),
							my.subject=tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1)),
						performance.level.cutscores=state,
						max.forward.progression.grade=sgp.projections.max.forward.progression.grade,
						max.forward.progression.years=sgp.iter[['sgp.projections.max.forward.progression.years']],
						max.order.for.progression=sgp.projections.baseline.max.order,
						percentile.trajectory.values=c(1, percentile.trajectory.values, 99),
						panel.data.vnames=getPanelDataVnames("sgp.projections.baseline", sgp.iter, sgp.data.names),
						grade.progression=sgp.iter[["sgp.projection.baseline.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.projection.baseline.content.areas"]],
						year_lags.progression=sgp.iter[["sgp.projection.baseline.panel.years.lags"]],
						grade.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						content_area.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						year_lags.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						calculate.sgps=!(tail(sgp.iter[["sgp.projection.baseline.panel.years"]], 1) %in%
							SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1)]] &
							is.null(sgp.projections.equated)),
						sgp.exact.grade.progression=sgp.iter[["sgp.exact.grade.progression"]],
						projcuts.digits=SGPstateData[[state]][["SGP_Configuration"]][["projcuts.digits"]],
						projection.unit=sgp.projections.projection.unit,
						projection.unit.label=sgp.projections.projection.unit.label,
						return.projection.group.identifier=sgp.iter[["sgp.projection.sequence"]],
						return.projection.group.scale.scores=return.projection.group.scale.scores,
						return.projection.group.dates=return.projection.group.dates,
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.projections.baseline"),
						...), mc.cores=par.start$workers, mc.preschedule=FALSE)

					tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
					if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
						tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
							sgp.projections.baseline.=getErrorReports(tmp, tmp.tf, par.sgp.config[['sgp.projections.baseline']]))
					}
				} # End MULTICORE
			} # END parallel flavors
			stopParallel(parallel.config, par.start)
		} ## END if sgp.projections.baseline


	#################################################
	###  LAGGED PROJECTIONS (COHORT Referenced)
	#################################################

		if (sgp.projections.lagged) {
			par.start <- startParallel(parallel.config, 'LAGGED_PROJECTIONS')

			###  FOREACH flavor
			if (toupper(parallel.config[["BACKEND"]]) == "FOREACH") {
				tmp <- foreach(sgp.iter=iter(par.sgp.config[['sgp.projections.lagged']]), .packages="SGP", .errorhandling = "pass", .inorder=FALSE,
					.options.multicore=par.start$foreach.options, .options.mpi=par.start$foreach.options, .options.redis=par.start$foreach.options) %dopar% {
					return(studentGrowthProjections(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.projections.lagged", sgp.iter, sgp.scale.score.equated=equate.variable, SGPt=SGPt),
							Coefficient_Matrices=selectCoefficientMatrices(tmp_sgp_object, coefficient.matrix.type),
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.projections.lagged")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1),
							my.extra.label="LAGGED"),
						use.my.coefficient.matrices=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1),
							my.subject=tail(sgp.iter[["sgp.content.areas"]], 1), my.extra.label=equate.label),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						performance.level.cutscores=state,
						max.order.for.progression=getMaxOrderForProgression(tail(sgp.iter[["sgp.panel.years"]], 1), tail(sgp.iter[["sgp.content.areas"]], 1), state,
							sgp.projections.equated),
                        percentile.trajectory.values=lagged.percentile.trajectory.values,
						max.forward.progression.grade=sgp.projections.max.forward.progression.grade,
						panel.data.vnames=getPanelDataVnames("sgp.projections.lagged", sgp.iter, sgp.data.names, equate.variable),
						achievement.level.prior.vname=paste("ACHIEVEMENT_LEVEL", tail(head(sgp.iter[["sgp.panel.years"]], -1), 1), tail(head(sgp.iter[["sgp.content.areas"]], -1), 1), sep="."),
						grade.progression=sgp.iter[["sgp.projection.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.projection.content.areas"]],
						year_lags.progression=sgp.iter[["sgp.projection.panel.years.lags"]],
						lag.increment=1,
						grade.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						content_area.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						year_lags.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						calculate.sgps=!(tail(sgp.iter[["sgp.panel.years"]], 1) %in%
							SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[tail(sgp.iter[["sgp.content.areas"]], 1)]] &
							is.null(sgp.projections.equated)),
						sgp.exact.grade.progression=sgp.iter[["sgp.exact.grade.progression"]],
						projcuts.digits=SGPstateData[[state]][["SGP_Configuration"]][["projcuts.digits"]],
						projection.unit=sgp.projections.projection.unit,
						projection.unit.label=sgp.projections.projection.unit.label,
						return.projection.group.identifier=sgp.iter[["sgp.projection.sequence"]],
						return.projection.group.scale.scores=return.projection.group.scale.scores,
						return.projection.group.dates=return.projection.group.dates,
						sgp.projections.equated=sgp.projections.equated,
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.projections.lagged"),
						...))
				}
				tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
				if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
					tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
						sgp.projections.lagged.=getErrorReports(tmp, tmp.tf, par.sgp.config[['sgp.projections.lagged']]))
				}
			} else {# END FOREACH
				###   SNOW flavor
				if (par.start$par.type == 'SNOW') {
					tmp <- clusterApplyLB(par.start$internal.cl, par.sgp.config[['sgp.projections.lagged']], function(sgp.iter) studentGrowthProjections(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.projections.lagged", sgp.iter, sgp.scale.score.equated=equate.variable, SGPt=SGPt),
							Coefficient_Matrices=selectCoefficientMatrices(tmp_sgp_object, coefficient.matrix.type),
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.projections.lagged")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1),
							my.extra.label="LAGGED"),
						use.my.coefficient.matrices=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1),
							my.subject=tail(sgp.iter[["sgp.content.areas"]], 1), my.extra.label=equate.label),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						performance.level.cutscores=state,
						max.order.for.progression=getMaxOrderForProgression(tail(sgp.iter[["sgp.panel.years"]], 1), tail(sgp.iter[["sgp.content.areas"]], 1), state,
							sgp.projections.equated),
                        percentile.trajectory.values=lagged.percentile.trajectory.values,
						max.forward.progression.grade=sgp.projections.max.forward.progression.grade,
						panel.data.vnames=getPanelDataVnames("sgp.projections.lagged", sgp.iter, sgp.data.names, equate.variable),
						achievement.level.prior.vname=paste("ACHIEVEMENT_LEVEL", tail(head(sgp.iter[["sgp.panel.years"]], -1), 1), tail(head(sgp.iter[["sgp.content.areas"]], -1), 1), sep="."),
						grade.progression=sgp.iter[["sgp.projection.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.projection.content.areas"]],
						year_lags.progression=sgp.iter[["sgp.projection.panel.years.lags"]],
						lag.increment=1,
						grade.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						content_area.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						year_lags.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						calculate.sgps=!(tail(sgp.iter[["sgp.panel.years"]], 1) %in%
							SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[tail(sgp.iter[["sgp.content.areas"]], 1)]] &
							is.null(sgp.projections.equated)),
						sgp.exact.grade.progression=sgp.iter[["sgp.exact.grade.progression"]],
						projcuts.digits=SGPstateData[[state]][["SGP_Configuration"]][["projcuts.digits"]],
						projection.unit=sgp.projections.projection.unit,
						projection.unit.label=sgp.projections.projection.unit.label,
						return.projection.group.identifier=sgp.iter[["sgp.projection.sequence"]],
						return.projection.group.scale.scores=return.projection.group.scale.scores,
						return.projection.group.dates=return.projection.group.dates,
						sgp.projections.equated=sgp.projections.equated,
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.projections.lagged"),
						...))

					tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
					if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
						tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
							sgp.projections.lagged.=getErrorReports(tmp, tmp.tf, par.sgp.config[['sgp.projections.lagged']]))
					}
				} # END SNOW

				###  MULTICORE flavor
				if (par.start$par.type == 'MULTICORE') {
					tmp <- mclapply(par.sgp.config[['sgp.projections.lagged']], function(sgp.iter) studentGrowthProjections(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.projections.lagged", sgp.iter, sgp.scale.score.equated=equate.variable, SGPt=SGPt),
							Coefficient_Matrices=selectCoefficientMatrices(tmp_sgp_object, coefficient.matrix.type),
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.projections.lagged")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1),
							my.extra.label="LAGGED"),
						use.my.coefficient.matrices=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1),
							my.subject=tail(sgp.iter[["sgp.content.areas"]], 1), my.extra.label=equate.label),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						performance.level.cutscores=state,
						max.order.for.progression=getMaxOrderForProgression(tail(sgp.iter[["sgp.panel.years"]], 1), tail(sgp.iter[["sgp.content.areas"]], 1), state,
							sgp.projections.equated),
                        percentile.trajectory.values=lagged.percentile.trajectory.values,
						max.forward.progression.grade=sgp.projections.max.forward.progression.grade,
						panel.data.vnames=getPanelDataVnames("sgp.projections.lagged", sgp.iter, sgp.data.names, equate.variable),
						achievement.level.prior.vname=paste("ACHIEVEMENT_LEVEL", tail(head(sgp.iter[["sgp.panel.years"]], -1), 1), tail(head(sgp.iter[["sgp.content.areas"]], -1), 1), sep="."),
						grade.progression=sgp.iter[["sgp.projection.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.projection.content.areas"]],
						year_lags.progression=sgp.iter[["sgp.projection.panel.years.lags"]],
						lag.increment=1,
						grade.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						content_area.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						year_lags.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						calculate.sgps=!(tail(sgp.iter[["sgp.panel.years"]], 1) %in%
							SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[tail(sgp.iter[["sgp.content.areas"]], 1)]] &
							is.null(sgp.projections.equated)),
						sgp.exact.grade.progression=sgp.iter[["sgp.exact.grade.progression"]],
						projcuts.digits=SGPstateData[[state]][["SGP_Configuration"]][["projcuts.digits"]],
						projection.unit=sgp.projections.projection.unit,
						projection.unit.label=sgp.projections.projection.unit.label,
						return.projection.group.identifier=sgp.iter[["sgp.projection.sequence"]],
						return.projection.group.scale.scores=return.projection.group.scale.scores,
						return.projection.group.dates=return.projection.group.dates,
						sgp.projections.equated=sgp.projections.equated,
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.projections.lagged"),
						...), mc.cores=par.start$workers, mc.preschedule=FALSE)

					tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
					if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
						tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
							sgp.projections.lagged.=getErrorReports(tmp, tmp.tf, par.sgp.config[['sgp.projections.lagged']]))
					}
				} # End MULTICORE
			} # END parallel flavors
			stopParallel(parallel.config, par.start)
		} ## END if sgp.projections.lagged


	#################################################
	###  LAGGED PROJECTIONS (BASELINE Referenced)
	#################################################

		if (sgp.projections.lagged.baseline) {
			par.start <- startParallel(parallel.config, 'LAGGED_PROJECTIONS')

			###  FOREACH flavor
			if (toupper(parallel.config[["BACKEND"]]) == "FOREACH") {
				tmp <- foreach(sgp.iter=iter(par.sgp.config[['sgp.projections.lagged.baseline']]), .packages="SGP", .errorhandling = "pass", .inorder=FALSE,
					.options.multicore=par.start$foreach.options, .options.mpi=par.start$foreach.options, .options.redis=par.start$foreach.options) %dopar% {
					return(studentGrowthProjections(
						panel.data=list(
							Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.projections.lagged.baseline", sgp.iter, SGPt=SGPt),
							Coefficient_Matrices=selectCoefficientMatrices(tmp_sgp_object, "BASELINE"),
							Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.projections.lagged")),
						sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1),
							my.extra.label="LAGGED.BASELINE"),
						use.my.coefficient.matrices=list(my.year="BASELINE", my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
						performance.level.cutscores=state,
						max.order.for.progression=sgp.projections.lagged.baseline.max.order,
                        percentile.trajectory.values=lagged.percentile.trajectory.values,
						max.forward.progression.grade=sgp.projections.max.forward.progression.grade,
						panel.data.vnames=getPanelDataVnames("sgp.projections.lagged", sgp.iter, sgp.data.names),
						achievement.level.prior.vname=paste("ACHIEVEMENT_LEVEL", tail(head(sgp.iter[["sgp.panel.years"]], -1), 1), tail(head(sgp.iter[["sgp.content.areas"]], -1), 1), sep="."),
						grade.progression=sgp.iter[["sgp.projection.baseline.grade.sequences"]],
						content_area.progression=sgp.iter[["sgp.projection.baseline.content.areas"]],
						year_lags.progression=sgp.iter[["sgp.projection.baseline.panel.years.lags"]],
						lag.increment=1,
						grade.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						content_area.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						year_lags.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
						calculate.sgps=!(tail(sgp.iter[["sgp.panel.years"]], 1) %in%
							SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[tail(sgp.iter[["sgp.content.areas"]], 1)]] &
							is.null(sgp.projections.equated)),
						sgp.exact.grade.progression=sgp.iter[["sgp.exact.grade.progression"]],
						projcuts.digits=SGPstateData[[state]][["SGP_Configuration"]][["projcuts.digits"]],
						projection.unit=sgp.projections.projection.unit,
						projection.unit.label=sgp.projections.projection.unit.label,
						return.projection.group.identifier=sgp.iter[["sgp.projection.sequence"]],
						return.projection.group.scale.scores=return.projection.group.scale.scores,
						return.projection.group.dates=return.projection.group.dates,
						SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.projections.lagged.baseline"),
						...))
				}
				tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
				if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
					tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
						sgp.projections.lagged.baseline.=getErrorReports(tmp, tmp.tf, par.sgp.config[['sgp.projections.lagged.baseline']]))
				}
			} else {# END FOREACH

			###  SNOW flavor
			if (par.start$par.type == 'SNOW') {
				tmp <- clusterApplyLB(par.start$internal.cl, par.sgp.config[['sgp.projections.lagged.baseline']], function(sgp.iter) studentGrowthProjections(
					panel.data=list(
						Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.projections.lagged.baseline", sgp.iter, SGPt=SGPt),
						Coefficient_Matrices=selectCoefficientMatrices(tmp_sgp_object, "BASELINE"),
						Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.projections.lagged")),
					sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1),
						my.extra.label="LAGGED.BASELINE"),
					use.my.coefficient.matrices=list(my.year="BASELINE", my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
					use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
					performance.level.cutscores=state,
					max.order.for.progression=sgp.projections.lagged.baseline.max.order,
                    percentile.trajectory.values=lagged.percentile.trajectory.values,
					max.forward.progression.grade=sgp.projections.max.forward.progression.grade,
					panel.data.vnames=getPanelDataVnames("sgp.projections.lagged", sgp.iter, sgp.data.names),
					achievement.level.prior.vname=paste("ACHIEVEMENT_LEVEL", tail(head(sgp.iter[["sgp.panel.years"]], -1), 1), tail(head(sgp.iter[["sgp.content.areas"]], -1), 1), sep="."),
					grade.progression=sgp.iter[["sgp.projection.baseline.grade.sequences"]],
					content_area.progression=sgp.iter[["sgp.projection.baseline.content.areas"]],
					year_lags.progression=sgp.iter[["sgp.projection.baseline.panel.years.lags"]],
					lag.increment=1,
					grade.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					content_area.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					year_lags.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					calculate.sgps=!(tail(sgp.iter[["sgp.panel.years"]], 1) %in%
						SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[tail(sgp.iter[["sgp.content.areas"]], 1)]] &
						is.null(sgp.projections.equated)),
					sgp.exact.grade.progression=sgp.iter[["sgp.exact.grade.progression"]],
					projcuts.digits=SGPstateData[[state]][["SGP_Configuration"]][["projcuts.digits"]],
					projection.unit=sgp.projections.projection.unit,
					projection.unit.label=sgp.projections.projection.unit.label,
					return.projection.group.identifier=sgp.iter[["sgp.projection.sequence"]],
					return.projection.group.scale.scores=return.projection.group.scale.scores,
					return.projection.group.dates=return.projection.group.dates,
					SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.projections.lagged.baseline"),
					...))

				tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
				if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
					tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
						sgp.projections.lagged.baseline.=getErrorReports(tmp, tmp.tf, par.sgp.config[['sgp.projections.lagged.baseline']]))
				}
			} # END SNOW

			###  MULTICORE flavor
			if (par.start$par.type == 'MULTICORE') {
				tmp <- mclapply(par.sgp.config[['sgp.projections.lagged.baseline']], function(sgp.iter) studentGrowthProjections(
					panel.data=list(
						Panel_Data=getPanelData(tmp_sgp_data_for_analysis, "sgp.projections.lagged.baseline", sgp.iter, SGPt=SGPt),
						Coefficient_Matrices=selectCoefficientMatrices(tmp_sgp_object, "BASELINE"),
						Knots_Boundaries=getKnotsBoundaries(sgp.iter, state, "sgp.projections.lagged")),
					sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1),
						my.extra.label="LAGGED.BASELINE"),
					use.my.coefficient.matrices=list(my.year="BASELINE", my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
					use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
					performance.level.cutscores=state,
					max.order.for.progression=sgp.projections.lagged.baseline.max.order,
                    percentile.trajectory.values=lagged.percentile.trajectory.values,
					max.forward.progression.grade=sgp.projections.max.forward.progression.grade,
					panel.data.vnames=getPanelDataVnames("sgp.projections.lagged", sgp.iter, sgp.data.names),
					achievement.level.prior.vname=paste("ACHIEVEMENT_LEVEL", tail(head(sgp.iter[["sgp.panel.years"]], -1), 1), tail(head(sgp.iter[["sgp.content.areas"]], -1), 1), sep="."),
					grade.progression=sgp.iter[["sgp.projection.baseline.grade.sequences"]],
					content_area.progression=sgp.iter[["sgp.projection.baseline.content.areas"]],
					year_lags.progression=sgp.iter[["sgp.projection.baseline.panel.years.lags"]],
					lag.increment=1,
					grade.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					content_area.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					year_lags.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					calculate.sgps=!(tail(sgp.iter[["sgp.panel.years"]], 1) %in%
						SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[tail(sgp.iter[["sgp.content.areas"]], 1)]] &
						is.null(sgp.projections.equated)),
					sgp.exact.grade.progression=sgp.iter[["sgp.exact.grade.progression"]],
					projcuts.digits=SGPstateData[[state]][["SGP_Configuration"]][["projcuts.digits"]],
					projection.unit=sgp.projections.projection.unit,
					projection.unit.label=sgp.projections.projection.unit.label,
					return.projection.group.identifier=sgp.iter[["sgp.projection.sequence"]],
					return.projection.group.scale.scores=return.projection.group.scale.scores,
					return.projection.group.dates=return.projection.group.dates,
					SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.projections.lagged.baseline"),
					...), mc.cores=par.start$workers, mc.preschedule=FALSE)

				tmp_sgp_object <- mergeSGP(Reduce(mergeSGP, tmp), tmp_sgp_object)
				if (any(tmp.tf <- sapply(tmp, function(x) identical(class(x), "try-error")))) {
					tmp_sgp_object[['Error_Reports']] <- c(tmp_sgp_object[['Error_Reports']],
						sgp.projections.lagged.baseline.=getErrorReports(tmp, tmp.tf, par.sgp.config[['sgp.projections.lagged.baseline']]))
				}
				} # End MULTICORE
			} # END parallel flavors
			stopParallel(parallel.config, par.start)
		} ## END if sgp.projections.lagged.baseline
	}  ## END if (!is.null(parallel.config))


	################################################################
	################################################################
	###	SEQUENTIAL OPTION (NON-Parallel Option)
	################################################################
	################################################################

	if (is.null(parallel.config)) {

		### sgp.percentiles

		if (sgp.percentiles) {
			for (sgp.iter in rev(par.sgp.config[['sgp.percentiles']])) {

				panel.data <- within(tmp_sgp_object, assign("Panel_Data", getPanelData(tmp_sgp_data_for_analysis, "sgp.percentiles", sgp.iter, csem.variable, SGPt=SGPt)))
				tmp.knots.boundaries <- getKnotsBoundaries(sgp.iter, state, "sgp.percentiles")
				panel.data[["Knots_Boundaries"]][[names(tmp.knots.boundaries)]] <- tmp.knots.boundaries[[names(tmp.knots.boundaries)]]

				tmp_sgp_object <- studentGrowthPercentiles(
					panel.data=panel.data,
					sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
					use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
					use.my.coefficient.matrices=sgp.use.my.coefficient.matrices,
					growth.levels=state,
					panel.data.vnames=getPanelDataVnames("sgp.percentiles", sgp.iter, sgp.data.names),
					additional.vnames.to.return=getPanelDataVnames("sgp.percentiles.to.return", sgp.iter, sgp.data.names),
					grade.progression=sgp.iter[["sgp.grade.sequences"]],
					content_area.progression=sgp.iter[["sgp.content.areas"]],
					year.progression=sgp.iter[["sgp.panel.years"]],
					max.order.for.percentile=SGPstateData[[state]][["SGP_Configuration"]][["max.order.for.percentile"]],
					percentile.cuts=SGPstateData[[state]][["SGP_Configuration"]][["percentile.cuts"]],
					calculate.confidence.intervals=get.simulate.sgps.arg(calculate.confidence.intervals, state, sgp.iter),
					drop.nonsequential.grade.progression.variables=FALSE,
					exact.grade.progression.sequence=sgp.iter[["sgp.exact.grade.progression"]],
					sgp.loss.hoss.adjustment=sgp.loss.hoss.adjustment,
					sgp.cohort.size=tmp.cohort.size,
					return.norm.group.scale.scores=return.norm.group.scale.scores,
					return.norm.group.dates=return.norm.group.dates,
					return.prior.scale.score.standardized=return.prior.scale.score.standardized,
					goodness.of.fit=goodness.of.fit.print.arg,
					goodness.of.fit.minimum.n=SGPstateData[[state]][["SGP_Configuration"]][["goodness.of.fit.minimum.n"]],
					verbose.output=verbose.output,
					print.other.gp=print.other.gp,
					parallel.config=lower.level.parallel.config,
					calculate.simex=get.calculate.simex.arg(sgp.iter[["sgp.calculate.simex"]], sgp.iter),
					max.n.for.coefficient.matrices=max.n.for.coefficient.matrices,
					SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.percentiles"),
                    SGPt.max.time=SGPt.max.time,
                    calculate.sgps=sgp.percentiles.calculate.sgps,
					...)
			}
		} ## END if sgp.percentiles


		### sgp.percentiles.equated

		if (sgp.percentiles.equated) {
			for (sgp.iter in rev(par.sgp.config[['sgp.percentiles.equated']])) {

				panel.data <- within(tmp_sgp_object, assign("Panel_Data", getPanelData(tmp_sgp_data_for_analysis, "sgp.percentiles", sgp.iter, csem.variable, equate.variable, SGPt=SGPt)))
				tmp.knots.boundaries <- getKnotsBoundaries(sgp.iter, state, "sgp.percentiles")
				panel.data[["Knots_Boundaries"]][[names(tmp.knots.boundaries)]] <- tmp.knots.boundaries[[names(tmp.knots.boundaries)]]

				tmp_sgp_object <- studentGrowthPercentiles(
					panel.data=panel.data,
					sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1),
						my.subject=tail(sgp.iter[["sgp.content.areas"]], 1), my.extra.label=equate.label),
					use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
					use.my.coefficient.matrices=sgp.use.my.coefficient.matrices,
					growth.levels=state,
					panel.data.vnames=getPanelDataVnames("sgp.percentiles", sgp.iter, sgp.data.names, equate.variable),
					additional.vnames.to.return=getPanelDataVnames("sgp.percentiles.to.return", sgp.iter, sgp.data.names),
					grade.progression=sgp.iter[["sgp.grade.sequences"]],
					content_area.progression=sgp.iter[["sgp.content.areas"]],
					year.progression=sgp.iter[["sgp.panel.years"]],
					max.order.for.percentile=SGPstateData[[state]][["SGP_Configuration"]][["max.order.for.percentile"]],
					percentile.cuts=SGPstateData[[state]][["SGP_Configuration"]][["percentile.cuts"]],
					calculate.confidence.intervals=get.simulate.sgps.arg(calculate.confidence.intervals, state, sgp.iter),
					drop.nonsequential.grade.progression.variables=FALSE,
					exact.grade.progression.sequence=sgp.iter[["sgp.exact.grade.progression"]],
					sgp.loss.hoss.adjustment=sgp.loss.hoss.adjustment,
					sgp.cohort.size=tmp.cohort.size,
					return.norm.group.scale.scores=return.norm.group.scale.scores,
					return.norm.group.dates=return.norm.group.dates,
					return.prior.scale.score.standardized=return.prior.scale.score.standardized,
					goodness.of.fit=goodness.of.fit.print.arg,
					goodness.of.fit.minimum.n=SGPstateData[[state]][["SGP_Configuration"]][["goodness.of.fit.minimum.n"]],
					verbose.output=verbose.output,
					print.other.gp=print.other.gp,
					parallel.config=lower.level.parallel.config,
					max.n.for.coefficient.matrices=max.n.for.coefficient.matrices,
					sgp.percentiles.equated=sgp.projections.equated,
					SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.percentiles.equated"),
                    SGPt.max.time=SGPt.max.time,
					...)
			}
		} ## END if sgp.percentiles.equated


		## sgp.percentiles.baseline

		if (sgp.percentiles.baseline) {
			for (sgp.iter in rev(par.sgp.config[['sgp.percentiles.baseline']])) {

				panel.data <- within(tmp_sgp_object, assign("Panel_Data", getPanelData(tmp_sgp_data_for_analysis, "sgp.percentiles", sgp.iter, csem.variable, SGPt=SGPt)))
				tmp.knots.boundaries <- getKnotsBoundaries(sgp.iter, state, "sgp.percentiles")
				panel.data[["Knots_Boundaries"]][[names(tmp.knots.boundaries)]] <- tmp.knots.boundaries[[names(tmp.knots.boundaries)]]

				tmp_sgp_object <- studentGrowthPercentiles(
					panel.data=panel.data,
					sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1),
						my.subject=tail(sgp.iter[["sgp.content.areas"]], 1), my.extra.label="BASELINE"),
					use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
					use.my.coefficient.matrices=list(my.year="BASELINE", my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
					growth.levels=state,
					panel.data.vnames=getPanelDataVnames("sgp.percentiles.baseline", sgp.iter, sgp.data.names),
					additional.vnames.to.return=getPanelDataVnames("sgp.percentiles.to.return", sgp.iter, sgp.data.names),
					grade.progression=sgp.iter[["sgp.baseline.grade.sequences"]],
					content_area.progression=sgp.iter[["sgp.baseline.content.areas"]],
					year_lags.progression=sgp.iter[["sgp.baseline.panel.years.lags"]],
					num.prior=min(sgp.iter[["sgp.baseline.max.order"]], sgp.percentiles.baseline.max.order),
					percentile.cuts=SGPstateData[[state]][["SGP_Configuration"]][["percentile.cuts"]],
					calculate.confidence.intervals=get.simulate.sgps.arg(calculate.confidence.intervals, state, sgp.iter),
					drop.nonsequential.grade.progression.variables=FALSE,
					exact.grade.progression.sequence=sgp.iter[["sgp.exact.grade.progression"]],
					sgp.loss.hoss.adjustment=sgp.loss.hoss.adjustment,
					return.norm.group.scale.scores=return.norm.group.scale.scores,
					return.norm.group.dates=return.norm.group.dates,
					return.prior.scale.score.standardized=return.prior.scale.score.standardized,
					goodness.of.fit=goodness.of.fit.print.arg,
					goodness.of.fit.minimum.n=SGPstateData[[state]][["SGP_Configuration"]][["goodness.of.fit.minimum.n"]],
					verbose.output=verbose.output,
					print.other.gp=print.other.gp,
					parallel.config=lower.level.parallel.config,
					calculate.simex=get.calculate.simex.arg(sgp.iter[["sgp.calculate.simex.baseline"]], sgp.iter),
					SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.percentiles.baseline"),
                    SGPt.max.time=SGPt.max.time,
					...)
			}
		} ## END if sgp.percentiles.baseline


		## sgp.projections

		if (sgp.projections) {
			for (sgp.iter in par.sgp.config[['sgp.projections']]) {

				panel.data <- within(tmp_sgp_object, assign("Panel_Data",
								getPanelData(tmp_sgp_data_for_analysis, "sgp.projections", sgp.iter, sgp.scale.score.equated=equate.variable, SGPt=SGPt)))
				tmp.knots.boundaries <- getKnotsBoundaries(sgp.iter, state, "sgp.projections")
				panel.data[["Knots_Boundaries"]][[names(tmp.knots.boundaries)]] <- tmp.knots.boundaries[[names(tmp.knots.boundaries)]]

				tmp_sgp_object <- studentGrowthProjections(
					panel.data=panel.data,
					sgp.labels=list(my.year=tail(sgp.iter[["sgp.projection.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.projection.content.areas"]], 1)),
					use.my.coefficient.matrices=list(my.year=tail(sgp.iter[["sgp.projection.panel.years"]], 1),
						my.subject=tail(sgp.iter[["sgp.projection.content.areas"]], 1), my.extra.label=equate.label),
					use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.projection.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.projection.content.areas"]], 1)),
					performance.level.cutscores=state,
					max.forward.progression.years=sgp.iter[['sgp.projections.max.forward.progression.years']],
					max.forward.progression.grade=sgp.projections.max.forward.progression.grade,
					max.order.for.progression=getMaxOrderForProgression(tail(sgp.iter[["sgp.projection.panel.years"]], 1), tail(sgp.iter[["sgp.projection.content.areas"]], 1), state,
						sgp.projections.equated),
					percentile.trajectory.values=c(1, percentile.trajectory.values, 99),
					panel.data.vnames=getPanelDataVnames("sgp.projections", sgp.iter, sgp.data.names, equate.variable),
					grade.progression=sgp.iter[["sgp.projection.grade.sequences"]],
					content_area.progression=sgp.iter[["sgp.projection.content.areas"]],
					year_lags.progression=sgp.iter[["sgp.projection.panel.years.lags"]],
					grade.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					content_area.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					year_lags.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					calculate.sgps=!(tail(sgp.iter[["sgp.projection.panel.years"]], 1) %in%
						SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[tail(sgp.iter[["sgp.projection.content.areas"]], 1)]] &
						is.null(sgp.projections.equated)),
					sgp.exact.grade.progression=sgp.iter[["sgp.exact.grade.progression"]],
					projcuts.digits=SGPstateData[[state]][["SGP_Configuration"]][["projcuts.digits"]],
					projection.unit=sgp.projections.projection.unit,
					projection.unit.label=sgp.projections.projection.unit.label,
					return.projection.group.identifier=sgp.iter[["sgp.projection.sequence"]],
					return.projection.group.scale.scores=return.projection.group.scale.scores,
					return.projection.group.dates=return.projection.group.dates,
					sgp.projections.equated=sgp.projections.equated,
					SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.projections"),
					...)
			}
		} ## END if sgp.projections


		## sgp.projections.baseline

		if (sgp.projections.baseline) {
			for (sgp.iter in par.sgp.config[['sgp.projections.baseline']]) {

				panel.data <- within(tmp_sgp_object, assign("Panel_Data", getPanelData(tmp_sgp_data_for_analysis, "sgp.projections.baseline", sgp.iter, SGPt=SGPt)))
				tmp.knots.boundaries <- getKnotsBoundaries(sgp.iter, state, "sgp.projections.baseline")
				panel.data[["Knots_Boundaries"]][[names(tmp.knots.boundaries)]] <- tmp.knots.boundaries[[names(tmp.knots.boundaries)]]

				tmp_sgp_object <- studentGrowthProjections(
					panel.data=panel.data,
					sgp.labels=list(my.year=tail(sgp.iter[["sgp.projection.baseline.panel.years"]], 1),
						my.subject=tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1), my.extra.label="BASELINE"),
					use.my.coefficient.matrices=list(my.year="BASELINE", my.subject=tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1)),
					use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.projection.baseline.panel.years"]], 1),
						my.subject=tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1)),
					performance.level.cutscores=state,
					max.forward.progression.years=sgp.iter[['sgp.projections.max.forward.progression.years']],
					max.forward.progression.grade=sgp.projections.max.forward.progression.grade,
					max.order.for.progression=sgp.projections.baseline.max.order,
					percentile.trajectory.values=c(1, percentile.trajectory.values, 99),
					panel.data.vnames=getPanelDataVnames("sgp.projections.baseline", sgp.iter, sgp.data.names),
					grade.progression=sgp.iter[["sgp.projection.baseline.grade.sequences"]],
					content_area.progression=sgp.iter[["sgp.projection.baseline.content.areas"]],
					year_lags.progression=sgp.iter[["sgp.projection.baseline.panel.years.lags"]],
					grade.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					content_area.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					year_lags.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					calculate.sgps=!(tail(sgp.iter[["sgp.projection.baseline.panel.years"]], 1) %in%
						SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[tail(sgp.iter[["sgp.projection.baseline.content.areas"]], 1)]] &
						is.null(sgp.projections.equated)),
					sgp.exact.grade.progression=sgp.iter[["sgp.exact.grade.progression"]],
					projcuts.digits=SGPstateData[[state]][["SGP_Configuration"]][["projcuts.digits"]],
					projection.unit=sgp.projections.projection.unit,
					projection.unit.label=sgp.projections.projection.unit.label,
					return.projection.group.identifier=sgp.iter[["sgp.projection.sequence"]],
					return.projection.group.scale.scores=return.projection.group.scale.scores,
					return.projection.group.dates=return.projection.group.dates,
					SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.projections.baseline"),
					...)
			}
		} ## END if sgp.projections.baseline


		## sgp.projections.lagged

		if (sgp.projections.lagged) {
			for (sgp.iter in par.sgp.config[['sgp.projections.lagged']]) {

				panel.data <- within(tmp_sgp_object, assign("Panel_Data", getPanelData(tmp_sgp_data_for_analysis, "sgp.projections.lagged", sgp.iter, sgp.scale.score.equated=equate.variable, SGPt=SGPt)))
				tmp.knots.boundaries <- getKnotsBoundaries(sgp.iter, state, "sgp.projections.lagged")
				panel.data[["Knots_Boundaries"]][[names(tmp.knots.boundaries)]] <- tmp.knots.boundaries[[names(tmp.knots.boundaries)]]

				tmp_sgp_object <- studentGrowthProjections(
					panel.data=panel.data,
					sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1),
						my.extra.label="LAGGED"),
					use.my.coefficient.matrices=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1),
						my.subject=tail(sgp.iter[["sgp.content.areas"]], 1), my.extra.label=equate.label),
					use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
					performance.level.cutscores=state,
					max.order.for.progression=getMaxOrderForProgression(tail(sgp.iter[["sgp.panel.years"]], 1), tail(sgp.iter[["sgp.content.areas"]], 1), state,
						sgp.projections.equated),
                    percentile.trajectory.values=lagged.percentile.trajectory.values,
					max.forward.progression.grade=sgp.projections.max.forward.progression.grade,
					panel.data.vnames=getPanelDataVnames("sgp.projections.lagged", sgp.iter, sgp.data.names, equate.variable),
					achievement.level.prior.vname=paste("ACHIEVEMENT_LEVEL", tail(head(sgp.iter[["sgp.panel.years"]], -1), 1), tail(head(sgp.iter[["sgp.content.areas"]], -1), 1), sep="."),
					grade.progression=sgp.iter[["sgp.projection.grade.sequences"]],
					content_area.progression=sgp.iter[["sgp.projection.content.areas"]],
					year_lags.progression=sgp.iter[["sgp.projection.panel.years.lags"]],
					lag.increment=1,
					grade.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					content_area.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					year_lags.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					calculate.sgps=!(tail(sgp.iter[["sgp.panel.years"]], 1) %in%
						SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[tail(sgp.iter[["sgp.content.areas"]], 1)]] &
						is.null(sgp.projections.equated)),
					sgp.exact.grade.progression=sgp.iter[["sgp.exact.grade.progression"]],
					projcuts.digits=SGPstateData[[state]][["SGP_Configuration"]][["projcuts.digits"]],
					projection.unit=sgp.projections.projection.unit,
					projection.unit.label=sgp.projections.projection.unit.label,
					return.projection.group.identifier=sgp.iter[["sgp.projection.sequence"]],
					return.projection.group.scale.scores=return.projection.group.scale.scores,
					return.projection.group.dates=return.projection.group.dates,
					sgp.projections.equated=sgp.projections.equated,
					SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.projections.lagged"),
					...)
			}
		} ## END sgp.projections.lagged


		## sgp.projections.lagged.baseline

		if (sgp.projections.lagged.baseline) {
			for (sgp.iter in par.sgp.config[['sgp.projections.lagged.baseline']]) {

				panel.data=within(tmp_sgp_object, assign("Panel_Data", getPanelData(tmp_sgp_data_for_analysis, "sgp.projections.lagged.baseline", sgp.iter, SGPt=SGPt)))
				tmp.knots.boundaries <- getKnotsBoundaries(sgp.iter, state, "sgp.projections.lagged")
				panel.data[["Knots_Boundaries"]][[names(tmp.knots.boundaries)]] <- tmp.knots.boundaries[[names(tmp.knots.boundaries)]]

				tmp_sgp_object <- studentGrowthProjections(
					panel.data=panel.data,
					sgp.labels=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1),
						my.extra.label="LAGGED.BASELINE"),
					use.my.coefficient.matrices=list(my.year="BASELINE", my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
					use.my.knots.boundaries=list(my.year=tail(sgp.iter[["sgp.panel.years"]], 1), my.subject=tail(sgp.iter[["sgp.content.areas"]], 1)),
					performance.level.cutscores=state,
					max.order.for.progression=sgp.projections.lagged.baseline.max.order,
                    percentile.trajectory.values=lagged.percentile.trajectory.values,
					max.forward.progression.grade=sgp.projections.max.forward.progression.grade,
					panel.data.vnames=getPanelDataVnames("sgp.projections.lagged", sgp.iter, sgp.data.names),
					achievement.level.prior.vname=paste("ACHIEVEMENT_LEVEL", tail(head(sgp.iter[["sgp.panel.years"]], -1), 1), tail(head(sgp.iter[["sgp.content.areas"]], -1), 1), sep="."),
					grade.progression=sgp.iter[["sgp.projection.baseline.grade.sequences"]],
					content_area.progression=sgp.iter[["sgp.projection.baseline.content.areas"]],
					year_lags.progression=sgp.iter[["sgp.projection.baseline.panel.years.lags"]],
					lag.increment=1,
					grade.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					content_area.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					year_lags.projection.sequence=SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]],
					calculate.sgps=!(tail(sgp.iter[["sgp.panel.years"]], 1) %in%
						SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]][[tail(sgp.iter[["sgp.content.areas"]], 1)]] &
						is.null(sgp.projections.equated)),
					sgp.exact.grade.progression=sgp.iter[["sgp.exact.grade.progression"]],
					projcuts.digits=SGPstateData[[state]][["SGP_Configuration"]][["projcuts.digits"]],
					projection.unit=sgp.projections.projection.unit,
					projection.unit.label=sgp.projections.projection.unit.label,
					return.projection.group.identifier=sgp.iter[["sgp.projection.sequence"]],
					return.projection.group.scale.scores=return.projection.group.scale.scores,
					return.projection.group.dates=return.projection.group.dates,
					SGPt=getSGPtNames(sgp.iter, SGPt, "sgp.projections.lagged.baseline"),
					...)
			}
		} ## END sgp.projections.lagged.baseline

		tmp_sgp_object[['Panel_Data']] <- NULL

	} ## END sequential analyzeSGP


	if (!keep.sqlite & sgp.sqlite) {if (del.dir) unlink("Data/tmp_data", recursive=TRUE, force=TRUE) else unlink("Data/tmp_data/TMP_SGP_Data.sqlite", recursive=TRUE)}
	sgp_object@SGP <- mergeSGP(tmp_sgp_object, sgp_object@SGP)

	if (goodness.of.fit.print) gof.print(sgp_object)
	setkeyv(sgp_object@Data, getKey(sgp_object)) # re-key data for combineSGP, etc.
	sgp_object@Version[["analyzeSGP"]][[as.character(gsub("-", "_", Sys.Date()))]] <- as.character(packageVersion("SGP"))
	messageSGP(paste("Finished analyzeSGP", date(), "in", convertTime(timetaken(started.at)), "\n"))
	return(sgp_object)
} ## END analyzeSGP Function
