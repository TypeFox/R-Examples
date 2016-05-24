`getSGPConfig` <-
function(sgp_object,
	state,
	tmp_sgp_object,
	content_areas,
	years,
	grades,
	sgp.config,
	trim.sgp.config,
	sgp.percentiles,
	sgp.projections,
	sgp.projections.lagged,
	sgp.percentiles.baseline,
	sgp.projections.baseline,
	sgp.projections.lagged.baseline,
	sgp.config.drop.nonsequential.grade.progression.variables,
	sgp.minimum.default.panel.years,
	sgp.projections.max.forward.progression.years,
	sgp.use.my.coefficient.matrices,
	calculate.simex=NULL,
	calculate.simex.baseline=NULL,
	year.for.equate=NULL,
	sgp.percentiles.equated=FALSE,
	SGPt=NULL) {

	YEAR <- CONTENT_AREA <- VALID_CASE <- NULL

	### Define variables

	sgp.config.list <- list()

	### Check arguments

	if (is.null(sgp.config) & !is.null(grades)) {
		grades <- type.convert(as.character(grades), as.is=TRUE)
		if (!is.numeric(grades)) {
			stop("\tNOTE: Automatic configuration of analyses is currently only available for integer grade levels. Manual specification of 'sgp.config' is required for non-traditional End of Course grade and course progressions.")
		}
	}

	if (is.null(sgp.use.my.coefficient.matrices)) sgp.use.my.coefficient.matrices <- FALSE

	###  If calculate.simex's are FALSE, set to NULL for consistent treatment
	if (is.null(calculate.simex)) calculate.simex <- FALSE
	if (is.null(calculate.simex.baseline)) calculate.simex.baseline <- FALSE

	if (!is.null(year.for.equate)) {
		grades.for.equate <- intersect(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][['Grades_Tested']],
			SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Assessment_Transition"]][[paste('Grades_Tested', year.for.equate, sep=".")]])
		if (!identical(sgp.percentiles.equated, FALSE)) sgp.percentiles.equated <- TRUE
	} else {
		sgp.percentiles.equated <- FALSE
	}

	year.for.scale.change <- SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["Scale_Change"]]


	### get.config function

	get.config <- function(content_area, year, grades) {

		### Data for Years & Grades
		tmp.unique.data <- lapply(sgp_object@Data[SJ("VALID_CASE", content_area), nomatch=0][, c("YEAR", "GRADE"), with=FALSE], function(x) sort(type.convert(unique(x), as.is=TRUE)))

		### Years (sgp.panel.years)
		sgp.panel.years <- as.character(tmp.unique.data$YEAR[1:which(tmp.unique.data$YEAR==year)])

		### Content Areas (sgp.content.areas)
		sgp.content.areas <- rep(content_area, length(sgp.panel.years))

		### Grades (sgp.grade.sequences)
		tmp.last.year.grades <- sort(type.convert(unique(subset(sgp_object@Data, YEAR==tail(sgp.panel.years, 1) & CONTENT_AREA==content_area & VALID_CASE=="VALID_CASE")[['GRADE']]), as.is=TRUE))
		if (!is.numeric(tmp.last.year.grades) | !is.numeric(tmp.unique.data[['GRADE']])) {
			stop("\tNOTE: Automatic 'sgp.config' calculation is only available for integer grade levels. Manual specification of 'sgp.config' is required for non-traditional grade and course progressions.")
		}
		tmp.sgp.grade.sequences <- lapply(tmp.last.year.grades, function(x) tail(tmp.unique.data$GRADE[tmp.unique.data$GRADE <= x], length(tmp.unique.data$YEAR)))
		if (!is.null(grades)) {
			tmp.sgp.grade.sequences <- tmp.sgp.grade.sequences[sapply(tmp.sgp.grade.sequences, function(x) tail(x,1)) %in% grades]
		}
		sgp.grade.sequences <- lapply(tmp.sgp.grade.sequences, function(x) if (length(x) > 1) x[(tail(x,1)-x) <= length(sgp.panel.years)-1])
		sgp.grade.sequences <- sgp.grade.sequences[!unlist(lapply(sgp.grade.sequences, function(x) !length(x) > 1))]
		sgp.grade.sequences <- lapply(sgp.grade.sequences, as.character)

		### Create and return sgp.config
		if ("YEAR_WITHIN" %in% names(sgp_object@Data)) {
			sgp.panel.years.within <- rep("LAST_OBSERVATION", length(sgp.content.areas))
			return(list(
				sgp.content.areas=sgp.content.areas,
				sgp.panel.years=sgp.panel.years,
				sgp.grade.sequences=sgp.grade.sequences,
				sgp.panel.years.within=sgp.panel.years.within))
		} else {
			return(list(
				sgp.content.areas=sgp.content.areas,
				sgp.panel.years=sgp.panel.years,
				sgp.grade.sequences=sgp.grade.sequences
				))
		}
	} ### END get.config


	### get.par.sgp.config function

	get.par.sgp.config <- function(sgp.config) {

		### Utility functions

		split.location <- function(years) sapply(strsplit(years, '_'), length)[1]

		### Set-up

		par.sgp.config <- list()

		### Loop over each element of sgp.config
		for (a in seq_along(sgp.config)) { # now seq_along names so that sgp.config lists can have same names for some elements

			### Convert sgp.grade.sequences to a list if supplied as a vector
			if (is.numeric(sgp.config[[a]][['sgp.grade.sequences']])) sgp.config[[a]][['sgp.grade.sequences']] <- list(sgp.config[[a]][['sgp.grade.sequences']])

			### Loop over grade distinct grade sequences
			b.iter <- seq(from=length(par.sgp.config)+1, length.out=length(sgp.config[[a]][['sgp.grade.sequences']]))
			for (b in seq_along(b.iter)) {

				### Create a per sgp.grade.sequence branch in par.sgp.config list
				par.sgp.config[[b.iter[b]]] <- sgp.config[[a]]
				par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']] <- as.character(sgp.config[[a]][['sgp.grade.sequences']][[b]])

				### Create sgp.exact.grade.progression
				if (!is.null(sgp.config[[a]][['sgp.exact.grade.progression']])) {
					par.sgp.config[[b.iter[b]]][['sgp.exact.grade.progression']] <- sgp.config[[a]][['sgp.exact.grade.progression']][[b]]
				} else {
					par.sgp.config[[b.iter[b]]][['sgp.exact.grade.progression']] <- FALSE
				}

				###  Set sgp.exact.grade.progression=TRUE if using multiple content areas in a single year as priors.
				if (any(duplicated(paste(par.sgp.config[[b.iter[b]]][['sgp.panel.years']], par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']], sep=".")))) {
					par.sgp.config[[b.iter[b]]][['sgp.exact.grade.progression']] <- TRUE
				} else {
					if (is.null(par.sgp.config[[b]][['sgp.exact.grade.progression']])) {
						par.sgp.config[[b.iter[b]]][['sgp.exact.grade.progression']] <- FALSE
					}
				}

				### Create index and re-specify years and content areas from sgp.panel.years and sgp.content.areas
				if (is.integer(type.convert(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']]))) {
					tmp.numeric.grades <- sort(type.convert(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']]))
					grade.span <- seq(min(tmp.numeric.grades), max(tmp.numeric.grades))
					index <- match(tmp.numeric.grades, grade.span)
					if (!sgp.config.drop.nonsequential.grade.progression.variables)  index <- seq_along(index)
					par.sgp.config[[b.iter[b]]][['sgp.panel.years']] <- tail(par.sgp.config[[b.iter[b]]][['sgp.panel.years']], max(index))[index]
					par.sgp.config[[b.iter[b]]][['sgp.content.areas']] <- tail(par.sgp.config[[b.iter[b]]][['sgp.content.areas']], max(index))[index]
					if ('sgp.panel.years.within' %in% names(sgp.config[[a]])) {
						par.sgp.config[[b.iter[b]]][['sgp.panel.years.within']] <- tail(par.sgp.config[[b.iter[b]]][['sgp.panel.years.within']], max(index))[index]
					}
				}

				### Create sgp.panel.years.lags (if NULL)
				if (is.null(sgp.config[[a]][['sgp.panel.years.lags']])) {
					par.sgp.config[[b.iter[b]]][['sgp.panel.years.lags']] <-
						diff(as.numeric(sapply(strsplit(par.sgp.config[[b.iter[b]]][['sgp.panel.years']], '_'), '[', split.location(par.sgp.config[[b.iter[b]]][['sgp.panel.years']]))))
				}

				### Create sgp.projection.grade.sequences (if NULL)
				if (sgp.projections | sgp.projections.lagged) {
					if (is.null(sgp.config[[a]][['sgp.projection.grade.sequences']])) {
						par.sgp.config[[b.iter[b]]][['sgp.projection.grade.sequences']] <- head(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']], -1)
					} else {
						par.sgp.config[[b.iter[b]]][['sgp.projection.grade.sequences']] <- as.character(sgp.config[[a]][['sgp.projection.grade.sequences']][[b]])
					}
				} else par.sgp.config[[b.iter[b]]][['sgp.projection.grade.sequences']] <- NULL # Make NULL here to feed into checkConfig() -> turns it to NA subsequently

				### Create sgp.projection.baseline.grade.sequences (if NULL)
				if (sgp.projections.baseline | sgp.projections.lagged.baseline) {
					if (is.null(sgp.config[[a]][['sgp.projection.baseline.grade.sequences']])) {
						par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.grade.sequences']] <- head(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']], -1)
					} else {
						par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.grade.sequences']] <- as.character(sgp.config[[a]][['sgp.projection.baseline.grade.sequences']][[b]])
					}
				} else par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.grade.sequences']] <- NULL # Make NULL here to feed into checkConfig() -> turns it to NA subsequently

				### Create sgp.projection.content.areas (if NULL)
				if (sgp.projections | sgp.projections.lagged) {
					if (is.null(sgp.config[[a]][['sgp.projection.content.areas']])) {
						par.sgp.config[[b.iter[b]]][['sgp.projection.content.areas']] <- head(par.sgp.config[[b.iter[b]]][['sgp.content.areas']], -1)
					} else {
						if (identical(par.sgp.config[[b.iter[b]]][['sgp.projection.grade.sequences']], "NO_PROJECTIONS")) {
							par.sgp.config[[b.iter[b]]][['sgp.projection.content.areas']] <- as.character(sgp.config[[a]][['sgp.projection.content.areas']])
						} else {
							par.sgp.config[[b.iter[b]]][['sgp.projection.content.areas']] <-
								tail(sgp.config[[a]][['sgp.projection.content.areas']], length(par.sgp.config[[b.iter[b]]][['sgp.projection.grade.sequences']]))
						}
					}
				} else par.sgp.config[[b.iter[b]]][['sgp.projection.content.areas']] <- NA

				### Create sgp.projection.baseline.content.areas (if NULL)
				if (sgp.projections.baseline | sgp.projections.lagged.baseline) {
					if (is.null(sgp.config[[a]][['sgp.projection.baseline.content.areas']])) {
						par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.content.areas']] <- head(par.sgp.config[[b.iter[b]]][['sgp.content.areas']], -1)
					} else {
						if (identical(par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.grade.sequences']], "NO_PROJECTIONS")) {
							par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.content.areas']] <- as.character(sgp.config[[a]][['sgp.projection.content.areas']])
						} else {
							par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.content.areas']] <-
								tail(sgp.config[[a]][['sgp.projection.baseline.content.areas']], length(par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.grade.sequences']]))
						}
					}
				} else par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.content.areas']] <- NA

				### Create sgp.projection.panel.years & sgp.projection.panel.years.lags (if NULL)
				if (sgp.projections | sgp.projections.lagged) {
					if (is.null(sgp.config[[a]][['sgp.projection.panel.years']])) {
						tmp.panel.years <- head(par.sgp.config[[b.iter[b]]][['sgp.panel.years']], -1)
						par.sgp.config[[b.iter[b]]][['sgp.projection.panel.years']] <-
							tail(as.character(sapply(tmp.panel.years, yearIncrement, tail(par.sgp.config[[b.iter[b]]][['sgp.panel.years.lags']], 1))),
								length(par.sgp.config[[b.iter[b]]][['sgp.projection.grade.sequences']]))
						par.sgp.config[[b.iter[b]]][['sgp.projection.panel.years.lags']] <- tail(head(par.sgp.config[[b.iter[b]]][['sgp.panel.years.lags']], -1),
							length(par.sgp.config[[b.iter[b]]][['sgp.projection.grade.sequences']])-1)
					} else {
						if (!identical(par.sgp.config[[b.iter[b]]][['sgp.projection.grade.sequences']], "NO_PROJECTIONS")) {
							par.sgp.config[[b.iter[b]]][['sgp.projection.panel.years.lags']] <-
								diff(as.numeric(sapply(strsplit(par.sgp.config[[b.iter[b]]][['sgp.projection.panel.years']], '_'), '[', split.location(par.sgp.config[[b.iter[b]]][['sgp.projection.panel.years']]))))
						}
					}
				} else par.sgp.config[[b.iter[b]]][['sgp.projection.panel.years']] <- NA

				### Create sgp.projection.baseline.panel.years & sgp.projection.baseline.panel.years.lags (if NULL)
				if (sgp.projections.baseline | sgp.projections.lagged.baseline) {
					if (is.null(sgp.config[[a]][['sgp.projection.baseline.panel.years']])) {
						tmp.panel.years <- head(par.sgp.config[[b.iter[b]]][['sgp.panel.years']], -1)
						par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.panel.years']] <-
							tail(as.character(sapply(tmp.panel.years, yearIncrement, tail(par.sgp.config[[b.iter[b]]][['sgp.panel.years.lags']], 1))),
								length(par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.grade.sequences']]))
						par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.panel.years.lags']] <- tail(head(par.sgp.config[[b.iter[b]]][['sgp.panel.years.lags']], -1),
							length(par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.grade.sequences']])-1)
					} else {
						if (!identical(par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.grade.sequences']], "NO_PROJECTIONS")) {
							par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.panel.years.lags']] <-
								diff(as.numeric(sapply(strsplit(par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.panel.years']], '_'), '[', split.location(par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.panel.years']]))))
						}
					}
				} else par.sgp.config[[b.iter[b]]][['sgp.projection.baseline.panel.years']] <- NA

				### Create sgp.projection.sequence (if NULL)
				if (sgp.projections | sgp.projections.lagged | sgp.projections.baseline | sgp.projections.lagged.baseline) {
					if (is.null(sgp.config[[a]][['sgp.projection.sequence']])) {
						par.sgp.config[[b.iter[b]]][['sgp.projection.sequence']] <- tail(par.sgp.config[[b.iter[b]]][['sgp.content.areas']], 1)
					}
				}

				### Create sgp.projections.max.forward.progression.years (if NULL)
				if (sgp.projections | sgp.projections.lagged | sgp.projections.baseline | sgp.projections.lagged.baseline) {
					if (is.null(sgp.config[[a]][['sgp.projections.max.forward.progression.years']])) {
						par.sgp.config[[b.iter[b]]][['sgp.projections.max.forward.progression.years']] <- sgp.projections.max.forward.progression.years
					}
				}

				### Create sgp.calculate.simex (if requested) - else leave NULL
				tmp.calculate.simex <- calculate.simex
				tmp.calculate.simex.baseline <- calculate.simex.baseline
				if (is.list(tmp.calculate.simex)) {
					par.sgp.config[[b.iter[b]]][['sgp.calculate.simex']] <- tmp.calculate.simex
					tmp.calculate.simex <- TRUE
				}

				if (!is.null(sgp.config[[a]][['sgp.calculate.simex']])) {
					if (is.logical(sgp.config[[a]][['sgp.calculate.simex']])) {
						if (sgp.config[[a]][['sgp.calculate.simex']]) {
							par.sgp.config[[b.iter[b]]][['sgp.calculate.simex']] <- list(
								state=state, lambda=seq(0,2,0.5), simulation.iterations=75, simex.sample.size=5000, extrapolation="linear", save.matrices=TRUE)
						} else par.sgp.config[[b.iter[b]]][['sgp.calculate.simex']] <- NULL # FALSE - same as NULL
					} else {
					# NOTE: if SIMEX configs are supplied to analyzeSGP in both a calculate.simex argument (as list) AND in a sgp.config element, the sgp.config is chosen here:
						par.sgp.config[[b.iter[b]]][['sgp.calculate.simex']] <- sgp.config[[a]][['sgp.calculate.simex']]
					}

					if (is.logical(par.sgp.config[[b.iter[b]]][['sgp.calculate.simex']][['simex.use.my.coefficient.matrices']])) {
						if (!par.sgp.config[[b.iter[b]]][['sgp.calculate.simex']][['simex.use.my.coefficient.matrices']]) {
							par.sgp.config[[b.iter[b]]][['sgp.calculate.simex']][['simex.use.my.coefficient.matrices']] <- NULL
						} else {
							sgp.use.my.coefficient.matrices <- TRUE # Looks for existing "naive" matrices if using pre-calc'd SIMEX matrices
						}
					}
					if (!is.null(par.sgp.config[[b.iter[b]]][['sgp.calculate.simex']][['simex.use.my.coefficient.matrices']])) {
						tmp.matrix.label <- paste(tail(par.sgp.config[[b.iter[b]]][['sgp.content.areas']], 1),
							tail(par.sgp.config[[b.iter[b]]][['sgp.panel.years']], 1), "SIMEX", sep=".")
						tmp.orders <- getsplineMatrices(
							my.matrices=tmp_sgp_object[['Coefficient_Matrices']][[tmp.matrix.label]],
							my.matrix.content_area.progression=par.sgp.config[[b.iter[b]]][['sgp.content.areas']],
							my.matrix.grade.progression=par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']],
							my.matrix.time.progression=par.sgp.config[[b.iter[b]]][['sgp.panel.years']],
							my.matrix.time.progression.lags=par.sgp.config[[b.iter[b]]][['sgp.panel.years.lags']],
							my.matrix.time.dependency=if (is.null(SGPt)) NULL else list(TIME="TIME", TIME_LAG="TIME_LAG"),
							what.to.return="ORDERS")
						tmp.max.order <- max(tmp.orders)

						for (L in par.sgp.config[[b.iter[b]]][['sgp.calculate.simex']][['lambda']][-1]) {
							if (par.sgp.config[[b.iter[b]]][['sgp.exact.grade.progression']]) ord.iter <- tmp.max.order else ord.iter <- seq_along(tmp.orders)
							for (k in ord.iter) {
							par.sgp.config[[b.iter[b]]][['sgp.matrices']][[paste(tmp.matrix.label, ".SIMEX", sep="")]][[
								paste("qrmatrices", tail(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']], 1), k, sep="_")]][[paste("lambda_", L, sep="")]] <-
							unlist(getsplineMatrices(
								my.matrices=tmp_sgp_object[['Coefficient_Matrices']][[paste(tmp.matrix.label, ".SIMEX", sep="")]][[
									paste("qrmatrices", tail(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']], 1), k, sep="_")]][[paste("lambda_", L, sep="")]],
								my.matrix.content_area.progression=par.sgp.config[[b.iter[b]]][['sgp.content.areas']],
								my.matrix.grade.progression=par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']],
								my.matrix.time.progression=par.sgp.config[[b.iter[b]]][['sgp.panel.years']],
								my.matrix.time.progression.lags=par.sgp.config[[b.iter[b]]][['sgp.panel.years.lags']],
								my.matrix.time.dependency=if (is.null(SGPt)) NULL else list(TIME="TIME", TIME_LAG="TIME_LAG"),
								my.exact.grade.progression.sequence=TRUE,
								return.multiple.matrices=TRUE,
								my.matrix.order=k), recursive=FALSE)
							}
						}
					}
				}

				if (sgp.use.my.coefficient.matrices) {
					tmp.matrix.label <- paste(tail(par.sgp.config[[b.iter[b]]][['sgp.content.areas']], 1), tail(par.sgp.config[[b.iter[b]]][['sgp.panel.years']], 1), sep=".")
					tmp.orders <- getsplineMatrices(
						my.matrices=tmp_sgp_object[['Coefficient_Matrices']][[tmp.matrix.label]],
						my.matrix.content_area.progression=par.sgp.config[[b.iter[b]]][['sgp.content.areas']],
						my.matrix.grade.progression=par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']],
						my.matrix.time.progression=par.sgp.config[[b.iter[b]]][['sgp.panel.years']],
						my.matrix.time.progression.lags=par.sgp.config[[b.iter[b]]][['sgp.panel.years.lags']],
						my.exact.grade.progression.sequence=par.sgp.config[[b.iter[b]]][['sgp.exact.grade.progression']],
						my.matrix.time.dependency=if (is.null(SGPt)) NULL else list(TIME="TIME", TIME_LAG="TIME_LAG"),
						what.to.return="ORDERS")

						if (length(tmp.orders)==0) {
							message("\tNOTE: Cohort Referenced Coefficient matrices are not available for ",
							paste(tmp.matrix.label, paste(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']], collapse=", "), sep=": "), ".", sep="")
						} else {
							tmp.max.order <- max(tmp.orders)
							if (par.sgp.config[[b.iter[b]]][['sgp.exact.grade.progression']]) ord.iter <- tmp.max.order else ord.iter <- seq_along(tmp.orders)
							for (k in ord.iter) {
								par.sgp.config[[b.iter[b]]][['sgp.matrices']][[tmp.matrix.label]][[
								paste("qrmatrices", tail(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']], 1), k, sep="_")]] <- getsplineMatrices(
									my.matrices=tmp_sgp_object[['Coefficient_Matrices']][[tmp.matrix.label]],
									my.matrix.content_area.progression=tail(par.sgp.config[[b.iter[b]]][['sgp.content.areas']], k+1),
									my.matrix.grade.progression=tail(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']], k+1),
									my.matrix.time.progression=tail(par.sgp.config[[b.iter[b]]][['sgp.panel.years']], k+1),
									my.matrix.time.progression.lags=tail(par.sgp.config[[b.iter[b]]][['sgp.panel.years.lags']], k+1),
									my.matrix.time.dependency=if (is.null(SGPt)) NULL else list(TIME="TIME", TIME_LAG="TIME_LAG"),
									my.exact.grade.progression.sequence=TRUE, my.matrix.order=k)[[1]]
						}
					}

					if (sgp.percentiles.equated) {
						tmp.matrix.label <- paste(tail(par.sgp.config[[b.iter[b]]][['sgp.content.areas']], 1), tail(par.sgp.config[[b.iter[b]]][['sgp.panel.years']], 1), "EQUATED", sep=".")
						tmp.orders <- getsplineMatrices(
							my.matrices=tmp_sgp_object[['Coefficient_Matrices']][[tmp.matrix.label]],
							my.matrix.content_area.progression=par.sgp.config[[b.iter[b]]][['sgp.content.areas']],
							my.matrix.grade.progression=par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']],
							my.matrix.time.progression=par.sgp.config[[b.iter[b]]][['sgp.panel.years']],
							my.matrix.time.progression.lags=par.sgp.config[[b.iter[b]]][['sgp.panel.years.lags']],
							my.exact.grade.progression.sequence=par.sgp.config[[b.iter[b]]][['sgp.exact.grade.progression']],
							my.matrix.time.dependency=if (is.null(SGPt)) NULL else list(TIME="TIME", TIME_LAG="TIME_LAG"),
							what.to.return="ORDERS")

						if (length(tmp.orders)==0) {
							message("\tNOTE: Equated Coefficient matrices are not available for ",
							paste(tmp.matrix.label, paste(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']], collapse=", "), sep=": "), ".", sep="")
						} else {
							tmp.max.order <- max(tmp.orders)
							if (par.sgp.config[[b.iter[b]]][['sgp.exact.grade.progression']]) ord.iter <- tmp.max.order else ord.iter <- seq_along(tmp.orders)
							for (k in ord.iter) {
								par.sgp.config[[b.iter[b]]][['sgp.equated.matrices']][[tmp.matrix.label]][[
								paste("qrmatrices", tail(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']], 1), k, sep="_")]] <- getsplineMatrices(
									my.matrices=tmp_sgp_object[['Coefficient_Matrices']][[tmp.matrix.label]],
									my.matrix.content_area.progression=tail(par.sgp.config[[b.iter[b]]][['sgp.content.areas']], k+1),
									my.matrix.grade.progression=tail(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']], k+1),
									my.matrix.time.progression=tail(par.sgp.config[[b.iter[b]]][['sgp.panel.years']], k+1),
									my.matrix.time.progression.lags=tail(par.sgp.config[[b.iter[b]]][['sgp.panel.years.lags']], k+1),
									my.matrix.time.dependency=if (is.null(SGPt)) NULL else list(TIME="TIME", TIME_LAG="TIME_LAG"),
									my.exact.grade.progression.sequence=TRUE, my.matrix.order=k)[[1]]
							}
						}
					}
				}

				### Create sgp.calculate.simex.baseline (if requested) - else leave NULL
				if (is.list(tmp.calculate.simex.baseline)) {
					par.sgp.config[[b.iter[b]]][['sgp.calculate.simex.baseline']] <- tmp.calculate.simex.baseline
					tmp.calculate.simex.baseline <- TRUE
				}

				if (!is.null(sgp.config[[a]][['sgp.calculate.simex.baseline']])) {
					if (is.logical(sgp.config[[a]][['sgp.calculate.simex.baseline']])) {
						if (sgp.config[[a]][['sgp.calculate.simex.baseline']]) {
							par.sgp.config[[b.iter[b]]][['sgp.calculate.simex.baseline']] <- list(state=state, lambda=seq(0,2,0.5), simulation.iterations=75,
								simex.sample.size=5000, extrapolation="linear", save.matrices=FALSE, simex.use.my.coefficient.matrices = TRUE)
						} else par.sgp.config[[b.iter[b]]][['sgp.calculate.simex.baseline']] <- NULL # FALSE - same as NULL
					} else {
						par.sgp.config[[b.iter[b]]][['sgp.calculate.simex.baseline']] <- sgp.config[[a]][['sgp.calculate.simex.baseline']]
					}
				}

				### Create baseline specific arguments
				if (sgp.percentiles.baseline | sgp.projections.baseline | sgp.projections.lagged.baseline) {
					tmp.matrix.label <- paste(strsplit(names(sgp.config)[a], "\\.")[[1]][1], ".BASELINE", sep="")
 					if (tmp.matrix.label %in% names(tmp_sgp_object[["Coefficient_Matrices"]])) {
						tmp.orders <- getsplineMatrices(
							my.matrices=tmp_sgp_object[['Coefficient_Matrices']][[tmp.matrix.label]],
							my.matrix.content_area.progression=par.sgp.config[[b.iter[b]]][['sgp.content.areas']],
							my.matrix.grade.progression=par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']],
							my.matrix.time.progression=rep("BASELINE", length(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']])),
							my.matrix.time.progression.lags=par.sgp.config[[b.iter[b]]][['sgp.panel.years.lags']],
							my.matrix.time.dependency=if (is.null(SGPt)) NULL else list(TIME="TIME", TIME_LAG="TIME_LAG"),
							what.to.return="ORDERS")

						if (length(tmp.orders) > 0) {
							tmp.matrices.tf <- TRUE
							if (!is.null(year.for.scale.change[[tail(par.sgp.config[[b.iter[b]]][['sgp.content.areas']], 1)]]) &&
								year.for.scale.change[[tail(par.sgp.config[[b.iter[b]]][['sgp.content.areas']], 1)]] < tail(par.sgp.config[[b.iter[b]]][['sgp.panel.years']], 1)) {
								tmp.max.order <-
									min(max(tmp.orders),
										as.numeric(tail(par.sgp.config[[b.iter[b]]][['sgp.panel.years']], 1))-
										as.numeric(tail(unlist(strsplit(year.for.scale.change[[tail(par.sgp.config[[b.iter[b]]][['sgp.content.areas']], 1)]], "_")), 1)))
							} else {
								tmp.max.order <- max(tmp.orders)
							}

							if (par.sgp.config[[b.iter[b]]][['sgp.exact.grade.progression']]) ord.iter <- tmp.max.order else ord.iter <- seq_along(tmp.orders)
							for (k in ord.iter) {
								par.sgp.config[[b.iter[b]]][['sgp.baseline.matrices']][[tmp.matrix.label]][[
									paste("qrmatrices", tail(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']], 1), k, sep="_")]] <- getsplineMatrices(
									my.matrices=tmp_sgp_object[['Coefficient_Matrices']][[tmp.matrix.label]],
									my.matrix.content_area.progression=tail(par.sgp.config[[b.iter[b]]][['sgp.content.areas']], k+1),
									my.matrix.grade.progression=tail(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']], k+1),
									my.matrix.time.progression=tail(rep("BASELINE", length(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']])), k+1),
									my.matrix.time.progression.lags=tail(par.sgp.config[[b.iter[b]]][['sgp.panel.years.lags']], k+1),
									my.matrix.time.dependency=if (is.null(SGPt)) NULL else list(TIME="TIME", TIME_LAG="TIME_LAG"),
									my.exact.grade.progression.sequence=TRUE, my.matrix.order=k)[[1]]
							}

							if (!is.null(par.sgp.config[[b.iter[b]]][['sgp.calculate.simex.baseline']])) {
								for (L in par.sgp.config[[b.iter[b]]][['sgp.calculate.simex.baseline']][['lambda']][-1]) {
									for (k in ord.iter) {
									par.sgp.config[[b.iter[b]]][['sgp.baseline.matrices']][[paste(tmp.matrix.label, ".SIMEX", sep="")]][[
										paste("qrmatrices", tail(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']], 1), k, sep="_")]][[paste("lambda_", L, sep="")]] <-
									unlist(getsplineMatrices(
										my.matrices=tmp_sgp_object[['Coefficient_Matrices']][[paste(tmp.matrix.label, ".SIMEX", sep="")]][[
											paste("qrmatrices", tail(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']], 1), k, sep="_")]][[paste("lambda_", L, sep="")]],
										my.matrix.content_area.progression=par.sgp.config[[b.iter[b]]][['sgp.content.areas']],
										my.matrix.grade.progression=par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']],
										my.matrix.time.progression=rep("BASELINE", length(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']])),
										my.matrix.time.progression.lags=par.sgp.config[[b.iter[b]]][['sgp.panel.years.lags']],
										my.matrix.time.dependency=if (is.null(SGPt)) NULL else list(TIME="TIME", TIME_LAG="TIME_LAG"),
										my.exact.grade.progression.sequence=TRUE,
										return.multiple.matrices=TRUE,
										my.matrix.order=k), recursive=FALSE)
									}
								}
							}
						} else {
							tmp.matrices.tf <- FALSE
						}
					} else tmp.matrices.tf <- FALSE

					if (!tmp.matrices.tf) {
						par.sgp.config[[b.iter[b]]][['sgp.baseline.grade.sequences']] <- "NO_BASELINE_COEFFICIENT_MATRICES"
						par.sgp.config[[b.iter[b]]][['sgp.baseline.max.order']] <- "NO_BASELINE_COEFFICIENT_MATRICES"
					} else {
						par.sgp.config[[b.iter[b]]][['sgp.baseline.grade.sequences']] <- as.character(tail(par.sgp.config[[b.iter[b]]][['sgp.grade.sequences']], tmp.max.order+1))
						par.sgp.config[[b.iter[b]]][['sgp.baseline.content.areas']] <- as.character(tail(par.sgp.config[[b.iter[b]]][['sgp.content.areas']], tmp.max.order+1))
						par.sgp.config[[b.iter[b]]][['sgp.baseline.max.order']] <- tmp.max.order
						par.sgp.config[[b.iter[b]]][['sgp.baseline.panel.years.lags']] <- tail(par.sgp.config[[b.iter[b]]][['sgp.panel.years.lags']], tmp.max.order)
					}
				} ### END if (sgp.percentiles.baseline | sgp.projections.baseline | sgp.projections.lagged.baseline

			} ### END b loop
		} ### END a loop
		return(par.sgp.config)
	} ## END get.par.sgp.config


	## test.projection.iter function

	test.projection.iter <- function(sgp.iter) {
		if (identical(sgp.iter[['sgp.projection.grade.sequences']], "NO_PROJECTIONS")) return(FALSE)
		if (!is.null(SGP::SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]]) & !is.null(sgp.iter[["sgp.projection.sequence"]])) {
			if (tail(sgp.iter[["sgp.grade.sequences"]], 1) == "EOCT") { # Only check EOCT configs/iters
				if (is.null(SGP::SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]])) return(FALSE)
				tmp.index <- match(tail(sgp.iter[["sgp.content.areas"]], 1),
					SGP::SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]])
				tmp.content_area.projection.sequence <-
					SGP::SGPstateData[[state]][["SGP_Configuration"]][["content_area.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]][1:tmp.index]
				tmp.grade.projection.sequence <-
					SGP::SGPstateData[[state]][["SGP_Configuration"]][["grade.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]][1:tmp.index]
				tmp.year_lags.projection.sequence <-
					SGP::SGPstateData[[state]][["SGP_Configuration"]][["year_lags.projection.sequence"]][[sgp.iter[["sgp.projection.sequence"]]]][1:(tmp.index-1)]
				if (all(identical(sgp.iter[["sgp.content.areas"]], tail(tmp.content_area.projection.sequence, length(sgp.iter[["sgp.content.areas"]]))) &
					identical(sgp.iter[["sgp.grade.sequences"]], tail(tmp.grade.projection.sequence, length(sgp.iter[["sgp.grade.sequences"]]))) &
					identical(as.numeric(sgp.iter[["sgp.panel.years.lags"]]), as.numeric(tail(tmp.year_lags.projection.sequence, length(sgp.iter[["sgp.panel.years.lags"]])))))) {
					 iter.test <- TRUE
				} else {
					iter.test <- FALSE
				}
			}	else iter.test <- TRUE
		}	else iter.test <- TRUE
		return(iter.test)
	}


	###
	### Construct sgp.config/par.sgp.config
	###

	if (is.null(sgp.config)) {
		tmp.sgp.config <- tmp.years <- list()
		if (is.null(content_areas)) {
			content_areas <- unique(sgp_object@Data["VALID_CASE"][['CONTENT_AREA']])
		}
		if (is.null(years)) {
			for (i in content_areas) {
				tmp.years[[i]] <- sort(tail(unique(sgp_object@Data[SJ("VALID_CASE", i)][['YEAR']]), - (sgp.minimum.default.panel.years-1)), decreasing=TRUE)
			}
		} else {
			for (i in content_areas) {
				tmp.years[[i]] <- years
			}
		}
		for (i in content_areas) {
			for (j in tmp.years[[i]]) {
				tmp.sgp.config[[paste(i,j,sep=".")]] <- get.config(i,j,grades)
			}
		}
		par.sgp.config <- checkConfig(get.par.sgp.config(tmp.sgp.config), "Standard")
	} else {
		par.sgp.config <- checkConfig(get.par.sgp.config(sgp.config), "Standard")
	}


	###
	### Extend sgp.config.list
	###

	if (sgp.percentiles | sgp.percentiles.equated) {
		tmp.config <- par.sgp.config
		for (i in 1:length(tmp.config)) tmp.config[[i]][['sgp.baseline.matrices']] <- NULL
		tmp.config <- tmp.config[which(sapply(tmp.config, function(x) !identical(x[['sgp.grade.sequences']], "NO_PERCENTILES")))]
		if (sgp.use.my.coefficient.matrices) {
			tmp.config <- tmp.config[which(sapply(tmp.config, function(x) !is.null(x[['sgp.matrices']])))]
		}
		if (sgp.percentiles) sgp.config.list[['sgp.percentiles']] <- tmp.config
		if (sgp.percentiles.equated) {
			tmp.config <- tmp.config[sapply(tmp.config, function(x) tail(x[['sgp.panel.years']], 1))==year.for.equate]
			if (!is.null(grades.for.equate)) {
				tmp.config <- tmp.config[which(sapply(tmp.config, function(x) tail(x[['sgp.grade.sequences']], 1) %in% grades.for.equate))]
			}
			sgp.config.list[['sgp.percentiles.equated']] <- tmp.config
		}
	}

	if (sgp.projections | sgp.projections.lagged) {
		tmp.config <- par.sgp.config[sapply(par.sgp.config, test.projection.iter)]
		if (length(tmp.config) > 0) {
			for (f in 1:length(tmp.config)) tmp.config[[f]][['sgp.exact.grade.progression']] <- FALSE
			while (any(sapply(tmp.config, function(x) length(x[['sgp.projection.sequence']])>1))) {
				tmp.index <- which(any(sapply(tmp.config, function(x) length(x[['sgp.projection.sequence']])>1)))[1]
				tmp.iter <- tmp.config[[tmp.index]]
				tmp.config <- tmp.config[-tmp.index]
				tmp.expand.config <- list()
				for (j in 1:length(tmp.iter[['sgp.projection.sequence']])) {
					tmp.expand.config[[j]] <- tmp.iter
					tmp.expand.config[[j]][['sgp.projection.sequence']] <- tmp.expand.config[[j]][['sgp.projection.sequence']][j]
				}
				tmp.config <- c(tmp.config, tmp.expand.config)
			}
			if (!is.null(year.for.equate)) tmp.config <- tmp.config[which(sapply(tmp.config, function(x) tail(x[['sgp.grade.sequences']], 1) %in% grades.for.equate))]
			if (sgp.projections) {
				sgp.config.list[['sgp.projections']] <- tmp.config
				for (i in 1:length(sgp.config.list[['sgp.projections']])) {
					sgp.config.list[['sgp.projections']][[i]][['sgp.matrices']] <- sgp.config.list[['sgp.projections']][[i]][['sgp.baseline.matrices']] <- NULL
				}
			}
			if (sgp.projections.lagged) {
				sgp.config.list[['sgp.projections.lagged']] <- tmp.config
				for (i in 1:length(sgp.config.list[['sgp.projections.lagged']])) {
					sgp.config.list[['sgp.projections.lagged']][[i]][['sgp.matrices']] <- sgp.config.list[['sgp.projections.lagged']][[i]][['sgp.baseline.matrices']] <- NULL
				}
			}
		} else {
			sgp.config.list[['sgp.projections']] <- sgp.config.list[['sgp.projections.lagged']] <- NULL
			message("\n NOTE:  No valid projections have been identified in the sgp.config lists provided -- ",
				"'sgp.projections'"[sgp.projections], " and "[sgp.projections & sgp.projections.lagged], "'sgp.projections.lagged'"[sgp.projections.lagged], " will NOT be calculated.",
				"  \n\tPlease check SGPstateData[['", state, "']][['SGP_Configuration']] for proper 'content_area.projection.sequence', 'grade.projection.sequence' and 'year_lags.projection.sequence' elements.\n")
		}
	}

	if (sgp.percentiles.baseline | sgp.projections.baseline | sgp.projections.lagged.baseline) {

		sgp.config.list[['sgp.percentiles.baseline']] <- par.sgp.config

		if (any(sapply(par.sgp.config, function(x) identical(x[['sgp.baseline.grade.sequences']], "NO_BASELINE_COEFFICIENT_MATRICES")))) {
			baseline.missings <- setdiff(which(sapply(par.sgp.config, function(x) identical(x[['sgp.baseline.grade.sequences']], "NO_BASELINE_COEFFICIENT_MATRICES"))),
						which(sapply(par.sgp.config, function(x) identical(x[['sgp.grade.sequences']], "NO_PERCENTILES"))))
			if (length(baseline.missings)>0) {
				baseline.missings <- paste(unlist(sapply(baseline.missings, function(x)
					paste(tail(par.sgp.config[[x]]$sgp.content.areas, 1), paste(par.sgp.config[[x]]$sgp.grade.sequences, collapse=", "), sep=": "))), collapse=";\n\t\t")
				message("\tNOTE: Baseline coefficient matrices are not available for:\n\t\t", baseline.missings, ".", sep="")

				sgp.config.list[['sgp.percentiles.baseline']] <-
					par.sgp.config[which(sapply(par.sgp.config, function(x) !identical(x[['sgp.baseline.grade.sequences']], "NO_BASELINE_COEFFICIENT_MATRICES")))]
			}
		}

		if (length(sgp.config.list[['sgp.percentiles.baseline']]) > 0) {
			for (i in seq_along(sgp.config.list[['sgp.percentiles.baseline']])) sgp.config.list[['sgp.percentiles.baseline']][[i]][['sgp.matrices']] <- NULL
			tmp.config <- sgp.config.list[['sgp.percentiles.baseline']][sapply(sgp.config.list[['sgp.percentiles.baseline']], test.projection.iter)]
			sgp.config.list[['sgp.percentiles.baseline']] <- tmp.config
		}

		if (sgp.projections.baseline | sgp.projections.lagged.baseline) {
			while (any(sapply(tmp.config, function(x) length(x$sgp.projection.sequence)>1))) {
				tmp.index <- which(any(sapply(tmp.config, function(x) length(x$sgp.projection.sequence)>1)))[1]
				tmp.iter <- tmp.config[[tmp.index]]
				tmp.config <- tmp.config[-tmp.index]
				tmp.expand.config <- list()
				for (j in seq_along(tmp.iter$sgp.projection.sequence)) {
					tmp.expand.config[[j]] <- tmp.iter
					tmp.expand.config[[j]]$sgp.projection.sequence <- tmp.expand.config[[j]]$sgp.projection.sequence[j]
				}
				tmp.config <- c(tmp.config, tmp.expand.config)
			}
		}

		if (!sgp.percentiles.baseline | length(sgp.config.list[['sgp.percentiles.baseline']]) == 0) sgp.config.list[['sgp.percentiles.baseline']] <- NULL
		if (sgp.projections.baseline) {
			sgp.config.list[['sgp.projections.baseline']] <- tmp.config
			for (i in seq_along(sgp.config.list[['sgp.projections.baseline']])) sgp.config.list[['sgp.projections.baseline']][[i]][['sgp.baseline.matrices']] <- NULL
		}
		if (sgp.projections.lagged.baseline) {
			sgp.config.list[['sgp.projections.lagged.baseline']] <- tmp.config
			for (i in seq_along(sgp.config.list[['sgp.projections.lagged.baseline']])) sgp.config.list[['sgp.projections.lagged.baseline']][[i]][['sgp.baseline.matrices']] <- NULL
		}
	}

	### Trim sgp.config if requested

	if (trim.sgp.config) {
		tmp.iter <- c('sgp.percentiles', 'sgp.percentiles.baseline', 'sgp.projections', 'sgp.projections.baseline', 'sgp.projections.lagged', 'sgp.projections.lagged.baseline')
		tmp.iter.tf <- c(sgp.percentiles, sgp.percentiles.baseline, sgp.projections, sgp.projections.baseline, sgp.projections.lagged, sgp.projections.lagged.baseline)
		for (i in tmp.iter[tmp.iter.tf]) {
			if (i=='sgp.projections') {
				tmp.content.areas.label <- 'sgp.projection.content.areas'
				tmp.panel.years.label <- 'sgp.projection.panel.years'
				tmp.grade.sequences.label <- 'sgp.projection.grade.sequences'
			}
			if (i=='sgp.projections.baseline') {
				tmp.content.areas.label <- 'sgp.projection.baseline.content.areas'
				tmp.panel.years.label <- 'sgp.projection.baseline.panel.years'
				tmp.grade.sequences.label <- 'sgp.projection.baseline.grade.sequences'
			}
			if (i %in% c('sgp.percentiles', 'sgp.percentiles.baseline', 'sgp.projections.lagged', 'sgp.projections.lagged.baseline')) {
				tmp.content.areas.label <- 'sgp.content.areas'
				tmp.panel.years.label <- 'sgp.panel.years'
				tmp.grade.sequences.label <- 'sgp.grade.sequences'
			}
			if (!is.null(content_areas)) {
				sgp.config.list[[i]] <- sgp.config.list[[i]][sapply(sgp.config.list[[i]], function(x) tail(x[[tmp.content.areas.label]], 1)) %in% content_areas]
			}
			if (!is.null(years)) {
				sgp.config.list[[i]] <- sgp.config.list[[i]][sapply(sgp.config.list[[i]], function(x) tail(x[[tmp.panel.years.label]], 1)) %in% years]
			}
			if (!is.null(grades)) {
				sgp.config.list[[i]] <- sgp.config.list[[i]][sapply(sgp.config.list[[i]], function(x) tail(x[[tmp.grade.sequences.label]], 1)) %in% grades]
			}
		}
	}


	### Clean up percentile configs for easier reading (don't do for projections - still depend on percentiles elements)

	for (p in grep('sgp.percentiles', names(sgp.config.list))) {
		for (l in seq_along(sgp.config.list[[p]])) {
			sgp.config.list[[p]][[l]] <- sgp.config.list[[p]][[l]][-grep("projection", names(sgp.config.list[[p]][[l]]))]
		}
	}

	### Return result

	return(sgp.config.list)
} ## END getSGPConfig
