`studentGrowthPercentiles` <-
function(panel.data,         ## REQUIRED
         sgp.labels,         ## REQUIRED
         panel.data.vnames=NULL,
         additional.vnames.to.return=NULL,
         grade.progression,
         content_area.progression=NULL,
         year.progression=NULL,
         year_lags.progression=NULL,
         num.prior,
         max.order.for.percentile=NULL,
         subset.grade,
         percentile.cuts=NULL,
         growth.levels,
         use.my.knots.boundaries,
         use.my.coefficient.matrices=NULL,
         calculate.confidence.intervals=NULL,
         print.other.gp=FALSE,
         print.sgp.order=FALSE,
         calculate.sgps=TRUE,
         rq.method="br",
         rq.method.for.large.n="fn",
         max.n.for.coefficient.matrices=NULL,
         knot.cut.percentiles=c(0.2,0.4,0.6,0.8),
         knots.boundaries.by.panel=FALSE,
         exact.grade.progression.sequence=FALSE,
         drop.nonsequential.grade.progression.variables=TRUE,
         convert.0and100=TRUE,
         sgp.quantiles="Percentiles",
         sgp.quantiles.labels=NULL,
         sgp.loss.hoss.adjustment=NULL,
         sgp.cohort.size=NULL,
         percuts.digits=0,
         isotonize=TRUE,
         convert.using.loss.hoss=TRUE,
         goodness.of.fit=TRUE,
         goodness.of.fit.minimum.n=NULL,
         goodness.of.fit.output.format="GROB",
         return.prior.scale.score=TRUE,
         return.prior.scale.score.standardized=TRUE,
         return.norm.group.identifier=TRUE,
         return.norm.group.scale.scores=NULL,
         return.norm.group.dates=NULL,
         return.panel.data=identical(parent.frame(), .GlobalEnv),
         print.time.taken=TRUE,
         parallel.config=NULL,
         calculate.simex=NULL,
         sgp.percentiles.set.seed=314159,
         sgp.percentiles.equated=NULL,
         SGPt=NULL,
         SGPt.max.time=NULL,
         verbose.output=FALSE) {

	started.at <- proc.time()
	started.date <- date()

	##########################################################
	###
	### Internal utility functions
	###
	##########################################################

	.smooth.bound.iso.row <- function(tmp.dt, iso=isotonize, sgp.loss.hoss.adjustment) {
		X <- NULL
		if (!is.null(sgp.loss.hoss.adjustment)) {
			my.path.knots.boundaries <- get.my.knots.boundaries.path(sgp.labels$my.subject, as.character(sgp.labels$my.year))
			bnd <- eval(parse(text=paste("Knots_Boundaries", my.path.knots.boundaries, "[['loss.hoss_", tmp.last, "']]", sep="")))
			tmp.dt[X < bnd[1], X:=bnd[1]]
			tmp.dt[X > bnd[2], X:=bnd[2]]
		}
		if (iso) setkey(tmp.dt, ID, X)
		return(tmp.dt[['X']])
	}

	.create.path <- function(labels, pieces=c("my.subject", "my.year", "my.extra.label")) {
		sub(' ', '_', toupper(sub('\\.+$', '', paste(unlist(sapply(labels[pieces], as.character)), collapse="."))))
	}

	.get.knots.boundaries <- function(data, by.grade) {
		num.panels <- (dim(data)[2]-1)/2

		if (knots.boundaries.by.panel) {
			tmp.years <- rep(yearIncrement(sgp.labels$my.year, (-num.panels+1):-1), each=dim(data)[1])
		} else {
			tmp.years <- rep(sgp.labels$my.year, dim(data)[1]*(num.panels-1))
		}

		if (by.grade) {
			tmp.grades <- unlist(lapply(data[,2:(2+num.panels-2), with=FALSE], as.character), use.names=FALSE)
		} else {
			tmp.grades <- rep(head(tmp.gp, -1), each=dim(data)[1])
		}

		tmp.stack <- data.table(
			VALID_CASE="VALID_CASE",
			CONTENT_AREA=rep(head(content_area.progression, -1), each=dim(data)[1]),
			GRADE=tmp.grades,
			SCALE_SCORE=unlist(lapply(data[,(2+num.panels):(2+2*num.panels-2), with=FALSE], as.numeric), use.names=FALSE),
			YEAR=tmp.years, key=c("VALID_CASE", "CONTENT_AREA", "GRADE"))

		createKnotsBoundaries(tmp.stack, knot.cut.percentiles)
	}

	.get.panel.data <- function(tmp.data, k, by.grade) {
		str1 <- str2 <- str3 <- NULL
		for (i in 0:k) {
			str1 <- paste(str1, " & !is.na(tmp.data[[", 1+2*num.panels-i, "]])", sep="")
			str2 <- paste(str2, " & tmp.data[[", 1+num.panels-i, "]]=='", rev(as.character(tmp.gp))[i+1], "'", sep="")
			str3 <- c(1+2*num.panels-i, str3)
		}
		if (by.grade) {
			tmp.data[eval(parse(text=paste(substring(str1, 4), str2, sep="")))][, c(1, str3), with=FALSE]
		} else {
			tmp.data[eval(parse(text=substring(str1, 4)))][, c(1, str3), with=FALSE]
		}
	}

	get.my.knots.boundaries.path <- function(content_area, year) {
		if (is.null(sgp.percentiles.equated)) {
			tmp.knots.boundaries.names <-
				names(Knots_Boundaries[[tmp.path.knots.boundaries]])[content_area==sapply(strsplit(names(Knots_Boundaries[[tmp.path.knots.boundaries]]), "[.]"), '[', 1)]
			if (length(tmp.knots.boundaries.names)==0) {
				return(paste("[['", tmp.path.knots.boundaries, "']]", sep=""))
			} else {
				tmp.knots.boundaries.years <- sapply(strsplit(tmp.knots.boundaries.names, "[.]"), function(x) x[2])
				tmp.sum <- sum(year >= sort(tmp.knots.boundaries.years), na.rm=TRUE)
				return(paste("[['", tmp.path.knots.boundaries, "']][['", paste(c(content_area, sort(tmp.knots.boundaries.years)[tmp.sum]), collapse="."), "']]", sep=""))
			}
		} else {
			return(paste("[['", tmp.path.knots.boundaries, "']][['", content_area, ".", sgp.percentiles.equated[['Year']], "']]", sep=""))
		}
	}

	get.prior.cutscore.path <- function(content_area, year) {
		if (is.null(sgp.percentiles.equated)) {
			tmp.cutscores <- grep(content_area, names(SGP::SGPstateData[[goodness.of.fit]][['Achievement']][['Cutscores']]), value=TRUE)
			if (length(tmp.cutscores) > 0) {
				tmp.cutscores.names <- tmp.cutscores[content_area==sapply(strsplit(tmp.cutscores, "[.]"), '[', 1)]
				tmp.cutscores.years <- sapply(strsplit(tmp.cutscores.names, "[.]"), function(x) x[2])
				tmp.sum <- sum(year >= sort(tmp.cutscores.years), na.rm=TRUE)
				return(paste(c(content_area, sort(tmp.cutscores.years)[tmp.sum]), collapse="."))
			} else return(content_area)
		} else {
			return(paste(content_area, sgp.percentiles.equated[['Year']], sep="."))
		}
	}

	.create.coefficient.matrices <- function(data, k, by.grade, max.n.for.coefficient.matrices) {
		tmp.data <- .get.panel.data(data, k, by.grade)
		if (dim(tmp.data)[1]==0) return(NULL)
		if (dim(tmp.data)[1] < sgp.cohort.size) return("Insufficient N")
		if (!is.null(max.n.for.coefficient.matrices) && dim(tmp.data)[1] > max.n.for.coefficient.matrices) tmp.data <- tmp.data[sample(seq.int(dim(tmp.data)[1]), max.n.for.coefficient.matrices)]
        if (is.null(max.n.for.coefficient.matrices) && dim(tmp.data)[1] >= 300000) rq.method <- rq.method.for.large.n
		tmp.num.variables <- dim(tmp.data)[2]
		mod <- character()
		s4Ks <- "Knots=list("
		s4Bs <- "Boundaries=list("
		tmp.gp.iter <- rev(tmp.gp)[2:(k+1)]
		for (i in seq_along(tmp.gp.iter)) {
			my.path.knots.boundaries <- get.my.knots.boundaries.path(rev(content_area.progression)[i+1], yearIncrement(rev(year.progression)[i+1], 0))
			.check.knots.boundaries(names(eval(parse(text=paste("Knots_Boundaries", my.path.knots.boundaries, sep="")))), tmp.gp.iter[i])
			knt <- paste("Knots_Boundaries", my.path.knots.boundaries, "[['knots_", tmp.gp.iter[i], "']]", sep="")
			bnd <- paste("Knots_Boundaries", my.path.knots.boundaries, "[['boundaries_", tmp.gp.iter[i], "']]", sep="")
			mod <- paste(mod, " + bs(tmp.data[[", tmp.num.variables-i, "]], knots=", knt, ", Boundary.knots=", bnd, ")", sep="")
			s4Ks <- paste(s4Ks, "knots_", tmp.gp.iter[i], "=", knt, ",", sep="")
			s4Bs <- paste(s4Bs, "boundaries_", tmp.gp.iter[i], "=", bnd, ",", sep="")
		}
		if (!is.null(SGPt)) {
			tmp.data <- data.table(Panel_Data[,c("ID", "TIME", "TIME_LAG"), with=FALSE], key="ID")[tmp.data][,c(names(tmp.data), "TIME", "TIME_LAG"), with=FALSE]
			mod <- paste(mod, " + I(tmp.data[['TIME']]) + I(tmp.data[['TIME_LAG']])", sep="")
		}

		if (!is.null(tmp.par.config <- parallel.config)) if (is.null(parallel.config[["WORKERS"]][["TAUS"]])) tmp.par.config <- NULL

		if (is.null(tmp.par.config)) {
			tmp.mtx <- eval(parse(text=paste("rq(tmp.data[[", tmp.num.variables, "]] ~ ", substring(mod,4), ", tau=taus, data=tmp.data, method=rq.method)[['coefficients']]", sep="")))
		} else {
			par.start <- startParallel(tmp.par.config, 'TAUS', qr.taus=taus)

			if (toupper(tmp.par.config[["BACKEND"]]) == "FOREACH") {
				# tmp.data <<- tmp.data
				tmp.mtx <- foreach(j = iter(par.start$TAUS.LIST), .export=c("tmp.data", "Knots_Boundaries", "rq.method"), .combine = "cbind", .packages="quantreg",
				.inorder=TRUE, .options.mpi = par.start$foreach.options, .options.multicore = par.start$foreach.options, .options.snow = par.start$foreach.options) %dopar% {
					eval(parse(text=paste("rq(tmp.data[[", tmp.num.variables, "]] ~ ", substring(mod,4), ", tau=j, data=tmp.data, method=rq.method)[['coefficients']]", sep="")))
				}
			} else {
				if (par.start$par.type == 'MULTICORE') {
					tmp.mtx <- mclapply(par.start$TAUS.LIST, function(x) eval(parse(text=paste("rq(tmp.data[[", tmp.num.variables, "]] ~ ",
						substring(mod,4), ", tau=x, data=tmp.data, method=rq.method)[['coefficients']]", sep=""))), mc.cores=par.start$workers, mc.preschedule = FALSE)
					tmp.mtx <- do.call(cbind, tmp.mtx)
				}

				if (par.start$par.type == 'SNOW') {
					tmp.mtx <- parLapplyLB(par.start$internal.cl, par.start$TAUS.LIST, function(x) eval(parse(text=paste("rq(tmp.data[[",
						tmp.num.variables, "]] ~ ", substring(mod,4), ", tau=x, data=tmp.data, method=rq.method)[['coefficients']]", sep=""))))
					tmp.mtx <- do.call(cbind, tmp.mtx)
				}
			}
			stopParallel(tmp.par.config, par.start)
		}

		tmp.version <- list(
			SGP_Package_Version=as.character(packageVersion("SGP")),
			Date_Prepared=date(),
			Matrix_Information=list(
				N=dim(tmp.data)[1],
				Model=paste("rq(tmp.data[[", tmp.num.variables, "]] ~ ", substring(mod,4), ", tau=taus, data=tmp.data, method=rq.method)[['coefficients']]", sep=""),
				SGPt=if (is.null(SGPt)) NULL else list(VARIABLES=unlist(SGPt), MAX_TIME=max(tmp.data$TIME, na.rm=TRUE), MAX_TIME_PRIOR=max(tmp.data$TIME-tmp.data$TIME_LAG, na.rm=TRUE), RANGE_TIME_LAG=range(tmp.data$TIME_LAG))))

		eval(parse(text=paste("new('splineMatrix', tmp.mtx, ", substring(s4Ks, 1, nchar(s4Ks)-1), "), ", substring(s4Bs, 1, nchar(s4Bs)-1), "), ",
			"Content_Areas=list(as.character(tail(content_area.progression, k+1))), ",
			"Grade_Progression=list(as.character(tail(tmp.slot.gp, k+1))), ",
			"Time=list(as.character(tail(year.progression, k+1))), ",
			"Time_Lags=list(as.numeric(tail(year_lags.progression, k))), ",
			"Version=tmp.version)", sep="")))

	} ### END .create.coefficient.matrices

	.check.knots.boundaries <- function(names, grade) {
		tmp <- do.call(rbind, strsplit(names, "_"))
		if (!grade %in% tmp[tmp[,1]=="knots", 2]) stop(paste("knots_", grade, " not found in Knots_Boundaries.", sep=""))
		if (!grade %in% tmp[tmp[,1]=="boundaries", 2]) stop(paste("boundaries_", grade, " not found in Knots_Boundaries.", sep=""))
	}

	.create_taus <- function(sgp.quantiles) {
		if (is.character(sgp.quantiles)) {
			taus <- switch(sgp.quantiles,
				PERCENTILES = (1:100-0.5)/100)
		}
		if (is.numeric(sgp.quantiles)) {
			taus <- sgp.quantiles
		}
		return(taus)
	}

	get.coefficient.matrix.name <- function(tmp.last, k) {
		return(paste("qrmatrix_", tmp.last, "_", k, sep=""))
	}

	.get.percentile.predictions <- function(my.data, my.matrix, SGPt.max.time=NULL) {
		SCORE <- TIME <- TIME_LAG <- NULL
		mod <- character()
		for (k in seq_along(my.matrix@Time_Lags[[1]])) {
			knt <- paste("my.matrix@Knots[[", k, "]]", sep="")
			bnd <- paste("my.matrix@Boundaries[[", k, "]]", sep="")
			mod <- paste(mod, ", bs(my.data[[", dim(my.data)[2]-k, "]], knots=", knt, ", Boundary.knots=", bnd, ")", sep="")
		}
		if (!is.null(SGPt)) {
			my.data <- data.table(Panel_Data[,c("ID", "TIME", "TIME_LAG"), with=FALSE], key="ID")[my.data][,c(names(my.data), "TIME", "TIME_LAG"), with=FALSE]
            tmp.time.shift.index <- getTimeShiftIndex(max(as.numeric(my.data[['TIME']])), my.matrix)
### TEMP TEST CODE
            if (max(my.data$TIME+365*-tmp.time.shift.index, na.rm=TRUE) > my.matrix@Version[['Matrix_Information']][['SGPt']][['MAX_TIME']]+30) stop("Matrix Misfit with TIME data!!!")
            if (is.null(SGPt.max.time) && tmp.time.shift.index != 0) my.data[,TIME:=TIME+365*-tmp.time.shift.index]
            if (!is.null(SGPt.max.time)) {
                my.data[,TIME_LAG:=TIME_LAG+my.matrix@Version[['Matrix_Information']][['SGPt']][['MAX_TIME']]-(TIME+365*-tmp.time.shift.index)]
                my.data[,TIME:=my.matrix@Version[['Matrix_Information']][['SGPt']][['MAX_TIME']]]
            }
			mod <- paste(mod, ", my.data[['TIME']], my.data[['TIME_LAG']]", sep="")
		}
		tmp <- eval(parse(text=paste("cbind(1L, ", substring(mod, 2), ") %*% my.matrix", sep="")))
		return(round(matrix(.smooth.bound.iso.row(data.table(ID=rep(seq.int(dim(tmp)[1]), each=length(taus)), X=as.vector(t(tmp))), isotonize, sgp.loss.hoss.adjustment),
			ncol=length(taus), byrow=TRUE), digits=5))
	}

	.get.quantiles <- function(data1, data2) {
		V1 <- NULL
		tmp <- as.data.table(max.col(cbind(data1 < data2, FALSE), "last"))[V1==101,V1:=0]
		if (!is.null(sgp.quantiles.labels)) {
			setattr(tmp[['V1']] <- as.factor(tmp[['V1']]), "levels", sgp.quantiles.labels)
			return(as.integer(levels(tmp[['V1']]))[tmp[['V1']]])
		} else {
			if (!is.null(sgp.loss.hoss.adjustment)) {
				my.path.knots.boundaries <- get.my.knots.boundaries.path(sgp.labels$my.subject, as.character(sgp.labels$my.year))
				tmp.hoss <- eval(parse(text=paste("Knots_Boundaries", my.path.knots.boundaries, "[['loss.hoss_", tmp.last, "']][2]", sep="")))
				if (length(tmp.index <- which(data2>=tmp.hoss)) > 0) {
					tmp[tmp.index, V1:=apply(data.table(data1 > data2, TRUE)[tmp.index], 1, function(x) which.max(x)-1L)]
				}
			}
			if (convert.0and100) {
				tmp[V1==0, V1:=1L]
				tmp[V1==100, V1:=99L]
			}
			return(tmp[['V1']])
		}
	}

	.get.percentile.cuts <- function(data1) {
		tmp <- round(data1[ , percentile.cuts+1, drop=FALSE], digits=percuts.digits)
		if (convert.using.loss.hoss) {
			my.path.knots.boundaries <- get.my.knots.boundaries.path(sgp.labels$my.subject, as.character(sgp.labels$my.year))
			bnd <- eval(parse(text=paste("Knots_Boundaries", my.path.knots.boundaries, "[['loss.hoss_", tmp.last, "']]", sep="")))
			tmp[tmp < bnd[1]] <- bnd[1]
			tmp[tmp > bnd[2]] <- bnd[2]
		}
		colnames(tmp) <- paste("PERCENTILE_CUT_", percentile.cuts, sep="")
		return(tmp)
	}

    .get.best.cuts <- function(list.of.cuts, label.suffix=NULL) {
        cuts.best <- data.table(rbindlist(list.of.cuts), key="ID")
        cuts.best <- cuts.best[c(which(!duplicated(cuts.best))[-1]-1, nrow(cuts.best))][,-1, with=FALSE]
        if (!is.null(label.suffix)) setnames(cuts.best, names(cuts.best), paste(names(cuts.best), label.suffix, sep="_"))
        return(cuts.best)
    }
	split.location <- function(years) sapply(strsplit(years, '_'), length)[1]

	###
	### SIMEX function
	###
	simex.sgp <- function(
		state, variable=NULL, csem.data.vnames=NULL, csem.loss.hoss=NULL,
		lambda, B, simex.sample.size, extrapolation, save.matrices, simex.use.my.coefficient.matrices=NULL, calculate.simex.sgps, dependent.var.error=FALSE, verbose=FALSE)
	{

		if (is.null(dependent.var.error)) dependent.var.error <- FALSE
		if (is.null(verbose)) verbose <- FALSE
		if (verbose) messageSGP(c("\n\tStarted SIMEX SGP calculation ", rev(content_area.progression)[1], " Grade ", rev(tmp.gp)[1], " ", date()))

		GRADE <- CONTENT_AREA <- YEAR <- V1 <- Lambda <- tau <- b <- .SD <- TEMP <- NULL ## To avoid R CMD check warnings
		my.path.knots.boundaries <- get.my.knots.boundaries.path(sgp.labels$my.subject, as.character(sgp.labels$my.year))
		if (is.logical(simex.use.my.coefficient.matrices) && !simex.use.my.coefficient.matrices) simex.use.my.coefficient.matrices <- NULL
		if (!is.null(state) & !is.null(variable)) stop("SIMEX config can not use both 'state' and 'variable' elements.")
		if (!is.null(state) & !is.null(csem.data.vnames)) stop("SIMEX config can not use both 'state' and 'csem.data.vnames' elements.")
		if (!is.null(csem.data.vnames) & !is.null(variable)) stop("SIMEX config can not use both 'csem.data.vnames' and 'variable' elements.")
		if (!is.null(parallel.config)) {
			if (is.null(parallel.config[["WORKERS"]][["SIMEX"]])) tmp.par.config <- NULL else tmp.par.config <- parallel.config
		} else tmp.par.config <- NULL

		rq.mtx <- function(tmp.gp.iter, lam, rqdata) {
			mod <- character()
			s4Ks <- "Knots=list("
			s4Bs <- "Boundaries=list("
			for (i in seq_along(tmp.gp.iter)) {
				knt <- paste("Knots_Boundaries", my.path.knots.boundaries, "[['Lambda_", lam, "']][['knots_", tmp.gp.iter[i], "']]", sep="")
				bnd <- paste("Knots_Boundaries", my.path.knots.boundaries, "[['Lambda_", lam, "']][['boundaries_", tmp.gp.iter[i], "']]", sep="")
				mod <- paste(mod, " + bs(prior_", i, ", knots=", knt, ", Boundary.knots=", bnd, ")", sep="")
				s4Ks <- paste(s4Ks, "knots_", tmp.gp.iter[i], "=", knt, ",", sep="")
				s4Bs <- paste(s4Bs, "boundaries_", tmp.gp.iter[i], "=", bnd, ",", sep="")
			}
			tmp.mtx <-eval(parse(text=paste("rq(final_yr ~", substring(mod,4), ", tau=taus, data = rqdata, method=rq.method)[['coefficients']]", sep="")))

			tmp.version <- list(
					SGP_Package_Version=as.character(packageVersion("SGP")),
					Date_Prepared=date(),
					Matrix_Information=list(
						N=dim(rqdata)[1],
						Model=paste("rq(tmp.data[[", tmp.num.variables, "]] ~ ", substring(mod,4), ", tau=taus, data=tmp.data, method=rq.method)[['coefficients']]", sep=""),
						SGPt=if (is.null(SGPt)) NULL else list(VARIABLES=unlist(SGPt), MAX_TIME=max(tmp.data$TIME, na.rm=TRUE), MAX_TIME_PRIOR=max(tmp.data$TIME-tmp.data$TIME_LAG, na.rm=TRUE), RANGE_TIME_LAG=range(tmp.data$TIME_LAG))))

			eval(parse(text=paste("new('splineMatrix', tmp.mtx, ", substring(s4Ks, 1, nchar(s4Ks)-1), "), ", substring(s4Bs, 1, nchar(s4Bs)-1), "), ",
				"Content_Areas=list(as.character(tail(content_area.progression, k+1))), ",
				"Grade_Progression=list(as.character(tail(tmp.slot.gp, k+1))), ",
				"Time=list(as.character(tail(year.progression, k+1))), ",
				"Time_Lags=list(as.numeric(tail(year_lags.progression, k))), ",
				"Version=tmp.version)", sep="")))
		}

		fitted <- extrap <- tmp.quantiles.simex <- simex.coef.matrices <- list()
		if (dependent.var.error) {
			loss.hoss <- matrix(nrow=2, ncol=length(tmp.gp))
			lh.ca <- rev(content_area.progression)
			lh.gp <- rev(tmp.gp)
		} else {
			loss.hoss <- matrix(nrow=2, ncol=length(tmp.gp)-1)
			lh.ca <- rev(content_area.progression)[-1]
			lh.gp <- rev(tmp.gp)[-1]
			if (!is.null(csem.data.vnames)) {
				if (length(content_area.progression) == length(csem.data.vnames)) csem.data.vnames <- head(csem.data.vnames, -1)
			}
		}
		if (!is.null(csem.loss.hoss)) {
			if (!is.list(csem.loss.hoss)) stop("SIMEX config element 'csem.loss.hoss' must be a 2 level nested list with LOSS/HOSS data for each subject (level 1) by grade (level 2).")
			for (g in 1:ncol(loss.hoss)) {
				loss.hoss[,g] <- csem.loss.hoss[[lh.ca[g]]][[paste("loss.hoss_", lh.gp[g], sep="")]]
			}}
		if (!is.null(state)) {
			for (g in 1:ncol(loss.hoss)) {
				loss.hoss[,g] <- SGP::SGPstateData[[state]][["Achievement"]][["Knots_Boundaries"]][[lh.ca[g]]][[paste("loss.hoss_", lh.gp[g], sep="")]]
			}}
		if (!is.null(variable)) {
			for (g in 1:ncol(loss.hoss)) {
				loss.hoss[,g] <- variable[[paste("loss.hoss_", lh.gp[g], sep="")]]
			}}

		if (!is.null(use.my.coefficient.matrices)) { # Passed implicitly from studentGrowthPercentiles arguments
			if (exact.grade.progression.sequence) {
				simex.matrix.priors <- num.prior
			} else {
				simex.matrix.priors <- seq(num.prior)
			}
		} else simex.matrix.priors <- coefficient.matrix.priors

		for (k in simex.matrix.priors) {
			tmp.data <- .get.panel.data(ss.data, k, by.grade)
			tmp.num.variables <- dim(tmp.data)[2]
			tmp.gp.iter <- rev(tmp.gp)[2:(k+1)]
			if (dependent.var.error) {
				perturb.var <- rev(tmp.gp)[1:(k+1)]
				start.index <- 1
				num.perturb.vars <- tmp.num.variables+1
			} else {
				perturb.var <- tmp.gp.iter
				start.index <- 2
				num.perturb.vars <- tmp.num.variables
			}
			tmp.ca.iter <- rev(content_area.progression)[start.index:(k+1)]
			tmp.yr.iter <- rev(year.progression)[start.index:(k+1)]
			if (is.null(csem.data.vnames)) {
				csem.int <- matrix(nrow=dim(tmp.data)[1], ncol=length(perturb.var)) # build matrix to store interpolated csem
				colnames(csem.int) <- paste("icsem", perturb.var, tmp.ca.iter, tmp.yr.iter, sep="")
			} else {
				csem.int <- data.table(Panel_Data[,c("ID", intersect(csem.data.vnames, names(Panel_Data))),with=FALSE], key="ID")[ID %in% tmp.data$ID]
				setnames(csem.int, csem.data.vnames, paste("icsem", head(tmp.gp, -1), head(content_area.progression, -1), head(year.progression, -1), sep=""))
			}

			# interpolate csem for all scale scores except that of the last grade
			if (!is.null(state)) {
				for (g in seq_along(perturb.var)) {
					if ("YEAR" %in% names(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["CSEM"]])) {
						CSEM_Data <- subset(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["CSEM"]],
							GRADE==perturb.var[g] & CONTENT_AREA== tmp.ca.iter[g] & YEAR==tmp.yr.iter[g])
					} else {
						CSEM_Data <- subset(SGP::SGPstateData[[state]][["Assessment_Program_Information"]][["CSEM"]],
							GRADE==perturb.var[g] & CONTENT_AREA== tmp.ca.iter[g])
					}
					if (dim(CSEM_Data)[1] == 0) stop(paste('CSEM data for', tmp.ca.iter[g], 'Grade', perturb.var[g], 'is required to use SIMEX functionality, but is not available in SGPstateData.  Please contact package administrators to add CSEM data.'))
					CSEM_Function <- splinefun(CSEM_Data[["SCALE_SCORE"]], CSEM_Data[["SCALE_SCORE_CSEM"]], method="natural")
					csem.int[, paste("icsem", perturb.var[g], tmp.ca.iter[g], tmp.yr.iter[g], sep="")] <- CSEM_Function(tmp.data[[num.perturb.vars-g]])
				}
			}

			if (!is.null(variable)){
				# the following is added by YS on 021915 to accommodate missing data when using "variable"
				csem.tmp <- vector()
				for (g in seq_along(perturb.var)) {
					csem.tmp <- cbind(csem.tmp, variable[[paste("CSEM.grade", perturb.var[g], ".", tmp.ca.iter[g], sep="")]])
				}
				csem.tmp <- na.omit(csem.tmp)

				for (g in seq_along(perturb.var)) {
					csem.int[, paste("icsem", perturb.var[g], tmp.ca.iter[g], tmp.yr.iter[g], sep="")] <- csem.tmp[,g]
				}
			}

			## naive model
			if (calculate.simex.sgps) {
				fitted[[paste("order_", k, sep="")]] <- matrix(0, nrow=length(lambda), ncol=dim(tmp.data)[1]*length(taus))
				tmp.matrix <- getsplineMatrices(
					Coefficient_Matrices[[tmp.path.coefficient.matrices]],
					tail(content_area.progression, k+1),
					tail(grade.progression, k+1),
					tail(year.progression, k+1),
					tail(year_lags.progression, k),
					my.matrix.order=k,
					my.matrix.time.dependency=SGPt)[[1]]

				fitted[[paste("order_", k, sep="")]][1,] <- as.vector(.get.percentile.predictions(tmp.data, tmp.matrix))
			}

			if (verbose) messageSGP(c("\t\t", rev(content_area.progression)[1], " Grade ", rev(tmp.gp)[1], " Order ", k, " Started simulation process ", date()))

			## perturb data
			if (!is.null(csem.data.vnames)) {
				tmp.data <- merge(tmp.data, csem.int, by="ID")
			}
			for (L in lambda[-1]) {
				big.data <- rbindlist(replicate(B, tmp.data, simplify = FALSE))
				# big.data[, Lambda := rep(L, each=dim(tmp.data)[1]*B)]
				big.data[, b := rep(1:B, each=dim(tmp.data)[1])]
				if (dependent.var.error) {
					tmp.names <- "b"
				} else {
					setnames(big.data, tmp.num.variables, "final_yr")
					tmp.names <- c("final_yr", "b")
				}
				for (g in seq_along(perturb.var)) {
					col.index <- num.perturb.vars-g
					if (is.null(csem.data.vnames)) {
						setkeyv(big.data, c(names(big.data)[col.index], tmp.names))
						big.data.uniques <- unique(big.data)
						big.data.uniques.indices <- which(!duplicated(big.data))
						big.data.uniques[, paste("icsem", perturb.var[g], tmp.ca.iter[g], tmp.yr.iter[g], sep="") :=
							rep(csem.int[, paste("icsem", perturb.var[g], tmp.ca.iter[g], tmp.yr.iter[g], sep="")], B)[big.data.uniques.indices]]
					} else {
						setkeyv(big.data, c(names(big.data)[col.index], tmp.names, paste("icsem", perturb.var[g], tmp.ca.iter[g], tmp.yr.iter[g], sep="")))
						big.data.uniques <- unique(big.data)
					}
					big.data.uniques[, TEMP := eval(parse(text=paste("big.data.uniques[[", num.perturb.vars-g, "]]+sqrt(L)*big.data.uniques[['icsem",
						perturb.var[g], tmp.ca.iter[g], tmp.yr.iter[g], "']] * rnorm(dim(big.data.uniques)[1])", sep="")))]
					big.data.uniques[big.data.uniques[[col.index]] < loss.hoss[1,g], col.index := loss.hoss[1,g], with=FALSE]
					big.data.uniques[big.data.uniques[[col.index]] > loss.hoss[2,g], col.index := loss.hoss[2,g], with=FALSE]
					if (is.null(key(big.data.uniques))) setkeyv(big.data.uniques, key(big.data))
					big.data[, num.perturb.vars-g := big.data.uniques[,c(key(big.data), "TEMP"), with=FALSE][big.data][['TEMP']]]

					if (is.null(simex.use.my.coefficient.matrices)) {
						ks <- big.data[, as.list(as.vector(unlist(round(quantile(big.data[[col.index]], probs=knot.cut.percentiles, na.rm=TRUE), digits=3))))] # Knots
						bs <- big.data[, as.list(as.vector(round(extendrange(big.data[[col.index]], f=0.1), digits=3)))] # Boundaries
						lh <- big.data[, as.list(as.vector(round(extendrange(big.data[[col.index]], f=0.0), digits=3)))] # LOSS/HOSS

						eval(parse(text=paste("Knots_Boundaries", my.path.knots.boundaries, "[['Lambda_", L, "']][['knots_", perturb.var[g],
																	"']] <- c(ks[,V1], ks[,V2], ks[,V3], ks[,V4])", sep="")))
						eval(parse(text=paste("Knots_Boundaries", my.path.knots.boundaries, "[['Lambda_", L, "']][['boundaries_", perturb.var[g],
																	"']] <- c(bs[,V1], bs[,V2])", sep="")))
						eval(parse(text=paste("Knots_Boundaries", my.path.knots.boundaries, "[['Lambda_", L, "']][['loss.hoss_", perturb.var[g],
																	"']] <- c(lh[,V1], lh[,V2])", sep="")))
					}

					if (dependent.var.error) setnames(big.data, num.perturb.vars-g, paste("prior_", g-1, sep="")) else setnames(big.data, num.perturb.vars-g, paste("prior_", g, sep=""))
					setkey(big.data, b, ID)
				}
				if (dependent.var.error) setnames(big.data, tmp.num.variables, "final_yr")

				if (!is.null(tmp.par.config)) { # Sequential
				## Write big.data to disk and remove from memory
				dir.create("tmp_data", recursive=TRUE, showWarnings=FALSE)
				if (!exists('year.progression.for.norm.group')) year.progression.for.norm.group <- year.progression # Needed during Baseline Matrix construction
				tmp.dbname <- paste("tmp_data/", paste(tail(paste(year.progression.for.norm.group,
					paste(content_area.progression, grade.progression, sep="_"), sep="_"), num.prior+1), collapse="-"), ".sqlite", sep="")
				con <- dbConnect(SQLite(), dbname = tmp.dbname)
				dbWriteTable(con, name = "simex_data", value=big.data, overwrite=TRUE, row.names=0)
				dbDisconnect(con)
				rm(big.data)
				}

				## Establish the simulation iterations - either 1) 1:B, or 2) a sample of either B or the number of previously computed matrices
				sim.iters <- 1:B

				if (!is.null(simex.use.my.coefficient.matrices)) { # Element from the 'calculate.simex' argument list.
					available.matrices <- unlist(getsplineMatrices(
						Coefficient_Matrices[[paste(tmp.path.coefficient.matrices, '.SIMEX', sep="")]][[
							paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]],
						tail(content_area.progression, k+1),
						tail(grade.progression, k+1),
						tail(year.progression, k+1),
						tail(year_lags.progression, k),
						my.exact.grade.progression.sequence=TRUE,
						return.multiple.matrices=TRUE,
						my.matrix.order=k,
						my.matrix.time.dependency=SGPt), recursive=FALSE)

					if (length(available.matrices) > B) sim.iters <- sample(1:length(available.matrices), B) # Stays as 1:B when length(available.matrices) == B
					if (length(available.matrices) < B) sim.iters <- sample(1:length(available.matrices), B, replace=TRUE)
				}

				if (is.null(tmp.par.config)) { # Sequential
					if (verbose) messageSGP(c("\t\t\tStarted coefficient matrix calculation, Lambda ", L, ": ", date()))
					if (is.null(simex.use.my.coefficient.matrices)) {
						for (z in seq_along(sim.iters)) {
							if (is.null(simex.sample.size) || dim(tmp.data)[1] <= simex.sample.size) {
								simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]][[z]] <-
									rq.mtx(tmp.gp.iter[1:k], lam=L, rqdata=  big.data[b==z][, b:=NULL])
							} else {
								simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]][[z]] <-
									rq.mtx(tmp.gp.iter[1:k], lam=L, rqdata=big.data[b==z][, b:=NULL])
							}
						}
					} else simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]] <- available.matrices[sim.iters]

					if (calculate.simex.sgps) {
						if (verbose) messageSGP(c("\t\t\tStarted percentile prediction calculation, Lambda ", L, ": ", date()))
						for (z in seq_along(sim.iters)) {
							fitted[[paste("order_", k, sep="")]][which(lambda==L),] <- fitted[[paste("order_", k, sep="")]][which(lambda==L),] +
								as.vector(.get.percentile.predictions(big.data[b==z][, b:=NULL],
									simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]][[z]])/B)
						}
					}
				} else {	# Parallel over sim.iters

					###  Always use FOREACH for coefficient matrix production -- need %dorng% to guarantee reproducibility across plateforms (also MUCH more efficient with SNOW/Windows).
					if (toupper(tmp.par.config[["BACKEND"]]) != "FOREACH") tmp.par.config[["BACKEND"]] <- "FOREACH"; tmp.par.config[["TYPE"]] <- "doParallel"

					par.start <- startParallel(tmp.par.config, 'SIMEX')

					## Note, that if you use the parallel.config for SIMEX here, you can also use it for TAUS in the naive analysis
					## Example parallel.config argument: '... parallel.config=list(BACKEND="FOREACH", TYPE="doParallel", WORKERS=list(SIMEX = 4, TAUS = 4))'

					## Calculate coefficient matricies (if needed/requested)
					if (is.null(simex.use.my.coefficient.matrices)) {
						if (verbose) messageSGP(c("\t\t\tStarted coefficient matrix calculation, Lambda ", L, ": ", date()))
						# if (toupper(tmp.par.config[["BACKEND"]]) == "FOREACH") {
							if (is.null(simex.sample.size) || dim(tmp.data)[1] <= simex.sample.size) {
								simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]] <-
									foreach(z=iter(sim.iters), .packages=c("quantreg", "data.table"),
										.export=c("Knots_Boundaries", "rq.method", "taus", "content_area.progression", "tmp.slot.gp", "year.progression", "year_lags.progression", "SGPt"),
										.options.mpi=par.start$foreach.options, .options.multicore=par.start$foreach.options, .options.snow=par.start$foreach.options) %dopar% {
											rq.mtx(tmp.gp.iter[1:k], lam=L, rqdata=dbGetQuery(dbConnect(SQLite(shared.cache = TRUE), dbname = tmp.dbname),
												paste("select * from simex_data where b in ('", z, "')", sep="")))
									}
							} else {
								simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]] <-
									foreach(z=iter(sim.iters), .packages=c("quantreg", "data.table"),
										.export=c("Knots_Boundaries", "rq.method", "taus", "content_area.progression", "tmp.slot.gp", "year.progression", "year_lags.progression", "SGPt"),
										.options.mpi=par.start$foreach.options, .options.multicore=par.start$foreach.options, .options.snow=par.start$foreach.options) %dorng% {
											rq.mtx(tmp.gp.iter[1:k], lam=L, rqdata=dbGetQuery(dbConnect(SQLite(shared.cache = TRUE), dbname = tmp.dbname),
												paste("select * from simex_data where b in ('", z, "')", sep=""))[sample(seq.int(dim(tmp.data)[1]), simex.sample.size),])
									}
							}
						# } else {
						# 	if (par.start$par.type == 'MULTICORE') {
						# 		if (is.null(simex.sample.size) || dim(tmp.data)[1] <= simex.sample.size) {
						# 			simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]] <-
						# 				mclapply(sim.iters, function(z) rq.mtx(tmp.gp.iter[1:k], lam=L, rqdata=dbGetQuery(dbConnect(SQLite(), dbname = tmp.dbname),
						# 					paste("select * from simex_data where b in ('", z, "')", sep=""))), mc.cores=par.start$workers)
						# 		} else {
						# 			simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]] <-
						# 				mclapply(sim.iters, function(z) rq.mtx(tmp.gp.iter[1:k], lam=L, rqdata=dbGetQuery(dbConnect(SQLite(), dbname = tmp.dbname),
						# 					paste("select * from simex_data where b in ('", z, "')", sep=""))[sample(seq.int(dim(tmp.data)[1]), simex.sample.size),]),
						# 					mc.cores=par.start$workers) # mc.preschedule = FALSE, mc.set.seed = FALSE,
						# 		}
						# 	}
						# 	if (par.start$par.type == 'SNOW') {
						# 		if (is.null(simex.sample.size) || dim(tmp.data)[1] <= simex.sample.size) {
						# 			simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]] <-
						# 				parLapply(par.start$internal.cl, sim.iters, function(z) rq.mtx(tmp.gp.iter[1:k], lam=L,
						# 					rqdata=dbGetQuery(dbConnect(SQLite(), dbname = tmp.dbname), paste("select * from simex_data where b in ('", z, "')", sep=""))))
						# 		} else {
						# 			simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]] <-
						# 				parLapply(par.start$internal.cl, sim.iters, function(z) rq.mtx(tmp.gp.iter[1:k], lam=L, rqdata=dbGetQuery(dbConnect(SQLite(), dbname = tmp.dbname),
						# 					paste("select * from simex_data where b in ('", z, "')", sep=""))[sample(seq.int(dim(tmp.data)[1]), simex.sample.size),]))
						# 		}
						# 	}
						# }
					} else {
						simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]] <- available.matrices[sim.iters]
					}

					if (!all(sapply(simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]], is.splineMatrix))) {
						recalc.index <- which(!sapply(simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]], is.splineMatrix))
						messageSGP(c("\n\t\t", rev(content_area.progression)[1], " Grade ", rev(tmp.gp)[1], " Order ", k, " Coefficient Matrix process(es) ", recalc.index, "FAILED!  Attempting to recalculate sequentially..."))
						for (z in recalc.index) {
							if (is.null(simex.sample.size) || dim(tmp.data)[1] <= simex.sample.size) {
								simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]][[z]] <-
									rq.mtx(tmp.gp.iter[1:k], lam=L, rqdata=dbGetQuery(dbConnect(SQLite(shared.cache = TRUE), dbname = tmp.dbname),
										paste("select * from simex_data where b in ('", z, "')", sep="")))
							} else {
								simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]][[z]] <-
									rq.mtx(tmp.gp.iter[1:k], lam=L, rqdata=dbGetQuery(dbConnect(SQLite(shared.cache = TRUE), dbname = tmp.dbname),
										paste("select * from simex_data where b in ('", z, "')", sep=""))[sample(seq.int(dim(tmp.data)[1]), simex.sample.size),])
							}
						}
					}

					## get percentile predictions from coefficient matricies
					if (calculate.simex.sgps) {
						if (verbose) messageSGP(c("\t\t\tStarted percentile prediction calculation, Lambda ", L, ": ", date()))
						# if (toupper(tmp.par.config[["BACKEND"]]) == "FOREACH") {
							mtx.subset <- simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]] # Save on memory copying to R SNOW workers
							environment(.get.percentile.predictions) <- environment()
							environment(.smooth.bound.iso.row) <- environment()
							fitted[[paste("order_", k, sep="")]][which(lambda==L),] <-
								foreach(z=iter(sim.iters), .combine="+", .export=c('tmp.gp', 'taus', 'sgp.loss.hoss.adjustment', 'isotonize', 'SGPt'),
									.options.multicore=par.start$foreach.options) %dopar% { # .options.snow=par.start$foreach.options
										as.vector(.get.percentile.predictions(my.matrix=mtx.subset[[z]], my.data=dbGetQuery(dbConnect(SQLite(), dbname = tmp.dbname),
											paste("select ", paste(c("ID", paste('prior_', k:1, sep=""), "final_yr"), collapse=", "), " from simex_data where b in ('",z,"')", sep="")))/B)
								}
						# } else {
						# 	if (par.start$par.type == 'MULTICORE') {
						# 		tmp.fitted <- mclapply(seq_along(sim.iters), function(z) {
						# 			as.vector(.get.percentile.predictions(dbGetQuery(dbConnect(SQLite(), dbname = tmp.dbname),
						# 				paste("select ", paste(c("ID", paste('prior_', k:1, sep=""), "final_yr"), collapse=", ")," from simex_data where b in ('",z,"')", sep="")),
						# 				simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]][[z]])/B)
						# 		}, mc.cores=par.start$workers)

						# 		fitted[[paste("order_", k, sep="")]][which(lambda==L),] <- Reduce('+', tmp.fitted)

						# 	}
						# 	if (par.start$par.type == 'SNOW') {
						# 		tmp.fitted <- parLapply(par.start$internal.cl, seq_along(sim.iters), function(z) {
						# 			as.vector(.get.percentile.predictions(dbGetQuery(dbConnect(SQLite(), dbname = tmp.dbname),
						# 				paste("select ", paste(c("ID", paste('prior_', k:1, sep=""), "final_yr"), collapse=", ")," from simex_data where b in ('",z,"')", sep="")),
						# 				simex.coef.matrices[[paste("qrmatrices", tail(tmp.gp,1), k, sep="_")]][[paste("lambda_", L, sep="")]][[z]])/B)
						# 		})

						# 		fitted[[paste("order_", k, sep="")]][which(lambda==L),] <- Reduce('+', tmp.fitted)

						# 	}
						# }
					}
					stopParallel(tmp.par.config, par.start)
				}
			} ### END for (L in lambda[-1])
            unlink("tmp_data", recursive=TRUE, force=TRUE)
			if (verbose) messageSGP(c("\t\t", rev(content_area.progression)[1], " Grade ", rev(tmp.gp)[1], " Order ", k, " Simulation process complete ", date()))

			if (calculate.simex.sgps) {
				switch(extrapolation,
					LINEAR = fit <- lm(fitted[[paste("order_", k, sep="")]] ~ lambda),
					QUADRATIC = fit <- lm(fitted[[paste("order_", k, sep="")]] ~ lambda + I(lambda^2)))
				extrap[[paste("order_", k, sep="")]] <-
					matrix(.smooth.bound.iso.row(data.table(ID=seq.int(dim(tmp.data)[1]), X=predict(fit, newdata=data.frame(lambda=-1))[1,]), isotonize, sgp.loss.hoss.adjustment),
						ncol=length(taus), byrow=TRUE)
				tmp.quantiles.simex[[k]] <- data.table(ID=tmp.data[["ID"]], SIMEX_ORDER=k,
					SGP_SIMEX=.get.quantiles(extrap[[paste("order_", k, sep="")]], tmp.data[[tmp.num.variables]]))
			}
		} ### END for (k in simex.matrix.priors)

		if (verbose) messageSGP(c("\tFinished SIMEX SGP calculation ", rev(content_area.progression)[1], " Grade ", rev(tmp.gp)[1], " ", date()))

		if (is.null(save.matrices)) simex.coef.matrices <- NULL
		if (calculate.simex.sgps) {
			quantile.data.simex <- data.table(rbindlist(tmp.quantiles.simex), key=c("ID", "SIMEX_ORDER"))
			setkey(quantile.data.simex, ID) # first key on ID and SIMEX_ORDER, then re-key on ID only to insure sorted order. Don't rely on rbindlist/k ordering...
		} else quantile.data.simex <- data.table("ID"=NA, "SIMEX_ORDER"=NA, "SGP_SIMEX"=NA) # set up empty data.table for ddcast and subsets below.
		if (print.other.gp) {
			quantile.data.simex <- ddcast(quantile.data.simex, ID~SIMEX_ORDER, value.var=setdiff(names(quantile.data.simex), c("ID", "SIMEX_ORDER")), sep="_SIMEX_ORDER_")
			setnames(quantile.data.simex, setdiff(names(quantile.data.simex), c("ID", "SGP_SIMEX", "SIMEX_ORDER")), paste("SGP_SIMEX_ORDER",
				setdiff(names(quantile.data.simex), c("ID", "SGP_SIMEX", "SIMEX_ORDER")), sep="_"))
			return(list(
				DT=data.table(quantile.data.simex,
				SGP_SIMEX=quantile.data.simex[c(which(!duplicated(quantile.data.simex))[-1]-1L, nrow(quantile.data.simex))][["SGP_SIMEX"]]),
				MATRICES = simex.coef.matrices))
		} else {
			if (print.sgp.order | return.norm.group.identifier) {
				return(list(
					DT=quantile.data.simex[c(which(!duplicated(quantile.data.simex))[-1]-1L, nrow(quantile.data.simex))],
					MATRICES=simex.coef.matrices))
			} else {
				return(list(
					DT=quantile.data.simex[c(which(!duplicated(quantile.data.simex))[-1]-1L, nrow(quantile.data.simex)), c("ID", "SGP_SIMEX"), with=FALSE],
					MATRICES=simex.coef.matrices))
			}
		}
		if (verbose) messageSGP(c("\tFinished SIMEX SGP calculation ", rev(content_area.progression)[1], " Grade ", rev(tmp.gp)[1], " ", date(), "\n"))
	} ### END simex.sgp function

	############################################################################
	###
	### Data Preparation & Checks
	###
	############################################################################

	ID <- tmp.messages <- ORDER <- SCALE_SCORE_PRIOR <- NULL

	if (missing(panel.data)) {
		stop("User must supply student achievement data for student growth percentile calculations. NOTE: data is now supplied to function using panel.data argument. See help page for details.")
	}
	if (!(is.matrix(panel.data) | is.list(panel.data))) {
		stop("Supplied panel.data not of a supported class. See help for details of supported classes")
	}
	if (identical(class(panel.data), "list")) {
		if (!("Panel_Data" %in% names(panel.data))) {
			stop("Supplied panel.data missing Panel_Data")
	}
	}
	if (identical(class(panel.data), "list")) {
		if (!is.data.frame(panel.data[["Panel_Data"]]) & !is.data.table(panel.data[["Panel_Data"]])) {
			stop("Supplied panel.data$Panel_Data is not a data.frame or a data.table")
		}
	}
	if (identical(class(panel.data), "list") && !is.null(panel.data[['Coefficient_Matrices']])) {
		panel.data[['Coefficient_Matrices']] <- checksplineMatrix(panel.data[['Coefficient_Matrices']])
	}

	if (!missing(sgp.labels)) {
		if (!is.list(sgp.labels)) {
			stop("Please specify an appropriate list of SGP function labels (sgp.labels). See help page for details.")
	}}
	if (!identical(names(sgp.labels), c("my.year", "my.subject")) &
		!identical(names(sgp.labels), c("my.year", "my.subject", "my.extra.label"))) {
		stop("Please specify an appropriate list for sgp.labels. See help page for details.")
	}
	sgp.labels <- lapply(sgp.labels, toupper)
	tmp.path <- .create.path(sgp.labels)

	if (!missing(growth.levels)) {
		tmp.growth.levels <- list()
		if (!is.list(growth.levels) & !is.character(growth.levels)) {
			tmp.messages <- c(tmp.messages, "\t\tNOTE: growth.levels must be supplied as a list or character abbreviation. See help page for details. studentGrowthPercentiles will be calculated without augmented growth.levels\n")
			tf.growth.levels <- FALSE
		}
		if (is.list(growth.levels)) {
			if (!identical(names(growth.levels), c("my.cuts", "my.levels"))) {
				tmp.messages <- c(tmp.messages, "\t\tNOTE: Please specify an appropriate list for growth.levels. See help page for details. Student growth percentiles will be calculated without augmented growth.levels\n")
				tf.growth.levels <- FALSE
			} else {
				tmp.growth.levels <- growth.levels
				tf.growth.levels <- TRUE
			}
		}
		if (is.character(growth.levels)) {
			if (is.null(SGP::SGPstateData[[growth.levels]][["Growth"]][["Levels"]])) {
				tmp.messages <- c(tmp.messages, "\t\tNOTE: Growth Levels are currently not specified for the indicated state. \n\tPlease contact the SGP package administrator to have your state's data included in the package. Student growth percentiles will be calculated without augmented growth levels\n")
				tf.growth.levels <- FALSE
			} else {
				tmp.growth.levels[["my.cuts"]] <- SGP::SGPstateData[[growth.levels]][["Growth"]][["Cutscores"]][["Cuts"]]
				tmp.growth.levels[["my.levels"]] <- SGP::SGPstateData[[growth.levels]][["Growth"]][["Levels"]]
				tf.growth.levels <- TRUE
			}
		}
	} else {
		tf.growth.levels <- FALSE
	}

	if (!missing(use.my.knots.boundaries)) {
		if (!is.list(use.my.knots.boundaries) & !is.character(use.my.knots.boundaries)) {
			stop("use.my.knots.boundaries must be supplied as a list or character abbreviation. See help page for details.")
		}
		if (is.list(use.my.knots.boundaries)) {
			if (!identical(class(panel.data), "list")) {
				stop("use.my.knots.boundaries is only appropriate when panel data is of class list. See help page for details.")
			}
			if (!identical(names(use.my.knots.boundaries), c("my.year", "my.subject")) &
				!identical(names(use.my.knots.boundaries), c("my.year", "my.subject", "my.extra.label"))) {
					stop("Please specify an appropriate list for use.my.knots.boundaries. See help page for details.")
			}
			tmp.path.knots.boundaries <- .create.path(use.my.knots.boundaries, pieces=c("my.subject", "my.year"))
			if (is.null(panel.data[["Knots_Boundaries"]]) | is.null(panel.data[["Knots_Boundaries"]][[tmp.path.knots.boundaries]])) {
				stop("Knots and Boundaries indicated by use.my.knots.boundaries are not included.")
			}
		}
		if (is.character(use.my.knots.boundaries)) {
			if (is.null(SGP::SGPstateData[[use.my.knots.boundaries]][["Achievement"]][["Knots_Boundaries"]])) {
				tmp.messages <- c(tmp.messages, paste("\t\tNOTE: Knots and Boundaries are currently not implemented for the state indicated (",
				use.my.knots.boundaries, "). Knots and boundaries will be calculated from the data.", "
				Please contact the SGP package administrator to have your Knots and Boundaries included in the package\n", sep=""))
			}
			tmp.path.knots.boundaries <- .create.path(sgp.labels, pieces=c("my.subject", "my.year"))
		}
	} else {
		tmp.path.knots.boundaries <- .create.path(sgp.labels, pieces=c("my.subject", "my.year"))
	}

	if (!is.null(use.my.coefficient.matrices) & !identical(use.my.coefficient.matrices, TRUE)) {
		if (!identical(class(panel.data), "list")) {
			stop("use.my.coefficient.matrices is only appropriate when panel data is of class list. See help page for details.")
		}
		if (!is.list(use.my.coefficient.matrices)) {
			stop("Please specify an appropriate list for argument 'use.my.coefficient.matrices'. See help page for details.")
		}
		if (!identical(names(use.my.coefficient.matrices), c("my.year", "my.subject")) &
			!identical(names(use.my.coefficient.matrices), c("my.year", "my.subject", "my.extra.label"))) {
				stop("Please specify an appropriate list for argument 'use.my.coefficient.matrices'. See help page for details.")
		}
		tmp.path.coefficient.matrices <- .create.path(use.my.coefficient.matrices, pieces=c("my.subject", "my.year"))
		if (is.null(panel.data[["Coefficient_Matrices"]]) | is.null(panel.data[["Coefficient_Matrices"]][[tmp.path.coefficient.matrices]])) {
			stop("Coefficient matrices indicated by argument 'use.my.coefficient.matrices' are not included.")
		}
	} else {
		tmp.path.coefficient.matrices <- tmp.path
	}

	if (is.character(sgp.quantiles)) {
		sgp.quantiles <- toupper(sgp.quantiles)
		if (sgp.quantiles != "PERCENTILES") {
			stop("Character options for sgp.quantiles include only Percentiles at this time. Other options available by specifying a numeric quantity. See help page for details.")
		}
		taus <- .create_taus(sgp.quantiles)
		sgp.quantiles.labels <- NULL
	}
	if (is.numeric(sgp.quantiles)) {
		if (!(all(sgp.quantiles > 0 & sgp.quantiles < 1))) {
			stop("Specify sgp.quantiles as as a vector of probabilities between 0 and 1.")
		}
		taus <- .create_taus(sgp.quantiles)
		if (!is.null(sgp.quantiles.labels)) {
			if (length(sgp.quantiles.labels)!=length(sgp.quantiles)+1) stop("Supplied 'sgp.quantiles.labels' must be 1 longer than supplied 'sgp.quantiles'.")
			if (any(is.na(as.integer(sgp.quantiles.labels)))) stop("Supplied 'sgp.quantiles.labels' must be integer values.")
			sgp.quantiles.labels <- as.integer(sgp.quantiles.labels)
		} else {
			sgp.quantiles.labels <- as.integer(c(100*taus, 100))
		}
	}
	if (!is.null(percentile.cuts)) {
		if (sgp.quantiles != "PERCENTILES") {
			stop("percentile.cuts only appropriate for growth percentiles. Set sgp.quantiles to Percentiles to produce requested percentile.cuts.")
		}
		if (!all(percentile.cuts %in% 0:100)) {
			stop("Specified percentile.cuts must be integers between 0 and 100.")
	}}
	if (!calculate.sgps & (is.character(goodness.of.fit) | goodness.of.fit==TRUE)) {
		tmp.messages <- c(tmp.messages, "\t\tNOTE: Goodness-of-Fit tables only produced when calculating SGPs.\n")
	}
	if (!is.null(calculate.confidence.intervals)) {
		csem.tf <- TRUE
		if (!is.character(calculate.confidence.intervals) & !is.list(calculate.confidence.intervals)) {
			tmp.messages <- c(tmp.messages, "\t\tNOTE: Please supply an appropriate state acronym, variable or list containing details to calculate.confidence.intervals. See help page for details. SGPs will be calculated without confidence intervals.\n")
			csem.tf <- FALSE
		}
		if (is.list(calculate.confidence.intervals)) {
			if (!(("state" %in% names(calculate.confidence.intervals)) | ("variable" %in% names(calculate.confidence.intervals)))) {
				tmp.messages <- c(tmp.messages, "\t\tNOTE: Please specify an appropriate list for calculate.confidence.intervals including state/csem variable, confidence.quantiles, simulation.iterations, distribution and round. See help page for details. SGPs will be calculated without confidence intervals.\n")
				csem.tf <- FALSE
			}
			if ("variable" %in% names(calculate.confidence.intervals) & is.null(panel.data.vnames)) {
				stop("To utilize a supplied CSEM variable for confidence interval calculation you must specify the variables to be used for student growth percentile calculations with the panel.data.vnames argument. See help page for details.")
			}
			if (all(c("state", "variable") %in% names(calculate.confidence.intervals))) {
				stop("Please specify EITHER a state OR a CSEM variable for SGP confidence interval calculation. See help page for details.")
			}
		}
		if (is.character(calculate.confidence.intervals)) {
			if (!calculate.confidence.intervals %in% c(objects(SGP::SGPstateData), names(panel.data[['Panel_Data']]))) {
				tmp.messages <- c(tmp.messages, "\t\tNOTE: Please provide an appropriate state acronym or variable name in supplied data corresponding to CSEMs. See help page for details. SGPs will be calculated without confidence intervals.\n")
				csem.tf <- FALSE
			}
			if (calculate.confidence.intervals %in% objects(SGP::SGPstateData)) {
				if ("YEAR" %in% names(SGP::SGPstateData[[calculate.confidence.intervals]][["Assessment_Program_Information"]][["CSEM"]])) {
					if (!sgp.labels$my.year %in% unique(SGP::SGPstateData[[calculate.confidence.intervals]][["Assessment_Program_Information"]][["CSEM"]][["YEAR"]])) {
						tmp.messages <- c(tmp.messages, "\t\tNOTE: SGPstateData contains year specific CSEMs but year requested is not available. Simulated SGPs and confidence intervals will not be calculated.\n")
						csem.tf <- FALSE
					}
				}
				if (!sgp.labels$my.subject %in% unique(SGP::SGPstateData[[calculate.confidence.intervals]][["Assessment_Program_Information"]][["CSEM"]][["CONTENT_AREA"]])) {
					tmp.messages <- c(tmp.messages, paste("\t\tNOTE: SGPstateData does not contain content area CSEMs for requested content area '", sgp.labels$my.subject, "'. Simulated SGPs and confidence intervals will not be calculated.\n", sep=""))
					csem.tf <- FALSE
				}
				calculate.confidence.intervals <- list(state=calculate.confidence.intervals)
			}
			if (calculate.confidence.intervals %in% names(panel.data[['Panel_Data']])) {
				calculate.confidence.intervals <- list(variable=calculate.confidence.intervals)
			}
		}
		if (sgp.quantiles != "PERCENTILES") {
			tmp.messages <- c(tmp.messages, "\t\tNOTE: When 'sgp.quantiles' is supplied and not equal to PERCENTILES, simulation based standard errors/confidences intervals for SGPs are not available.\n")
			csem.tf <- FALSE
		}
	} else {
		csem.tf <- FALSE
	}

	if (is.logical(calculate.simex) && !calculate.simex) calculate.simex <- NULL # check for calculate.simex=FALSE - same as calculate.simex=NULL
	if (!is.null(calculate.simex)) {
		simex.tf <- TRUE
		if (!is.character(calculate.simex) & !is.list(calculate.simex)) {
			tmp.messages <- c(tmp.messages, "\t\tNOTE: Please supply an appropriate state acronym, variable or list containing details to calculate.simex. See help page for details. SGPs will be calculated without measurement error correction.\n")
			simex.tf <- FALSE
		}
		if (is.list(calculate.simex)) {
			if (!("state" %in% names(calculate.simex)) & !("variable" %in% names(calculate.simex)) & !("csem.data.vnames" %in% names(calculate.simex))) {
				tmp.messages <- c(tmp.messages, "\t\tNOTE: Please specify an appropriate list for calculate.simex including state/csem variable, simulation.iterations, lambda and extrapolation. See help page for details. SGPs will be calculated without measurement error correction.\n")
				simex.tf <- FALSE
			}
			if (all(c("state", "variable") %in% names(calculate.simex))) {
				stop("Please specify EITHER a state OR a CSEM variable for SGP measurement error correction. See help page for details.")
			}
			if (!is.null(calculate.simex$lambda)) {
				if (!is.numeric(calculate.simex$lambda)) {
					tmp.messages <- c(tmp.messages, "\t\tNOTE: Please supply numeric values to lambda. See help page for details. SGPs will be calculated without measurement error correction.\n")
					simex.tf <- FALSE
				}
				if (any(calculate.simex$lambda < 0)) {
					messageSGP("lambda should not contain negative values. Negative values will be ignored.")
					lambda <- calculate.simex$lambda[calculate.simex$lambda >= 0]
				} else lambda=calculate.simex$lambda
				if (is.null(panel.data.vnames) & !is.null(calculate.simex$csem.data.vnames)) stop("Use of csem.data.vnames in SIMEX requires panel.data.vnames be provided.")
			}
		}
		if (is.character(calculate.simex)) {
			if (!calculate.simex %in% c(objects(SGP::SGPstateData), names(panel.data))) {
				tmp.messages <- c(tmp.messages, "\t\tNOTE: Please provide an appropriate state acronym or variable name in supplied data corresponding to CSEMs. See help page for details. SGPs will be calculated without measurement error correction.\n")
				simex.tf <- FALSE
			}
			if (calculate.simex %in% objects(SGP::SGPstateData)) {
				if ("YEAR" %in% names(SGP::SGPstateData[[calculate.simex]][["Assessment_Program_Information"]][["CSEM"]])) {
					if (!sgp.labels$my.year %in% unique(SGP::SGPstateData[[calculate.simex]][["Assessment_Program_Information"]][["CSEM"]][["YEAR"]])) {
						tmp.messages <- c(tmp.messages, "\t\tNOTE: SGPstateData contains year specific CSEMs but year requested is not available. SGPs will be calculated without measurement error correction.\n")
						simex.tf <- FALSE
					}
				}
				if (!sgp.labels$my.subject %in% unique(SGP::SGPstateData[[calculate.simex]][["Assessment_Program_Information"]][["CSEM"]][["CONTENT_AREA"]])) {
					tmp.messages <- c(tmp.messages, paste("\t\tNOTE: SGPstateData does not contain content area CSEMs for requested content area '",
						sgp.labels$my.subject, "'. SGPs will be calculated without measurement error correction.\n", sep=""))
					simex.tf <- FALSE
				}
				calculate.simex <- list(state=calculate.simex)
			}
			if (calculate.simex %in% names(panel.data)) {
				calculate.simex <- list(variable=calculate.simex)
			}
		}
		if (is.null(calculate.simex$simulation.iterations)) calculate.simex$simulation.iterations <- 20
		if (!is.null(calculate.simex$simex.sample.size) && !is.numeric(calculate.simex$simex.sample.size)) calculate.simex$simulation.sample.size <- NULL
		if (is.null(calculate.simex$lambda)) calculate.simex$lambda <- seq(0,2,0.5)
		if (is.null(calculate.simex$extrapolation)) {
			calculate.simex$extrapolation <- "LINEAR"
		} else {
			calculate.simex$extrapolation <- toupper(calculate.simex$extrapolation)
		}
		if (!any(calculate.simex$extrapolation == c("QUADRATIC", "LINEAR", "NATURAL"))) {
			messageSGP("\t\tNOTE: Extrapolation not implemented. Using: linear")
			calculate.simex$extrapolation <- "LINEAR"
		}

	} else {
		simex.tf <- FALSE
	}

	if (!is.null(additional.vnames.to.return)) {
		if (!all(names(additional.vnames.to.return) %in% names(panel.data[["Panel_Data"]]))) {
			tmp.messages <- c(tmp.messages, "\t\tNOTE: Supplied 'additional.vnames.to.return' are not all contained in supplied panel.data. No additional variables will be returned.\n")
			additional.vnames.to.return <- NULL
		}
	}

	if (is.null(sgp.cohort.size)) sgp.cohort.size <- 0

	if (is.null(goodness.of.fit.minimum.n)) goodness.of.fit.minimum.n <- 250

	if (!is.null(sgp.percentiles.set.seed)) set.seed(as.integer(sgp.percentiles.set.seed))

	if (!is.null(SGPt)) {
		if (identical(SGPt, TRUE)) SGPt <- list(TIME="TIME", TIME_LAG="TIME_LAG")
		if (is.list(SGPt) && !all(c("TIME", "TIME_LAG") %in% names(SGPt))) {
			tmp.messages <- c(tmp.messages, "\t\tNOTE: 'TIME' and 'TIME_LAG' are not contained in list supplied to 'SGPt' argument. SGPt is set to NULL")
			SGPt <- NULL
		} else {
			if (!((all(unlist(SGPt) %in% names(panel.data))) | (all(unlist(SGPt) %in% names(panel.data$Panel_Data))))) {
				tmp.messages <- c(tmp.messages, "\t\tNOTE: Variables", paste(unlist(SGPt), collapse=", "), "are not all contained in the supplied 'panel.data'. 'SGPt' is set to NULL.\n")
				SGPt <- NULL
			}
		}
	}

	if (is.null(SGPt) && !is.null(return.norm.group.dates)) {
		return.norm.group.dates <- NULL
	}

    if (!is.null(SGPt.max.time) && is.null(SGPt)) {
        SGPt.max.time <- NULL
    }

	if (identical(return.norm.group.dates, TRUE)) {
		return.norm.group.dates <- "TIME[!_]"
	}

	if (identical(return.norm.group.scale.scores, FALSE)) {
		return.norm.group.scale.scores <- NULL
	}


	### Create object to store the studentGrowthPercentiles objects

	tmp.objects <- c("Coefficient_Matrices", "Cutscores", "Goodness_of_Fit", "Knots_Boundaries", "Panel_Data", "SGPercentiles", "SGProjections", "Simulated_SGPs")
	Coefficient_Matrices <- Cutscores <- Goodness_of_Fit <- Knots_Boundaries <- Panel_Data <- SGPercentiles <- SGProjections <- Simulated_SGPs <- SGP_STANDARD_ERROR <- Verbose_Messages <- NULL
	SGP_SIMEX <- SGP_NORM_GROUP_SCALE_SCORES <- SGP_NORM_GROUP_DATES <- SGP_NORM_GROUP <- NULL

	if (identical(class(panel.data), "list")) {
		for (i in tmp.objects) {
			if (!is.null(panel.data[[i]])) {
				assign(i, panel.data[[i]])
			}
		}

		## Check class and construction of coefficient matrices

		if (!is.null(panel.data[['Coefficient_Matrices']])) {
			tmp.matrices <- Coefficient_Matrices; tmp.changes <- FALSE
			for (i in names(tmp.matrices)) {
				splineMatrix.tf <- sapply(tmp.matrices[[i]], validObject, test=TRUE)==TRUE
				if (!any(splineMatrix.tf)) {
					tmp.changes <- TRUE
					tmp.content_area <- unlist(strsplit(i, "[.]"))[1]; tmp.year <- unlist(strsplit(i, "[.]"))[2]
					for (j in names(tmp.matrices[[i]])[!splineMatrix.tf]) {
						messageSGP(paste("\t\tUpdating Existing Coefficient Matrix", i, j, "to new splineMatrix class."))
						tmp.matrices[[i]][[j]] <- as.splineMatrix(matrix_argument=tmp.matrices[[i]][[j]],
							matrix_argument_name=j, content_area=tmp.content_area, year=tmp.year, sgp_object=panel.data)
					}
				}
			}
			if (tmp.changes) {
				Coefficient_Matrices <- tmp.matrices
			}
		}
	} ### if (identical(class(panel.data), "list"))


	### Create Panel_Data based upon class of input data

	if (is.matrix(panel.data)) {
		Panel_Data <- panel.data <- as.data.table(panel.data)
	}
	if (is.data.frame(panel.data)) {
		Panel_Data <- as.data.table(panel.data)
	}
	if (identical(class(panel.data), "list") && !is.data.table(panel.data[["Panel_Data"]])) {
			Panel_Data <- as.data.table(panel.data[["Panel_Data"]])
	}

	### Create ss.data from Panel_Data

	if (dim(Panel_Data)[1]==0 | dim(Panel_Data)[2]<3) {
		tmp.messages <- paste("\t\tNOTE: Supplied data together with grade progression contains no data (dim = ", paste(dim(Panel_Data), collapse=", "), "). Check data, function arguments and see help page for details.\n", sep="")
		messageSGP(paste("\tStarted studentGrowthPercentiles", started.date))
		messageSGP(paste("\t\tSubject: ", sgp.labels$my.subject, ", Year: ", sgp.labels$my.year, ", Grade Progression: ",
			paste(grade.progression, collapse=", "), " ", sgp.labels$my.extra.label, sep=""))
		messageSGP(paste(tmp.messages, "\tFinished SGP Student Growth Percentile Analysis", date(), "in", convertTime(timetaken(started.at)), "\n"))

		return(
			list(Coefficient_Matrices=Coefficient_Matrices,
				Cutscores=Cutscores,
				Goodness_of_Fit=Goodness_of_Fit,
				Knots_Boundaries=Knots_Boundaries,
				Panel_Data=Panel_Data,
				SGPercentiles=SGPercentiles,
				SGProjections=SGProjections,
				Simulated_SGPs=Simulated_SGPs))
	}

	if (!is.null(SGPt)) {
		setnames(Panel_Data, unlist(SGPt), c("TIME", "TIME_LAG"))
	}

	if (!is.null(panel.data.vnames)) {
		if (!all(panel.data.vnames %in% names(Panel_Data))) {
			tmp.messages <- c(tmp.messages, "\t\tNOTE: Supplied 'panel.data.vnames' are not all in the supplied Panel_Data.\n\t\t\tAnalyses will continue with the intersection names contain in Panel_Data.\n")
		}
		ss.data <- Panel_Data[,intersect(panel.data.vnames, names(Panel_Data)), with=FALSE]
	} else {
		ss.data <- Panel_Data
	}

	if (dim(ss.data)[2] %% 2 != 1) {
		stop(paste("Number of columns of supplied panel data (", dim(ss.data)[2], ") does not conform to data requirements. See help page for details."))
	}

	num.panels <- (dim(ss.data)[2]-1)/2

	### Rename variables in ss.data based upon grade progression

	if (!missing(grade.progression)) {
		tmp.gp <- grade.progression
		by.grade <- TRUE

		if (length(tmp.gp[!is.na(tmp.gp)]) > num.panels) {
			tmp.messages <- c(tmp.messages, paste("\t\tNOTE: Supplied 'grade progression', grade.progression=c(", paste(grade.progression, collapse=","), "), exceeds number of panels (", num.panels, ") in provided data.\n\t\t\tAnalyses will utilize maximum number of priors supplied by the data.\n", sep=""))
		tmp.gp <- tail(grade.progression, num.panels)
	}}
	if (!missing(subset.grade) & missing(grade.progression)) {
		tmp.gp <- (subset.grade-num.panels+1):subset.grade
		by.grade <- TRUE
	}
	if (missing(subset.grade) & missing(grade.progression)) {
		tmp.gp <- 1:num.panels
		by.grade <- FALSE
	}
	if (!missing(num.prior) & !exact.grade.progression.sequence) {
		if (length(num.prior) > 1 | !((num.prior-round(num.prior)) < .Machine$double.eps^0.5) | num.prior <= 0) {
			stop("Specified num.prior not positive integer(s)")
		}
		if (num.prior > length(tmp.gp[!is.na(tmp.gp)])-1) {
			tmp.messages <- c(tmp.messages, paste("\t\tNOTE: Specified argument num.prior (", num.prior, ") exceeds number of panels of data supplied.\n\t\t\tAnalyses will utilize maximum number of priors possible.\n", sep=""))
			num.prior <- length(tmp.gp[!is.na(tmp.gp)])-1
		} else {
			tmp.gp <- grade.progression <- tail(tmp.gp[!is.na(tmp.gp)], num.prior+1)
			if (!is.null(content_area.progression) && length(content_area.progression > num.prior+1)) content_area.progression <- tail(content_area.progression, num.prior+1)

	}} else {
		num.prior <- length(tmp.gp[!is.na(tmp.gp)])-1
	}

	if (exact.grade.progression.sequence){
		tmp.gp <- grade.progression
		by.grade <- TRUE
		num.prior <- length(tmp.gp[!is.na(tmp.gp)])-1
	}

	if (!is.null(max.order.for.percentile)) {
		tmp.gp <- tail(tmp.gp, max.order.for.percentile+1)
		num.prior <- min(num.prior, max.order.for.percentile)
		if (!is.null(content_area.progression)) content_area.progression <- tail(content_area.progression, length(tmp.gp))
		if (!is.null(year.progression)) year.progression <- year.progression.for.norm.group <- tail(year.progression, length(tmp.gp))
	}

	if (is.character(tmp.gp)) {
		tmp.slot.gp <- tmp.gp
		tmp.gp <- tmp.gp[!is.na(tmp.gp)]
	} else {
		tmp.slot.gp <- grade.progression
	}

	if (is.numeric(tmp.gp) & drop.nonsequential.grade.progression.variables && any(diff(tmp.gp) > 1)) {
		ss.data <- ss.data[,c(1, (num.panels+1)-rev(c(1, cumsum(rev(diff(tmp.gp)))+1)-1), (2*num.panels+1)-rev(c(1, cumsum(rev(diff(tmp.gp)))+1)-1)), with=FALSE]
		num.panels <- (dim(ss.data)[2]-1)/2
	}

	## Run this check before the setup of ss.data - otherwise function chokes on negative subscripts
	if (exact.grade.progression.sequence & num.prior > num.panels) {
		tmp.messages <- paste("\t\tNOTE: Supplied data together with EXACT grade progression contains fewer panel years than required. \n\t\t
			Check data, function arguments and see help page for details.\n")
		messageSGP(paste("\tStarted studentGrowthPercentiles", started.date))
		messageSGP(paste("\t\tSubject: ", sgp.labels$my.subject, ", Year: ", sgp.labels$my.year, ", Grade Progression: ",
			paste(tmp.slot.gp, collapse=", "), " ", sgp.labels$my.extra.label, sep=""))
		messageSGP(paste(tmp.messages, "\tStudent Growth Percentile Analysis NOT RUN", date(), "\n"))

		return(
			list(Coefficient_Matrices=Coefficient_Matrices,
				Cutscores=Cutscores,
				Goodness_of_Fit=Goodness_of_Fit,
				Knots_Boundaries=Knots_Boundaries,
				Panel_Data=Panel_Data,
				SGPercentiles=SGPercentiles,
				SGProjections=SGProjections,
				Simulated_SGPs=Simulated_SGPs))
	}

	### Create ss.data

	tmp.last <- tail(tmp.gp, 1)
	ss.data <- data.table(ss.data[,c(1, (1+num.panels-num.prior):(1+num.panels), (1+2*num.panels-num.prior):(1+2*num.panels)), with=FALSE], key=names(ss.data)[1])
	num.panels <- (dim(ss.data)[2]-1)/2
	if (is.factor(ss.data[[1]])) ss.data[[1]] <- as.character(ss.data[[1]])
	if (exact.grade.progression.sequence) tmp.num.prior <- num.prior else tmp.num.prior <- 1

	max.cohort.size <- dim(.get.panel.data(ss.data, tmp.num.prior, by.grade))[1]
	if (max.cohort.size == 0) {
		tmp.messages <- "\t\tNOTE: Supplied data together with grade progression contains no data. Check data, function arguments and see help page for details.\n"
		messageSGP(paste("\tStarted studentGrowthPercentiles", started.date))
		messageSGP(paste("\t\tSubject: ", sgp.labels$my.subject, ", Year: ", sgp.labels$my.year, ", Grade Progression: ",
			paste(tmp.slot.gp, collapse=", "), " ", sgp.labels$my.extra.label, sep=""))
		messageSGP(paste(tmp.messages, "\tFinished SGP Student Growth Percentile Analysis", date(), "in", convertTime(timetaken(started.at)), "\n"))

		return(
			list(Coefficient_Matrices=Coefficient_Matrices,
				Cutscores=Cutscores,
				Goodness_of_Fit=Goodness_of_Fit,
				Knots_Boundaries=Knots_Boundaries,
				Panel_Data=Panel_Data,
				SGPercentiles=SGPercentiles,
				SGProjections=SGProjections,
				Simulated_SGPs=Simulated_SGPs))
	}

	if (max.cohort.size < sgp.cohort.size) {
		tmp.messages <- paste("\t\tNOTE: Supplied data together with grade progression contains fewer than the minimum cohort size.\n\t\tOnly", max.cohort.size,
			"valid cases provided with", sgp.cohort.size, "indicated as minimum cohort N size. Check data, function arguments and see help page for details.\n")
		messageSGP(paste("\tStarted studentGrowthPercentiles", started.date))
		messageSGP(paste("\t\tSubject: ", sgp.labels$my.subject, ", Year: ", sgp.labels$my.year, ", Grade Progression: ",
			paste(tmp.slot.gp, collapse=", "), " ", sgp.labels$my.extra.label, sep=""))
		messageSGP(paste(tmp.messages, "\tStudent Growth Percentile Analysis NOT RUN", date(), "\n"))

		return(
			list(Coefficient_Matrices=Coefficient_Matrices,
				Cutscores=Cutscores,
				Goodness_of_Fit=Goodness_of_Fit,
				Knots_Boundaries=Knots_Boundaries,
				Panel_Data=Panel_Data,
				SGPercentiles=SGPercentiles,
				SGProjections=SGProjections,
				Simulated_SGPs=Simulated_SGPs))
	}

	### PROGRESSION variable creation:

	grade.progression <- tmp.gp
	if (is.null(content_area.progression)) {
		content_area.progression <- rep(sgp.labels$my.subject, length(tmp.gp))
	} else {
		if (!identical(class(content_area.progression), "character")) {
			stop("The 'content_area.progression' vector/argument should be a character vector. See help page for details.")
		}
		if (!identical(tail(content_area.progression, 1), sgp.labels[['my.subject']])) {
			stop("The last element in the 'content_area.progression' vector/argument must be identical to 'my.subject' of the sgp.labels. See help page for details.")
		}
		if (length(content_area.progression) != length(tmp.gp)) {
			tmp.messages <- c(tmp.messages, "\t\tNOTE: The 'content_area.progression' vector/argument does not have the same number of elements as the 'grade.progression' vector/argument.\n\t\t\t'content_area.progression' will be trimmed based upon the length of 'grade.progression'.\n")
			content_area.progression <- tail(content_area.progression, length(tmp.gp))
		}
	}

	if (is.null(year.progression) & is.null(year_lags.progression)) {
		if (is.character(type.convert(as.character(grade.progression), as.is=TRUE))) {
			stop("\tNOTE: Non-numeric grade progressions must be accompanied by arguments 'year.progression' and 'year_lags.progression'")
		} else {
			year.progression <- year.progression.for.norm.group <-
				tail(rev(yearIncrement(sgp.labels[['my.year']], c(0, -cumsum(rev(diff(type.convert(as.character(grade.progression)))))))), length(tmp.gp))
		}
	}

	if (is.null(year.progression) & !is.null(year_lags.progression)) {
		if (!identical(sgp.labels[['my.extra.label']], "BASELINE")) {
			year.progression <- year.progression.for.norm.group <- tail(rev(yearIncrement(sgp.labels[['my.year']], c(0, -cumsum(rev(year_lags.progression))))), length(tmp.gp))
		}
		if (identical(sgp.labels[['my.extra.label']], "BASELINE")) {
			year.progression <- rep("BASELINE", length(tmp.gp))
			year.progression.for.norm.group <- tail(rev(yearIncrement(sgp.labels[['my.year']], c(0, -cumsum(rev(year_lags.progression))))), length(tmp.gp))
		}
		if (!identical(class(year.progression), "character")) {
			stop("year.area.progression should be a character vector. See help page for details.")
		}
		if (!identical(sgp.labels[['my.extra.label']], "BASELINE") & !identical(tail(year.progression, 1), sgp.labels[['my.year']])) {
			stop("The last element in the year.progression must be identical to 'my.year' of the sgp.labels. See help page for details.")
		}
		if (length(year.progression) != length(tmp.gp)) {
			tmp.messages <- c(tmp.messages, "\t\tNOTE: The year.progression vector does not have the same number of elements as the grade.progression vector.\n")
		}
	}

	if (!is.null(year.progression) & is.null(year_lags.progression)) {
		year.progression <- tail(year.progression, length(tmp.gp))
		if (year.progression[1] == "BASELINE") {
			year_lags.progression <- rep(1, length(year.progression)-1)
			year.progression.for.norm.group <- year.progression
		} else {
			year_lags.progression <- diff(as.numeric(sapply(strsplit(year.progression, '_'), '[', split.location(year.progression))))
			year.progression.for.norm.group <- year.progression
		}
	}

	### Create Knots and Boundaries if requested (uses only grades in tmp.gp)

	if (missing(use.my.knots.boundaries)) {
		tmp.knots <- c(Knots_Boundaries[[tmp.path.knots.boundaries]], .get.knots.boundaries(ss.data, by.grade))
		Knots_Boundaries[[tmp.path.knots.boundaries]] <- tmp.knots[!duplicated(names(tmp.knots))]
	} else {
		if (is.character(use.my.knots.boundaries)) {
			if (!is.null(SGP::SGPstateData[[use.my.knots.boundaries]][["Achievement"]][["Knots_Boundaries"]])) {
				for (h in unique(content_area.progression)) {
					for (i in grep(h, names(SGP::SGPstateData[[use.my.knots.boundaries]][["Achievement"]][["Knots_Boundaries"]]), value=TRUE)) {
						Knots_Boundaries[[tmp.path.knots.boundaries]][[i]] <- SGP::SGPstateData[[use.my.knots.boundaries]][["Achievement"]][["Knots_Boundaries"]][[i]]
					}
				}
			}
		}
	}

	### QR Calculations: coefficient matrices are saved/read into/from panel.data[["Coefficient_Matrices"]]

	if (is.null(use.my.coefficient.matrices)) {
		if (exact.grade.progression.sequence) {
			coefficient.matrix.priors <- num.prior
		} else {
			coefficient.matrix.priors <- seq(num.prior)
		}
		for (k in coefficient.matrix.priors) {
			Coefficient_Matrices[[tmp.path.coefficient.matrices]][['TMP_NAME']] <- .create.coefficient.matrices(ss.data, k, by.grade, max.n.for.coefficient.matrices)
			if (identical(Coefficient_Matrices[[tmp.path.coefficient.matrices]][['TMP_NAME']], "Insufficient N")) {
				tmp.messages <- c(tmp.messages, paste("\t\tNOTE: Some grade progressions contain fewer than the minimum cohort size.",
					"\n\t\tOnly analyses with MAX grade progression", paste(rev(rev(tmp.gp)[1:k]), collapse = ', '), "will be produced given", sgp.cohort.size,
					"indicated as minimum cohort N size. \n\t\tCheck data, function arguments and see help page for details.\n"))
				Coefficient_Matrices[[tmp.path.coefficient.matrices]][['TMP_NAME']] <- NULL
				grade.progression <- tmp.gp <- rev(rev(tmp.gp)[1:k])
				# num.prior <- length(tmp.gp[2:k]) # Force lots of warnings (?)
				break
			}
			names(Coefficient_Matrices[[tmp.path.coefficient.matrices]])[length(Coefficient_Matrices[[tmp.path.coefficient.matrices]])] <- get.coefficient.matrix.name(tmp.last, k)

			if (verbose.output) {
				tmp.coefficient.matrix.name <- get.coefficient.matrix.name(tmp.last, k)
				tmp.grade.names <- paste("Grade",
					rev(head(unlist(Coefficient_Matrices[[tmp.path.coefficient.matrices]][[tmp.coefficient.matrix.name]]@Grade_Progression), -1)))
				Verbose_Messages <- paste("\t\tNOTE: Coefficient Matrix ", tmp.coefficient.matrix.name, " created using:", sep="")
				for (l in seq_along(tmp.grade.names)) {
					tmp.knots <- paste(tmp.grade.names[l], Coefficient_Matrices[[tmp.path.coefficient.matrices]][[tmp.coefficient.matrix.name]]@Knots[l])
					tmp.boundaries <- paste(tmp.grade.names[l], Coefficient_Matrices[[tmp.path.coefficient.matrices]][[tmp.coefficient.matrix.name]]@Boundaries[l])
					Verbose_Messages <- c(Verbose_Messages, paste("\n\t\t\tKnots: ", tmp.knots, " and Boundaries: ", tmp.boundaries, ".", sep=""))
				}
			}
		}
	}

	### Calculate SIMEX corrected coefficient matrices and percentiles (if requested)

	if (simex.tf) {
		quantile.data.simex <- simex.sgp(
			state=calculate.simex$state,
			variable=calculate.simex$variable,
			csem.data.vnames=calculate.simex$csem.data.vnames,
			csem.loss.hoss=calculate.simex$csem.loss.hoss,
			lambda=calculate.simex$lambda,
			B=calculate.simex$simulation.iterations,
			simex.sample.size=calculate.simex$simex.sample.size,
			extrapolation=calculate.simex$extrapolation,
			save.matrices=calculate.simex$save.matrices,
			simex.use.my.coefficient.matrices=calculate.simex$simex.use.my.coefficient.matrices,
			calculate.simex.sgps=calculate.sgps,
			dependent.var.error=calculate.simex$dependent.var.error,
			verbose=calculate.simex$verbose)

		if (!is.null(quantile.data.simex[['MATRICES']])) {
			tmp_sgp_1 <- list(Coefficient_Matrices = list(TMP_SIMEX=Coefficient_Matrices[[paste(tmp.path.coefficient.matrices, '.SIMEX', sep="")]]))
			tmp_sgp_2 <- list(Coefficient_Matrices = list(TMP_SIMEX=quantile.data.simex[['MATRICES']]))
			tmp_sgp_combined <- mergeSGP(tmp_sgp_1, tmp_sgp_2)
			Coefficient_Matrices[[paste(tmp.path.coefficient.matrices, '.SIMEX', sep="")]] <- tmp_sgp_combined[["Coefficient_Matrices"]][["TMP_SIMEX"]]
		}
	}

	### Calculate growth percentiles (if requested), percentile cuts (if requested), and simulated confidence intervals (if requested)

	if (calculate.sgps) {

		tmp.matrices <- getsplineMatrices(
					Coefficient_Matrices[[tmp.path.coefficient.matrices]],
					content_area.progression,
					grade.progression,
					year.progression,
					year_lags.progression,
					exact.grade.progression.sequence,
					my.matrix.time.dependency=SGPt)

		tmp.orders <- sapply(tmp.matrices, function(x) length(x@Grade_Progression[[1]])-1)
		max.order <- max(tmp.orders)

		if (max.order < num.prior) {
			tmp.messages <- c(tmp.messages, paste("\t\tNOTE: Requested number of prior scores (num.prior=", num.prior, ") exceeds maximum matrix order (max.order=",
			max.order, "). Only matrices of order up to max.order=", max.order, " will be used.\n", sep=""))
		}
		if (max.order > num.prior) {
			tmp.messages <- c(tmp.messages, paste("\t\tNOTE: Maximum coefficient matrix order (max.order=", max.order, ") exceeds that of specified number of priors,
				(num.prior=", num.prior, "). Only matrices of order up to num.prior=", num.prior, " will be used.\n", sep=""))
			tmp.matrices <- tmp.matrices[tmp.orders <= max.order]
		}


		tmp.quantiles <- tmp.percentile.cuts <- tmp.csem.quantiles <- list()

		for (j in seq_along(tmp.orders)) {
			tmp.data <- .get.panel.data(ss.data, tmp.orders[j], by.grade)
			if (dim(tmp.data)[1] > 0) {
				tmp.matrix <- tmp.matrices[[j]]
				tmp.predictions <- .get.percentile.predictions(tmp.data, tmp.matrix)
				tmp.quantiles[[j]] <- data.table(ID=tmp.data[[1]], ORDER=tmp.orders[j], SGP=.get.quantiles(tmp.predictions, tmp.data[[dim(tmp.data)[2]]]))
				if (csem.tf) {
					if (is.null(calculate.confidence.intervals[['simulation.iterations']])) calculate.confidence.intervals[['simulation.iterations']] <- 100
					if (!is.null(calculate.confidence.intervals[['variable']])) {
						tmp.csem.variable <- Panel_Data[Panel_Data[[1]] %in% ss.data[list(tmp.data[[1]])][[1]]][[calculate.confidence.intervals[['variable']]]]
					} else {
						tmp.csem.variable <- NULL
					}
					if (!is.null(additional.vnames.to.return)) {
						tmp.id.etc <- panel.data[["Panel_Data"]][,c("ID", names(additional.vnames.to.return)), with=FALSE][tmp.data[, names(tmp.data)[1], with=FALSE]]
						setnames(tmp.id.etc, names(additional.vnames.to.return), unlist(additional.vnames.to.return))
					}	else tmp.id.etc <- tmp.data[, names(tmp.data)[1], with=FALSE]

					tmp.csem.quantiles[[j]] <- data.table(
									tmp.id.etc,
									matrix(replicate(calculate.confidence.intervals[['simulation.iterations']],
												.get.quantiles(
													tmp.predictions,
													csemScoreSimulator(
													scale_scores=tmp.data[[dim(tmp.data)[2]]],
													grade=tmp.last,
													content_area=sgp.labels[['my.subject']],
													year=sgp.labels[['my.year']],
													state=calculate.confidence.intervals[['state']],
													variable=tmp.csem.variable,
													distribution=calculate.confidence.intervals[['distribution']],
													round.digits=calculate.confidence.intervals[['round']]))), ncol=calculate.confidence.intervals[['simulation.iterations']]))
					setnames(tmp.csem.quantiles[[j]], paste("V", seq(calculate.confidence.intervals[['simulation.iterations']]), sep=""),
										paste("SGP_SIM", seq(calculate.confidence.intervals[['simulation.iterations']]), sep="_"))
				} ## END CSEM analysis

				if (!is.null(percentile.cuts)) {
                    tmp.percentile.cuts[[paste("ORDER", j, sep="_")]] <- data.table(ID=tmp.data[[1]], .get.percentile.cuts(tmp.predictions))
                    if (!is.null(SGPt.max.time)) tmp.percentile.cuts[[paste("ORDER", j, "MAX_TIME", sep="_")]] <- data.table(ID=tmp.data[[1]], .get.percentile.cuts(.get.percentile.predictions(tmp.data, tmp.matrix, SGPt.max.time)))
				}
				if ((is.character(goodness.of.fit) | goodness.of.fit==TRUE | return.prior.scale.score) & j==1) prior.ss <- tmp.data[[dim(tmp.data)[2]-1]]
				if (exact.grade.progression.sequence & return.prior.scale.score) prior.ss <- tmp.data[[dim(tmp.data)[2]-1]]
			} ### END if (dim(tmp.data)[1] > 0)
		} ## END j loop

		quantile.data <- data.table(rbindlist(tmp.quantiles), key="ID")

		if (print.other.gp) {
			quantile.data <- data.table(ddcast(quantile.data, ID ~ ORDER, value.var=setdiff(names(quantile.data), c("ID", "ORDER"))),
				SGP=quantile.data[c(which(!duplicated(quantile.data))[-1]-1L, nrow(quantile.data))][["SGP"]],
				ORDER=as.integer(quantile.data[c(which(!duplicated(quantile.data))[-1]-1L, nrow(quantile.data))][["ORDER"]]))
			setnames(quantile.data, setdiff(names(quantile.data), c("ID", "SGP", "ORDER")), paste("SGP_ORDER", setdiff(names(quantile.data), c("ID", "SGP", "ORDER")), sep="_"))
		} else {
			if (print.sgp.order | return.norm.group.identifier) {
				quantile.data <- quantile.data[c(which(!duplicated(quantile.data))[-1]-1L, nrow(quantile.data))]
			} else {
				quantile.data <- quantile.data[c(which(!duplicated(quantile.data))[-1]-1L, nrow(quantile.data)), c("ID", "SGP"), with=FALSE]
			}
		}

		quantile.data[,SCALE_SCORE_PRIOR:=prior.ss]

		if (return.prior.scale.score.standardized) {
			SCALE_SCORE_PRIOR_STANDARDIZED <- NULL
			quantile.data[,SCALE_SCORE_PRIOR_STANDARDIZED:=round(as.numeric(scale(prior.ss)), digits=3)]
		}

		if (tf.growth.levels) {
			SGP_LEVEL <- NULL
			quantile.data[, SGP_LEVEL:=factor(findInterval(quantile.data[["SGP"]], tmp.growth.levels[["my.cuts"]]),
				levels=seq(length(tmp.growth.levels[["my.levels"]]))-1, ## Necessary in case the full range of SGPs isn't present
				labels=tmp.growth.levels[["my.levels"]], ordered=TRUE)]
		}

		if (csem.tf) {
			simulation.data <- data.table(rbindlist(tmp.csem.quantiles), key="ID")
			simulation.data <- simulation.data[c(which(!duplicated(simulation.data))[-1]-1, nrow(simulation.data))]

			if (is.character(calculate.confidence.intervals) | is.list(calculate.confidence.intervals)) {
				if (is.null(calculate.confidence.intervals$confidence.quantiles) | identical(toupper(calculate.confidence.intervals$confidence.quantiles), "STANDARD_ERROR")) {
					quantile.data[,SGP_STANDARD_ERROR:=round(apply(simulation.data[, -1, with=FALSE], 1, sd, na.rm=TRUE), digits=2)]
				} else {
					if (!(is.numeric(calculate.confidence.intervals$confidence.quantiles) & all(calculate.confidence.intervals$confidence.quantiles < 1) &
						all(calculate.confidence.intervals$confidence.quantiles > 0))) {
						stop("Argument to 'calculate.confidence.intervals$confidence.quantiles' must be numeric and consist of quantiles.")
					}
					tmp.cq <- data.table(round(t(apply(simulation.data[, -1, with=FALSE], 1, quantile, probs = calculate.confidence.intervals$confidence.quantiles))))
					quantile.data[,paste("SGP_", calculate.confidence.intervals$confidence.quantiles, "_CONFIDENCE_BOUND", sep=""):=tmp.cq]
				}
			}
			if (!is.null(calculate.confidence.intervals$confidence.quantiles)) {
				tmp.cq <- data.table(round(t(apply(simulation.data[, -1, with=FALSE], 1, quantile, probs = calculate.confidence.intervals$confidence.quantiles))))
				quantile.data[,paste("SGP_", calculate.confidence.intervals$confidence.quantiles, "_CONFIDENCE_BOUND", sep=""):=tmp.cq]
			}
			Simulated_SGPs[[tmp.path]] <- rbindlist(list(simulation.data, Simulated_SGPs[[tmp.path]]), fill=TRUE)
		}

		if (simex.tf) quantile.data[, SGP_SIMEX:=quantile.data.simex[['DT']][["SGP_SIMEX"]]]

		if (!is.null(percentile.cuts)){
			quantile.data <- data.table(quantile.data, .get.best.cuts(tmp.percentile.cuts[grep("MAX_TIME", names(tmp.percentile.cuts), invert=TRUE)]))
            if (!is.null(SGPt.max.time)) quantile.data <- data.table(quantile.data, .get.best.cuts(tmp.percentile.cuts[grep("MAX_TIME", names(tmp.percentile.cuts))], "MAX_TIME"))
		}

		if (print.sgp.order | return.norm.group.identifier) {
			if (exact.grade.progression.sequence) {
				norm.groups <- paste(tail(paste(year.progression.for.norm.group, paste(content_area.progression, grade.progression, sep="_"), sep="/"), num.prior+1), collapse="; ")
			} else {
				norm.groups <- sapply(seq_along(year.progression.for.norm.group)[-1][1:(num.panels-1)],
				function(x) paste(tail(paste(year.progression.for.norm.group, paste(content_area.progression, grade.progression, sep="_"), sep="/"), x), collapse="; "))
			}
			if (!print.sgp.order) { # Return only SGP_NORM_GROUP
				if (exact.grade.progression.sequence) {
					quantile.data[, SGP_NORM_GROUP:=factor(factor(ORDER, labels=norm.groups))]
				} else {
					quantile.data[, SGP_NORM_GROUP:=factor(factor(ORDER, levels=seq_along(norm.groups), labels=norm.groups))]
				}
				quantile.data[, ORDER:=NULL]
			} else { # Return both ORDER and SGP_NORM_GROUP
				if (exact.grade.progression.sequence) {
					quantile.data[, SGP_NORM_GROUP:=factor(factor(ORDER, labels=norm.groups))]
				} else {
					quantile.data[, SGP_NORM_GROUP:=factor(factor(ORDER, levels=seq_along(norm.groups), labels=norm.groups))]
				}
				setnames(quantile.data, "ORDER", "SGP_ORDER")
			}
		}

		if (!is.null(return.norm.group.dates)) {
			my.tmp <- data.table(Panel_Data[,c("ID", setdiff(grep("TIME", names(Panel_Data), value=TRUE), grep("TIME_LAG", names(Panel_Data), value=TRUE))), with=FALSE],
				key="ID")[list(quantile.data$ID),-1,with=FALSE]
            my.tmp <- my.tmp[,tail(seq(dim(my.tmp)[2]), length(tmp.gp)),with=FALSE]
			quantile.data[,SGP_NORM_GROUP_DATES:=gsub("NA; ", "", do.call(paste, c(as.data.table(lapply(my.tmp, function(x) as.Date(x, origin="1970-01-01"))), list(sep="; "))))]
		}

		if (!is.null(return.norm.group.scale.scores)) {
			my.tmp <- data.table(ss.data[,c("ID", names(tmp.data)[-1]), with=FALSE], key="ID")[list(quantile.data$ID),-1,with=FALSE]
			quantile.data[,SGP_NORM_GROUP_SCALE_SCORES:=gsub("NA; ", "", do.call(paste, c(my.tmp, list(sep="; "))))]
		}

		if ((is.character(goodness.of.fit) | goodness.of.fit==TRUE) & dim(quantile.data)[1] <= goodness.of.fit.minimum.n) {
			messageSGP(c("\tNOTE: Due to small number of cases (", dim(quantile.data)[1], ") no goodness of fit plots produced."))
			goodness.of.fit <- FALSE
		}

		if (is.character(goodness.of.fit) | goodness.of.fit==TRUE) {
			if (simex.tf) {
				sgps.for.gof <- c("SGP", "SGP_SIMEX")
				sgps.for.gof.path <- c(tmp.path, paste(tmp.path, "SIMEX", sep="."))
			} else {
				sgps.for.gof <- "SGP"
				sgps.for.gof.path <- tmp.path
			}
			if (is.character(goodness.of.fit) & goodness.of.fit %in% objects(SGP::SGPstateData) &&
                !identical(sgp.labels$my.extra.label, "EQUATED") &&
				!is.null(SGP::SGPstateData[[goodness.of.fit]][['Achievement']][['Cutscores']][[get.prior.cutscore.path(rev(content_area.progression)[2], yearIncrement(rev(year.progression)[2], 1, year_lags.progression[1]))]][[paste("GRADE_", rev(tmp.gp)[2], sep="")]])) {
				GRADE <- YEAR <- CONTENT_AREA <- NULL
				tmp.gof.data <- getAchievementLevel(
							sgp_data=data.table(
								ID=quantile.data[['ID']],
								SCALE_SCORE=quantile.data[['SCALE_SCORE_PRIOR']],
								quantile.data[, c(sgps.for.gof, "SGP_NORM_GROUP"), with=FALSE],
								VALID_CASE="VALID_CASE",
								CONTENT_AREA=rev(content_area.progression)[2],
								YEAR=rev(year.progression.for.norm.group)[2],
								GRADE=rev(tmp.gp)[2],
								CONTENT_AREA_CURRENT=sgp.labels[['my.subject']],
								YEAR_CURRENT=sgp.labels[['my.year']],
								GRADE_CURRENT=tmp.last),
							state=goodness.of.fit,
							year=rev(year.progression.for.norm.group)[2],
							content_area=rev(content_area.progression)[2],
							grade=tail(tmp.gp, 2)[1])[,!"GRADE", with=FALSE]

				setnames(tmp.gof.data, c("SCALE_SCORE", "ACHIEVEMENT_LEVEL", "CONTENT_AREA", "CONTENT_AREA_CURRENT", "YEAR", "YEAR_CURRENT", "GRADE_CURRENT"),
					c("SCALE_SCORE_PRIOR", "ACHIEVEMENT_LEVEL_PRIOR", "CONTENT_AREA_PRIOR", "CONTENT_AREA", "YEAR_PRIOR", "YEAR", "GRADE"))
				setnames(ss.data, dim(ss.data)[2], "SCALE_SCORE")
				setkeyv(tmp.gof.data, "ID")
				tmp.gof.data <- ss.data[, c("ID", "SCALE_SCORE"), with=FALSE][tmp.gof.data]

				### Rename SGP_NORM_GROUP_BASELINE for gofSGP - expecting consistent name to establish norm.group.var in that function
				if ("SGP_NORM_GROUP_BASELINE" %in% names(tmp.gof.data)) setnames(tmp.gof.data, "SGP_NORM_GROUP_BASELINE", "SGP_NORM_GROUP")

				for (gof.iter in seq_along(sgps.for.gof)) {
					Goodness_of_Fit[[sgps.for.gof.path[gof.iter]]][['TMP_NAME']] <- gofSGP(
						sgp_object=tmp.gof.data,
						state=goodness.of.fit,
						years=sgp.labels[['my.year']],
						content_areas=sgp.labels[['my.subject']],
						content_areas_prior=tmp.gof.data[['CONTENT_AREA_PRIOR']][1],
						grades=tmp.last,
						use.sgp=sgps.for.gof[gof.iter],
						output.format=goodness.of.fit.output.format)

                    if (!is.null(Goodness_of_Fit[[sgps.for.gof.path[gof.iter]]][['TMP_NAME']])) {
                        tmp.gof.plot.name <-
                            paste(tail(paste(year.progression.for.norm.group, paste(content_area.progression, grade.progression, sep="_"), sep="/"), num.prior+1), collapse="; ")
                        tmp.gof.plot.name <- gsub("MATHEMATICS", "MATH", tmp.gof.plot.name)
                        names(Goodness_of_Fit[[sgps.for.gof.path[gof.iter]]])[length(Goodness_of_Fit[[tmp.path]])] <-
                            gsub("/", "_", paste(gsub(";", "", rev(unlist(strsplit(tail(tmp.gof.plot.name, 1), " ")))), collapse=";"))
                    }
				}
			} else {
				tmp.gof.data <- data.table(
					ID=quantile.data[['ID']],
					SCALE_SCORE_PRIOR=quantile.data[['SCALE_SCORE_PRIOR']],
					quantile.data[, sgps.for.gof, with=FALSE],
					VALID_CASE="VALID_CASE",
					CONTENT_AREA=sgp.labels[['my.subject']],
					YEAR=sgp.labels[['my.year']],
					GRADE=tmp.last, key="ID")

				setnames(ss.data, dim(ss.data)[2], "SCALE_SCORE")
				tmp.gof.data <- ss.data[, c("ID", "SCALE_SCORE"), with=FALSE][tmp.gof.data]

				for (gof.iter in seq_along(sgps.for.gof)) {
					Goodness_of_Fit[[sgps.for.gof.path[gof.iter]]][['TMP_NAME']] <- gofSGP(
						sgp_object=tmp.gof.data,
						years=sgp.labels[['my.year']],
						content_areas=sgp.labels[['my.subject']],
						grades=tmp.last,
						use.sgp=sgps.for.gof[gof.iter],
						output.format=goodness.of.fit.output.format)

                    if (!is.null(Goodness_of_Fit[[sgps.for.gof.path[gof.iter]]][['TMP_NAME']])) {
                        tmp.gof.plot.name <-
                            paste(tail(paste(year.progression.for.norm.group, paste(content_area.progression, grade.progression, sep="_"), sep="/"), num.prior+1), collapse="; ")
                        tmp.gof.plot.name <- gsub("MATHEMATICS", "MATH", tmp.gof.plot.name)
                        names(Goodness_of_Fit[[sgps.for.gof.path[gof.iter]]])[length(Goodness_of_Fit[[tmp.path]])] <-
                            gsub("/", "_", paste(gsub(";", "", rev(unlist(strsplit(tail(tmp.gof.plot.name, 1), " ")))), collapse=";"))
                    }
				}
			}
		}

		if (identical(sgp.labels[['my.extra.label']], "BASELINE")) setnames(quantile.data, "SGP", "SGP_BASELINE")
		if (identical(sgp.labels[['my.extra.label']], "BASELINE") & tf.growth.levels) setnames(quantile.data, "SGP_LEVEL", "SGP_LEVEL_BASELINE")
		if (identical(sgp.labels[['my.extra.label']], "BASELINE")) setnames(quantile.data, gsub("SGP_STANDARD_ERROR", "SGP_BASELINE_STANDARD_ERROR", names(quantile.data)))
		if (identical(sgp.labels[['my.extra.label']], "BASELINE")) setnames(quantile.data, gsub("SGP_ORDER", "SGP_BASELINE_ORDER", names(quantile.data)))
		if (identical(sgp.labels[['my.extra.label']], "BASELINE")) setnames(quantile.data, gsub("SGP_NORM_GROUP", "SGP_NORM_GROUP_BASELINE", names(quantile.data)))
		if (identical(sgp.labels[['my.extra.label']], "BASELINE") & simex.tf) setnames(quantile.data, "SGP_SIMEX", "SGP_SIMEX_BASELINE")
		if (identical(sgp.labels[['my.extra.label']], "EQUATED")) setnames(quantile.data, "SGP", "SGP_EQUATED")
		if (identical(sgp.labels[['my.extra.label']], "EQUATED") & tf.growth.levels) setnames(quantile.data, "SGP_LEVEL", "SGP_LEVEL_EQUATED")
		if (identical(sgp.labels[['my.extra.label']], "EQUATED")) setnames(quantile.data, gsub("SGP_NORM_GROUP", "SGP_NORM_GROUP_EQUATED", names(quantile.data)))

		if (!is.null(additional.vnames.to.return)) {
			quantile.data <- data.table(panel.data[["Panel_Data"]][,c("ID", names(additional.vnames.to.return)), with=FALSE], key="ID")[quantile.data]
			setnames(quantile.data, names(additional.vnames.to.return), unlist(additional.vnames.to.return))
		}

		if (!return.prior.scale.score) {
			quantile.data[,SCALE_SCORE_PRIOR:=NULL]
		}

		SGPercentiles[[tmp.path]] <- rbindlist(list(quantile.data, SGPercentiles[[tmp.path]]), fill=TRUE)

	} ## End if calculate.sgps


	### Start/Finish Message & Return SGP Object

	if (print.time.taken) {
		messageSGP(paste("\tStarted studentGrowthPercentiles:", started.date))
		if (calculate.sgps) {
			messageSGP(paste("\t\tContent Area: ", sgp.labels$my.subject, ", Year: ", sgp.labels$my.year, ", Grade Progression: ",
				paste(tmp.slot.gp, collapse=", "), " ", sgp.labels$my.extra.label, " (N=", format(dim(quantile.data)[1], big.mark=","), ")", sep=""))
		} else {
			messageSGP(paste("\t\tContent Area: ", sgp.labels$my.subject, ", Year: ", sgp.labels$my.year, ", Grade Progression: ",
				paste(tmp.slot.gp, collapse=", "), " ", sgp.labels$my.extra.label, " (N=", format(max.cohort.size, big.mark=","), ")", sep=""))
		}
		if (verbose.output) messageSGP(Verbose_Messages)
		messageSGP(c(tmp.messages, "\tFinished SGP Student Growth Percentile Analysis: ", date(), " in ", convertTime(timetaken(started.at)), "\n"))
	}

	list(Coefficient_Matrices=Coefficient_Matrices,
		Cutscores=Cutscores,
		Goodness_of_Fit=Goodness_of_Fit,
		Knots_Boundaries=Knots_Boundaries,
		Panel_Data = if (return.panel.data) Panel_Data else NULL,
		SGPercentiles=SGPercentiles,
		SGProjections=SGProjections,
		Simulated_SGPs=Simulated_SGPs)
} ### END studentGrowthPercentiles Function
