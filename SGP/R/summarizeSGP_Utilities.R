`sgpSummary` <- 
function(sgp.groups.to.summarize,
	produce.confidence.interval,
	tmp.simulation.dt,
	state,
	sgp.summaries,
	confidence.interval.groups,
	my.sgp,
	sgp_key,
	variables.for.summaries,
	sim.info) {

	WEIGHT <- MEDIAN_SGP_with_SHRINKAGE <- NULL

	tmp.sgp.summaries <- sgp.summaries
	sgp.summaries.names <- unlist(strsplit(names(sgp.summaries), "[.]"))
	if (produce.confidence.interval) {
		if ("Bootstrap_CI" %in% confidence.interval.groups$TYPE) {
			tmp.list <- list()
			tmp.quantiles <- paste("c(", paste(confidence.interval.groups$QUANTILES, collapse=", "), ")", sep="")
			for (i in confidence.interval.groups$VARIABLES) {
				tmp.list[[paste("MEDIAN_", i, "_QUANTILES", sep="")]] <- paste("boot.sgp(", i, ", ", tmp.quantiles, ")", sep="")
			}
			tmp.sgp.summaries <- c(tmp.sgp.summaries, tmp.list)
			sgp.summaries.names <- c(sgp.summaries.names, do.call(paste, c(data.table(expand.grid("MEDIAN", my.sgp, confidence.interval.groups$QUANTILES, "CONFIDENCE_BOUND_BOOTSTRAP"), key="Var2"), sep="_")))
		} 
		if ("Bootstrap_SE" %in% confidence.interval.groups$TYPE) {
			tmp.list <- list()
			for (i in confidence.interval.groups$VARIABLES) {
				tmp.list[[paste("MEDIAN_", i, "_SE", sep="")]] <- paste("boot.sgp(", i, ")", sep="")
			}
			tmp.sgp.summaries <- c(tmp.sgp.summaries, tmp.list)
			sgp.summaries.names <- c(sgp.summaries.names, do.call(paste, c(data.table(expand.grid("MEDIAN", my.sgp, "STANDARD_ERROR_BOOTSTRAP"), key="Var2"), sep="_"))) 
		}
	}
 
	ListExpr <- parse(text=paste("as.list(c(", paste(unlist(tmp.sgp.summaries), collapse=", "),"))",sep="")) 
	ByExpr <- parse(text=paste("list(", paste(sgp.groups.to.summarize, collapse=", "), ")", sep=""))
	
	pull.vars <- c(unlist(sapply(dbListFields(dbConnect(SQLite(), dbname = "Data/tmp_data/TMP_Summary_Data.sqlite"), "summary_data"), 
		function(p) if (any(grepl(p, tmp.sgp.summaries))) return(p)), use.names=FALSE), strsplit(sgp.groups.to.summarize, ", ")[[1]])
	
	tmp <- pullData(tmp.simulation.dt, state, pull.vars, variables.for.summaries, sgp.groups.to.summarize, sgp_key)[, eval(ListExpr), keyby=eval(ByExpr)]
	setnames(tmp, (dim(tmp)[2]-length(sgp.summaries.names)+1):dim(tmp)[2], sgp.summaries.names)

	if (produce.confidence.interval & "CSEM" %in% confidence.interval.groups[['TYPE']]) {
		pull.vars <- c(sgp_key, unlist(strsplit(sgp.groups.to.summarize, ", ")))
		tmp <- pullData(tmp.simulation.dt, state, pull.vars, variables.for.summaries, sgp.groups.to.summarize, sgp_key, tmp_key = key(tmp), sim.info=sim.info)[tmp]
		setcolorder(tmp, c(grep("CSEM", names(tmp), invert=TRUE), grep("CSEM", names(tmp))))
	}

	if ('MEDIAN_SGP_STANDARD_ERROR' %in% names(tmp)) {
		constant <- var(tmp[['MEDIAN_SGP']], na.rm=TRUE) - mean(tmp[['MEDIAN_SGP_STANDARD_ERROR']]^2, na.rm=TRUE)
		tmp[,MEDIAN_SGP_with_SHRINKAGE := round(50 + ((tmp[['MEDIAN_SGP']]-50) * (constant/(constant+tmp[['MEDIAN_SGP_STANDARD_ERROR']]^2))))]
	}

	messageSGP(paste("\tFinished with", sgp.groups.to.summarize))
	return(tmp)
} ### END sgpSummary function


`pullData` <- 
function(tmp.simulation.dt,
	state,
	pull.vars,
	variables.for.summaries,
	sgp.groups.to.summarize,
	sgp_key,
	tmp_key,
	sim.info=NULL) {

	SGP_SIM <- V1 <- V2 <- SIM_NUM <- WEIGHT <- ACHIEVEMENT_LEVEL <- ACHIEVEMENT_LEVEL_PRIOR <- CATCH_UP_KEEP_UP_STATUS <- MOVE_UP_STAY_UP_STATUS <- NULL
	CATCH_UP_KEEP_UP_STATUS_BASELINE <- MOVE_UP_STAY_UP_STATUS_BASELINE <- NULL

	con <- dbConnect(SQLite(), dbname = "Data/tmp_data/TMP_Summary_Data.sqlite")

	if (!is.null(sim.info)) {
		tmp.list.1 <- list()
		tmp_data <- data.table(dbGetQuery(con, paste("select", paste(pull.vars, collapse = ","), "from summary_data")), key = sgp_key)
		if (is.data.frame(tmp.simulation.dt)) {
			tmp.list.1 <- lapply(seq.int(sim.info[['n.simulated.sgps']]), function(i) {
					tmp_data[,c(key(tmp_data), unlist(strsplit(sgp.groups.to.summarize, ", "))), with=FALSE][
					tmp.simulation.dt[seq.int(i, length.out=sim.info[['n.unique.cases']], by=sim.info[['n.simulated.sgps']])], allow.cartesian=TRUE][,
					list(median_na(SGP_SIM, NULL), mean_na(SGP_SIM, NULL)), keyby=c(unlist(strsplit(sgp.groups.to.summarize, ", ")), "BASELINE")]})
		} else {
			tmp.list.1 <- lapply(seq.int(sim.info[['n.simulated.sgps']]), function(i) {
					tmp_data[data.table(dbGetQuery(con, paste("select * from sim_data where SIM_NUM =", i)), key = sgp_key), allow.cartesian=TRUE][,
					list(median_na(SGP_SIM, NULL), mean_na(SGP_SIM, NULL)), keyby=c(unlist(strsplit(sgp.groups.to.summarize, ", ")), "BASELINE")]})
			dbDisconnect(con)
		}

		tmp.csem <- data.table(reshape(rbindlist(tmp.list.1)[,list(sd_na(V1), sd_na(V2)), keyby=c(unlist(strsplit(sgp.groups.to.summarize, ", ")), "BASELINE")],
			idvar=c(unlist(strsplit(sgp.groups.to.summarize, ", "))), timevar="BASELINE", direction="wide"), key = tmp_key)
		if (length(grep("BASELINE", names(tmp.csem)))==0) {
			setnames(tmp.csem, c("V1.COHORT", "V2.COHORT"), c("MEDIAN_SGP_STANDARD_ERROR_CSEM", "MEAN_SGP_STANDARD_ERROR_CSEM"))
		} else {
			setnames(tmp.csem, c("V1.COHORT", "V2.COHORT", "V1.BASELINE", "V2.BASELINE"), 
				c("MEDIAN_SGP_STANDARD_ERROR_CSEM", "MEAN_SGP_STANDARD_ERROR_CSEM", "MEDIAN_SGP_BASELINE_STANDARD_ERROR_CSEM", "MEAN_SGP_BASELINE_STANDARD_ERROR_CSEM")) 
		}
		return(tmp.csem)
	}
	
	tmp_data <- data.table(dbGetQuery(con, paste("select", paste(pull.vars, collapse = ","), "from summary_data"))) 
	if (all((my.key <- intersect(sgp_key, variables.for.summaries)) %in% names(tmp_data))) setkeyv(tmp_data, my.key)
	if ("CATCH_UP_KEEP_UP_STATUS" %in% names(tmp_data)) {
		tmp_data[, CATCH_UP_KEEP_UP_STATUS := factor(CATCH_UP_KEEP_UP_STATUS)]
	}
	if ("CATCH_UP_KEEP_UP_STATUS_BASELINE" %in% names(tmp_data)) {
		tmp_data[, CATCH_UP_KEEP_UP_STATUS_BASELINE := factor(CATCH_UP_KEEP_UP_STATUS_BASELINE)]
	}
	if ("MOVE_UP_STAY_UP_STATUS" %in% names(tmp_data)) {
		tmp_data[, MOVE_UP_STAY_UP_STATUS := factor(MOVE_UP_STAY_UP_STATUS)]
	}
	if ("MOVE_UP_STAY_UP_STATUS_BASELINE" %in% names(tmp_data)) {
		tmp_data[, MOVE_UP_STAY_UP_STATUS_BASELINE := factor(MOVE_UP_STAY_UP_STATUS_BASELINE)]
	}
	dbDisconnect(con)
	return(tmp_data)
} ### END pullData function
 

`median_na` <- 
function(x,
	weight) {
	if (is.null(weight)) {
		median(as.numeric(x), na.rm=TRUE)
	} else {
		weightedMedian(as.numeric(x), w=weight, na.rm=TRUE)
	}
}


`mean_na` <- 
function(x,
	weight,
	result.digits=2) {

	if (is.null(weight)) {
		tmp.x <- x[!is.na(x)]
		round(sum(tmp.x)/length(tmp.x), digits=result.digits)
	} else {
		round(weighted.mean(as.numeric(x), w=weight, na.rm=TRUE), digits=result.digits)
	}
}


`sd_na` <- function(x, result.digits=2) round(sd(as.numeric(x), na.rm=TRUE), digits=result.digits)


`num_non_missing` <- function(x) length(x[!is.na(x)])


`sgp_standard_error` <- function(x,y=1) round(y*sd(x, na.rm=TRUE)/sqrt(sum(!is.na(as.numeric(x)))), digits=2)


`percent_in_category` <- 
function(x,
	in.categories,
	of.categories,
	result.digits=1) {

	if (!is.list(in.categories)) in.categories <- list(in.categories)
	if (!is.list(of.categories)) of.categories <- list(of.categories)
	tmp <- table(x[!is.na(x)])
	return(unlist(lapply(seq_along(in.categories), function(i) round(100*sum(tmp[in.categories[[i]]], na.rm=TRUE)/sum(tmp[of.categories[[i]]], na.rm=TRUE), digits=result.digits))))
} ### END percent_in_category function


`percent_at_above_target` <- 
function(sgp,
	target,
	result.digits=1) {

	tmp.logical <- sgp >= target
	tmp.pct <- round(sum(tmp.logical, na.rm=TRUE)/sum(!is.na(tmp.logical))*100, digits=result.digits)
	return(tmp.pct)
} ### END percent_at_above_target function


`boot.sgp` <- 
function(dat,
	conf.quantiles=NULL,
	nboot=100) {

	ID <- SCORE <- NULL
	CI <- c(NA,NA); SE <- NA
	if (!all(is.na(dat))) {
		dat.no.na <- dat[!is.na(dat)]
		out <- data.table(ID=rep(seq.int(nboot)), SCORE=dat.no.na[sample.int(length(dat.no.na), length(dat.no.na)*nboot, replace=TRUE)])[,median(SCORE), by=ID][['V1']]
		if (!is.null(conf.quantiles)) CI <- round(quantile(out, conf.quantiles), digits=1) else SE <- round(sd(out), digits=1)
	}
	if (!is.null(conf.quantiles)) return(CI) else return(SE)
} ### END boot.sgp function
