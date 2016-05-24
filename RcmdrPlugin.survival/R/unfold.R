# last modified 27 January 2009 by J. Fox

Unfold <- function(){
	# require(survival)
	if (!activeDataSetP()) return()
	initializeDialog(title=gettext("Reshape Wide Survival Data to Long", 
			domain="R-RcmdrPlugin.survival"))
	.activeDataSet <- ActiveDataSet()
	dsname <- tclVar(paste(.activeDataSet, ".long", sep=""))
	dsnameFrame <- tkframe(top)
	entryDsname <- ttkentry(dsnameFrame, width="20", textvariable=dsname)
	nCovSets <- 0
	.CovSets <- list()
	.CovNames <- character(0)
	onOK <- function(){
		if (nCovSets == 0){
			errorCondition(recall=Unfold,
				message=gettext("No time-varying covariates specified.", 
					domain="R-RcmdrPlugin.survival"))
			return()
		}
		dsnameValue <- trim.blanks(tclvalue(dsname))
		if (dsnameValue == "") {
			errorCondition(recall=Unfold,
				message=gettext("You must enter the name of a data set.", 
					domain="R-RcmdrPlugin.survival"))
			return()
		}
		if (!is.valid.name(dsnameValue)) {
			errorCondition(recall=Unfold,
				message=paste('"', dsnameValue, '" ', gettext("is not a valid name.", 
						domain="R-RcmdrPlugin.survival"), sep=""))
			return()
		}
		if (is.element(dsnameValue, listDataSets())) {
			if ("no" == tclvalue(checkReplace(dsnameValue, gettext("Data set", 
						domain="R-RcmdrPlugin.survival")))){
				tkdestroy(top)
				Unfold()
				return()
			}
		}
		time <- getSelection(timeBox)
		if (length(time) == 0){
			errorCondition(recall=Unfold, 
				message=gettext("You must select a time-to-event variable.", domain="R-RcmdrPlugin.survival"))
			return()
		}
		event <- getSelection(eventBox)
		if (length(event) == 0){
			errorCondition(recall=Unfold, message=gettext("You must select an event indicator.", 
					domain="R-RcmdrPlugin.survival"))
			return()
		}
		lag <- tclvalue(lagSliderValue)
		lag <- if (lag == "0") "" else paste(", lag=", lag, sep="")
		closeDialog()
		con <- textConnection("cov", open="w", local=TRUE)
		dump(".CovSets", file=con)
		close(con)
		cov <- paste(cov, collapse="")
		doItAndPrint(cov)
		command <- paste(dsnameValue, " <- unfold(", .activeDataSet,', time="', time, '", event="', event,
			'", cov=.CovSets, cov.names=c(', 
			paste(paste('"', names(.CovSets), '"', sep=""), collapse=","), ')', lag,  ')' , sep="")
		doItAndPrint(command)
		logger("remove(.CovSets)")
		remove(.CovSets, envir=.GlobalEnv)
		tkfocus(CommanderWindow())
		activeDataSet(dsnameValue)
	}
	onCovSelect <- function(){
		covs <- sortVarNames(getSelection(covariateBox))
		if (nCovSets > 0){
			nTimes <- length(.CovSets[[1]])
			if (length(covs) != nTimes) errorCondition(recall=Unfold,
					message=sprintf(gettext("Covariate set has %d entries; should have %d entries", 
							domain="R-RcmdrPlugin.survival"), length(covs), nTimes))
			nCovSets <<- nCovSets + 1
		} else {
			if (length(covs) < 2) errorCondition(recall=Unfold,
					message=gettext("Covariate set must have at least 2 entries.", domain="R-RcmdrPlugin.survival"))
			else nCovSets <<- 1
		}
		name <- trim.blanks(tclvalue(covVariableName))
		if (!is.valid.name(name)){
			errorCondition(recall=Unfold,
				message=paste('"', newVar, '" ',
					gettext("is not a valid name.", domain="R-RcmdrPlugin.survival"), sep=""))
			return()
		}
		if (is.element(name, Variables())) {
			if ("no" == tclvalue(checkReplace(name))){
				tkdestroy(top)
				Unfold()
				return()
			}
		}
		tkconfigure(lagSlider, to=round(length(covs)/4))
		covs <- list(covs)
		names(covs) <- name
		.CovSets <<- c(.CovSets, covs)
		.CovNames <<- c(.CovNames, name)
		tkdelete(newCovBox$listbox, "0", "end")
		for (cov in .CovNames) tkinsert(newCovBox$listbox, "end", cov)
		newCovBox$varlist <<- .CovNames
		tkselection.clear(covariateBox$listbox, "0", "end")
		tclvalue(covVariableName) <- paste("covariate.", nCovSets + 1, sep="")
	}
	OKCancelHelp(helpSubject="Unfold", model=TRUE)
	survFrame <- tkframe(top)
	.activeDataSet <- ActiveDataSet()
	.numeric <- NumericOrDate()
	.factors <- Factors()
	time1 <- eval(parse(text=paste('attr(', .activeDataSet, ', "time1")', sep="")))
	time1 <- if (!is.null(time1)) which(time1 == .numeric) - 1 
	event <- eval(parse(text=paste('attr(', .activeDataSet, ', "event")', sep="")))
	event <- if (!is.null(event)) which(event == Numeric()) - 1 
	timeBox <- variableListBox(survFrame, NumericOrDate(), 
		title=gettext("Time to event\n(select one)", domain="R-RcmdrPlugin.survival"),
		initialSelection=if(is.null(time1)) NULL else time1)
	eventBox <- variableListBox(survFrame, Numeric(), title=gettext("Event indicator\n(select one)", 
			domain="R-RcmdrPlugin.survival"), initialSelection=event)
	covFrame <- tkframe(top)
	covSelectFrame <- tkframe(covFrame)
	covariateBox <- variableListBox(covSelectFrame, Variables(), 
		title=gettext("Select set of\ntime-dependent covariates", domain="R-RcmdrPlugin.survival"),
		selectmode="multiple")
	covSelectButton <- buttonRcmdr(covSelectFrame, 
		text=gettext("Select", domain="R-RcmdrPlugin.survival"), command=onCovSelect)	
	covVariableName <- tclVar("covariate.1")
	newCovFrame <- tkframe(covFrame)
	newCovariate <- ttkentry(newCovFrame, width="20", textvariable=covVariableName)
	newCovBox <- variableListBox(covFrame, c(gettext("<none defined>", domain="R-RcmdrPlugin.survival"), rep("", 4)), 
		title=gettext("Time-dependent covariates", domain="R-RcmdrPlugin.survival"), initialSelection=-1)
	lagSliderValue <- tclVar("0")
	lagSlider <- tkscale(newCovFrame, from=0, to=10,
		showvalue=TRUE, variable=lagSliderValue,
		resolution=1, orient="horizontal")
	tkgrid(labelRcmdr(dsnameFrame, text=gettext("Enter name for data set:", 
				domain="R-RcmdrPlugin.survival")), entryDsname, sticky="w")
	tkgrid(dsnameFrame, sticky="w")
	tkgrid(getFrame(timeBox), labelRcmdr(survFrame, text="  "), getFrame(eventBox), sticky="sw")
	tkgrid(labelRcmdr(survFrame, text=""))
	tkgrid(survFrame, sticky="w")
	tkgrid(labelRcmdr(newCovFrame, text=gettext("Name for covariate", domain="R-RcmdrPlugin.survival"), 
			fg="blue"), sticky="nw")
	tkgrid(newCovariate, sticky="nw")
	tkgrid(labelRcmdr(newCovFrame, text=""))
	tkgrid(labelRcmdr(newCovFrame, text="Lag covariates", fg="blue"), sticky="w")
	tkgrid(lagSlider, sticky="nw")
	tkgrid(getFrame(covariateBox), sticky="nw")
	tkgrid(covSelectButton, sticky="ew")
	tkgrid(covSelectFrame, labelRcmdr(covFrame, text="   "), newCovFrame, labelRcmdr(covFrame, text="   "),
		getFrame(newCovBox), sticky="nw")
	tkgrid(covFrame, sticky="w")
	tkgrid(labelRcmdr(top, text=""))
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=9, columns=1)
}

unfold <- function(data, ...){
	UseMethod("unfold")
}

unfold.data.frame <- function(data, time, event, cov,
	cov.names=paste('covariate', '.', 1:ncovs, sep=""),
	suffix='.time', cov.times=0:ncov, common.times=TRUE, lag=0, 
	show.progress=TRUE, ...){
	# if (show.progress && !require(tcltk)) stop("tcltk package missing")
	vlag <- function(x, lag) c(rep(NA, lag), x[1:(length(x) - lag)])
	xlag <- function(x, lag) apply(as.matrix(x), 2, vlag, lag=lag)
	all.cov <- unlist(cov)
	if (!is.numeric(all.cov)) all.cov <- which(is.element(names(data), all.cov))
	if (!is.list(cov)) cov <- list(cov)
	ncovs <- length(cov)
	nrow <- nrow(data)
	ncol <- ncol(data)
	ncov <- length(cov[[1]])
	nobs <- nrow*ncov
	if (length(unique(c(sapply(cov, length), length(cov.times) - 1))) > 1)
		stop(paste(
				"all elements of cov must be of the same length and \n",
				"cov.times must have one more entry than each element of cov."))
	var.names <- names(data)
	subjects <- rownames(data)
	omit.cols <- if (!common.times) c(all.cov, cov.times) else all.cov
	keep.cols <- (1:ncol)[-omit.cols]
	factors <- names(data)[keep.cols][sapply(data[keep.cols], is.factor)]
	levels <- lapply(data[factors], levels)
	first.covs <- sapply(cov, function(x) x[1])
	factors.covs <- which(sapply(data[first.covs], is.factor))
	levels.covs <- lapply(data[names(factors.covs)], levels)
	nkeep <- length(keep.cols)
	if (is.numeric(event)) event <- var.names[event]
	events <- sort(unique(data[[event]]))
	if (length(events) > 2 || (!is.numeric(events) && !is.logical(events))) 
		stop("event indicator must have values {0, 1}, {1, 2} or {FALSE, TRUE}")
	if (!(all(events == 0:1) || all(events == c(FALSE, TRUE)))){
		if (all(events = 1:2)) data[[event]] <- data[[event]] - 1
		else stop("event indicator must have values {0, 1}, {1, 2} or {FALSE, TRUE}")
	}
	times <- if (common.times) matrix(cov.times, nrow, ncov + 1, byrow=TRUE)
		else data[, cov.times]
	new.data <- matrix(Inf, nobs, 3 + ncovs + nkeep)
	rownames <- rep("", nobs)
	colnames(new.data) <- c('start', 'stop', paste(event, suffix, sep=""),
		var.names[-omit.cols], cov.names)
	end.row <- 0
	if (show.progress){
		progress <- myTkProgressBar(title = "Progress", label = "",
			min = 0, max = 1, initial = 0, width = 300)
		position <- if (is.element("Rcmdr", loadedNamespaces())) 
				paste("+", paste(10 + commanderPosition(), collapse="+"), sep="")
			else "-20+20"
	tkwm.geometry(progress$window, position)
	}
	data <- as.matrix(as.data.frame(lapply(data, as.numeric)))
	for (i in 1:nrow){
		if (show.progress){
			info <- sprintf("%d%% percent done", round(100*i/nrow))
			setTkProgressBar(progress, value=i/nrow, label=info)
		}
		start.row <- end.row + 1
		end.row <- end.row + ncov
		start <- times[i, 1:ncov]
		stop <- times[i, 2:(ncov+1)]
		event.time <- ifelse (stop == data[i, time] & data[i, event] == 1, 1, 0)
#		keep <- matrix(unlist(data[i, -omit.cols]), ncov, nkeep, byrow=TRUE)
		keep <- matrix(data[i, -omit.cols], ncov, nkeep, byrow=TRUE)
		select <- apply(matrix(!is.na(data[i, all.cov]), ncol=ncovs), 1, all)
		rows <- start.row:end.row
#		cov.mat <- xlag(matrix(unlist(data[i, all.cov]), nrow=length(rows)), lag)
		cov.mat <- xlag(matrix(data[i, all.cov], nrow=length(rows)), lag)
		new.data[rows[select], ] <-
			cbind(start, stop, event.time, keep, cov.mat)[select,]
		rownames[rows] <- paste(subjects[i], '.', seq(along=rows), sep="")
	}
	row.names(new.data) <- rownames
	new.data <- as.data.frame(new.data[new.data[, 1] != Inf &
				apply(as.matrix(!is.na(new.data[, cov.names])), 1, all), ])
	for (fac in factors){
		new.data[[fac]] <- factor(levels[[fac]][new.data[[fac]]])
	}
	fcv <- 0
	for (cv in factors.covs){
		fcv <- fcv + 1
		new.data[[cov.names[cv]]] <- factor(levels.covs[[fcv]][new.data[[cov.names[cv]]]])
	}
	attr(new.data, "time1") <- "start"
	attr(new.data, "time2") <- "stop"
	attr(new.data, "event") <- paste(event, suffix, sep="")
	close(progress)
	new.data
}

# the following is a modified version of tkProgressBar() from tcltk:

myTkProgressBar <-
	function (title = "R progress bar", label = "", min = 0, max = 1, 
		initial = 0, width = 300) 
{
	useText <- FALSE
	have_ttk <- as.character(tcl("info", "tclversion")) >= "8.5"
	if (!have_ttk && as.character(tclRequire("PBar")) == "FALSE") 
		useText <- TRUE
	.win <- tktoplevel()
	.val <- initial
	.killed <- FALSE
	tkwm.geometry(.win, sprintf("%dx80", width + 40))
	tkwm.title(.win, title)
#	fn <- tkfont.create(family = "helvetica", size = 12)
	if (useText) {
#		.lab <- tklabel(.win, text = label, font = fn, padx = 20)
		.lab <- tklabel(.win, text = label, padx = 20)
		tkpack(.lab, side = "left")
		fn2 <- tkfont.create(family = "helvetica", size = 16)
		.vlab <- tklabel(.win, text = "0%", font = fn2, padx = 20)
		tkpack(.vlab, side = "right")
		up <- function(value) {
			if (!is.finite(value) || value < min || value > max) 
				return()
			.val <<- value
			tkconfigure(.vlab, text = sprintf("%d%%", round(100 * 
							(value - min)/(max - min))))
		}
	}
	else {
#		.lab <- tklabel(.win, text = label, font = fn, pady = 10)
		.lab <- tklabel(.win, text = label, pady = 10)
		.tkval <- tclVar(0)
		tkpack(.lab, side = "top")
#		tkpack(tklabel(.win, text = "", font = fn), side = "bottom")
		tkpack(tklabel(.win, text = ""), side = "bottom")
		pBar <- if (have_ttk) 
				ttkprogressbar(.win, length = width, variable = .tkval)
			else tkwidget(.win, "ProgressBar", width = width, variable = .tkval)
		tkpack(pBar, side = "bottom")
		up <- function(value) {
			if (!is.finite(value) || value < min || value > max) 
				return()
			.val <<- value
			tclvalue(.tkval) <<- 100 * (value - min)/(max - min)
		}
	}
	getVal <- function() .val
	kill <- function() if (!.killed) {
			tkdestroy(.win)
			.killed <<- TRUE
		}
	title <- function(title) tkwm.title(.win, title)
	lab <- function(label) tkconfigure(.lab, text = label)
	tkbind(.win, "<Destroy>", kill)
	up(initial)
	structure(list(getVal = getVal, up = up, title = title, label = lab, 
			kill = kill, window=.win), class = "tkProgressBar")
}


