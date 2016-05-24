# last modified 2015-08-31 by J. Fox

startStop <- function(time){
	times <- na.omit(eval(parse(text=paste(ActiveDataSet(), '[,c("', time[1], '", "', time[2],'")]', sep=""))))
	if (all(times[[time[1]]] <= times[[time[2]]])){
		return(list(start=time[1], stop=time[2], error=FALSE))
	} else if (all(times[[time[2]]] <= times[[time[1]]])){
		return(list(start=time[2], stop=time[1], error=FALSE))
	}
	else return(list(start="", stop="", error=TRUE))
}

padNA <- function(X, res){
	XX <- matrix(NA, length(res), ncol(X))
	colnames(XX) <- colnames(X)
	XX[!is.na(res), ] <- X
	XX
}

SurvivalData <- function(){
	if (!activeDataSetP()) return()
	initializeDialog(title=gettext("Survival Data Definition", domain="R-RcmdrPlugin.survival"))
	onOK <- function(){
		activeDataSet <- ActiveDataSet()
		time <- getSelection(timeBox)
		if (length(time) == 1){
			time1 <- time
			time2 <- numeric(0)
		}
		else if (length(time) == 2){
			ss <- startStop(time)
			if (ss$error) errorCondition(recall=SurvivalData, 
					message=gettext("Start and stop times must be ordered.", 
						domain="R-RcmdrPlugin.survival"), model=TRUE)
			time1 <- ss$start
			time2 <- ss$stop
		}
		else {
			errorCondition(recall=SurvivalData, message=gettext("You must select one or two time variables.", 
					domain="R-RcmdrPlugin.survival"), model=TRUE)
			return()
		}
		event <- getSelection(eventBox)
		survtype <- as.character(tclvalue(survtypeVariable))
		if (survtype == "interval" && length(event) == 0){
		  errorCondition(recall=SurvivalData, 
		                 message=gettext("You must select an event indicator if censoring is 'interval'.", 
		                                 domain="R-RcmdrPlugin.survival"))
		  return()
		}
		if (survtype == "interval2" && length(event) != 0){
		  errorCondition(recall=SurvivalData, 
		                 message=gettext("You should not select an event indicator if censoring is 'interval2'.", 
		                                 domain="R-RcmdrPlugin.survival"))
		  return()
		}
		if (length(time) == 2 && (! survtype %in% c("interval", "interval2", "counting"))){
		  errorCondition(recall=SurvivalData,
		                 message=gettext("start-end times only for interval or counting-process censoring\nselect Interval, Interval type 2, or Counting process.",
		                                 domain="R-RcmdrPlugin.survival"))
		  return()
		}
		if (length(time) == 1 && survtype %in% c("interval", "interval2", "counting")){
		  errorCondition(recall=SurvivalData,
		                 message=gettext("start-end times required for interval or counting-process censoring.",
		                                 domain="R-RcmdrPlugin.survival"))
		  return()
		}
		strata <- getSelection(strataBox)
		cluster <- getSelection(clusterBox)
		closeDialog()
		command <- paste("attr(", activeDataSet, ', "time1") <- "', time1, '"', sep="")
		doItAndPrint(command)
		if (length(time2) > 0){
			command <- paste("attr(", activeDataSet, ', "time2") <- "', time2, '"', sep="")
			doItAndPrint(command)
		}
		if (length(event) > 0){
  		command <- paste("attr(", activeDataSet, ', "event") <- "', event, '"', sep="")
  		doItAndPrint(command)
		}
		if (length(survtype) > 0 && survtype != "default"){
		  command <- paste("attr(", activeDataSet, ', "survtype") <- "', survtype, '"', sep="")
		  doItAndPrint(command)
		}
		
		if (length(strata) > 0){
			command <- paste("attr(", activeDataSet, ', "strata") <- c(', paste(paste('"', strata, '"', sep=""), collapse=","), ')', sep="")
			doItAndPrint(command)
		}
		if (length(cluster) > 0 && nchar(cluster) > 0){
			command <- paste("attr(", activeDataSet, ', "cluster") <- "', cluster, '"', sep="")
			doItAndPrint(command)
		} 
		tkfocus(CommanderWindow())
	}
	onRefresh <- function(type){
		vars <- if (type == "all") Variables() else Factors()
		tkdelete(clusterBox$listbox, "0", "end")
		for (var in vars) tkinsert(clusterBox$listbox, "end", var)
		clusterBox$varlist <<- vars
		cmd <- paste('options(clusters="', if (type == "all") "all.variables" else "factors.only", '")', sep="")
		doItAndPrint(cmd)
		tkfocus(top)
	}
	OKCancelHelp(helpSubject="SurvivalData")
	survFrame <- tkframe(top)
	timeBox <- variableListBox(survFrame, NumericOrDate(), title=gettext("Time or start/end times\n(select one or two)", 
			domain="R-RcmdrPlugin.survival"), selectmode="multiple")
	eventBox <- variableListBox(survFrame, Variables(), title=gettext("Event indicator\n(select one)", 
			domain="R-RcmdrPlugin.survival"))
	radioButtons(survFrame, name="survtype",
	             buttons=c("default", "right", "left", "interval", "interval2", "counting"),
	             labels=gettext(c("Default", "Right", "Left", "Interval", "Interval type 2", "Counting process")),
	             initialValue="default", title=gettext("Type of Censoring", domain="R-RcmdrPlugin.survival"))
	strataBox <- variableListBox(survFrame, Factors(), title=gettext("Strata\n(select zero or more)", 
			domain="R-RcmdrPlugin.survival"), initialSelection=-1, selectmode="multiple")
	cl.vars <- if (allVarsClusters()) Variables() else Factors()
	len.cl.vars <- length(cl.vars)
	if (len.cl.vars < 5) cl.vars <- c(cl.vars, rep("", 5 - len.cl.vars)) 
	clusterBox <- variableListBox(survFrame, cl.vars, 
		title=gettext("Clusters\n(optional)", domain="R-RcmdrPlugin.survival"), initialSelection=-1)
	radioButtons(survFrame, name="clusterButtons",
		buttons=c("factors", "all"), initialValue=if (allVarsClusters()) "all" else "factors",
		labels=gettext(c("Factors only", "All variables"), domain="R-RcmdrPlugin.survival"), 
		title=gettext("Candidates for clusters", domain="R-RcmdrPlugin.survival"))
	tkbind(factorsButton, "<Button-1>", function() onRefresh("factors"))
	tkbind(allButton, "<Button-1>", function() onRefresh("all"))
	tkgrid(getFrame(timeBox), labelRcmdr(survFrame, text="  "), getFrame(eventBox), 
	       labelRcmdr(survFrame, text="  "), survtypeFrame, sticky="nw")
	tkgrid(labelRcmdr(survFrame, text=""))
	tkgrid(getFrame(strataBox), labelRcmdr(survFrame, text="  "), getFrame(clusterBox), sticky="nw")
	tkgrid(labelRcmdr(survFrame, text=""), labelRcmdr(survFrame, text=""), clusterButtonsFrame, sticky="w")
	tkgrid(survFrame, sticky="w")
	tkgrid(labelRcmdr(top, text=""))
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=12, columns=1)
}

allVarsClusters <- function(){
	opt <- match.arg(getOption("clusters"), c("factors.only", "all.variables"))
	opt == "all.variables"
	
}

#toDate <- function(){
#	dataSet <- activeDataSet()
#	initializeDialog(title=gettext("Date Conversion", domain="R-RcmdrPlugin.survival"))
#	oldVariableFrame <- tkframe(top)
#	variableBox <- variableListBox(oldVariableFrame, Factors(), 
#		title=gettext("Variable to convert", domain="R-RcmdrPlugin.survival"))
#	newVariableFrame <- tkframe(top)
#	newVariableName <- tclVar(gettext("date", domain="R-RcmdrPlugin.survival"))
#	newVariable <- ttkentry(newVariableFrame, width="20", textvariable=newVariableName)
#	dateFormatVar <- tclVar("%Y-%m-%d")
#	formatFrame <- tkframe(oldVariableFrame)
#	dateFormat <- ttkentry(formatFrame, width="20", textvariable=dateFormatVar)
#	radioButtonsFrame <- tkframe(newVariableFrame)
#	dateButton <- tkradiobutton(radioButtonsFrame)
#	DateButton <- tkradiobutton(radioButtonsFrame)
#	tkbind(dateButton, "<Button-1>", function() tclvalue(dateFormatVar) <- "mdy")
#	tkbind(DateButton, "<Button-1>", function() tclvalue(dateFormatVar) <- "%Y-%m-%d")
#	dateValue <- tclVar("Date")
#	tkconfigure(dateButton, variable=dateValue, value="date")
#	tkconfigure(DateButton, variable=dateValue, value="Date")
#	onOK <- function(){
#		x <- getSelection(variableBox)
#		if (length(x) == 0){
#			errorCondition(recall=toDate, message=gettextRcmdr("You must select a variable."))
#			return()
#		}
#		newVar <- trim.blanks(tclvalue(newVariableName))
#		if (!is.valid.name(newVar)){
#			errorCondition(recall=toDate,
#				message=paste('"', newVar, '" ', gettextRcmdr("is not a valid name."), sep=""))
#			return()
#		}
#		if (is.element(newVar, Variables())) {
#			if ("no" == tclvalue(checkReplace(newVar, gettextRcmdr("Variable")))){
#				toDate()
#				return()
#			}
#		}
#		fmt <- trim.blanks(tclvalue(dateFormatVar))
#		whichDate <- tclvalue(dateValue)
#		closeDialog()
#		command <-  if (whichDate == "Date") 
#				paste(dataSet,"$",newVar, " <- as.Date(", dataSet, "$", x, ', format="', fmt, '")', sep="")
#			else paste(dataSet,"$",newVar, " <- as.date(as.character(", dataSet, "$", x, '), order="', fmt, '")', sep="")
#		logger(command)
#		result <- justDoIt(command)
#		if (class(result)[1] !=  "try-error") activeDataSet(dataSet, flushModel=FALSE)
#		tkfocus(CommanderWindow())
#	}
#	OKCancelHelp(helpSubject="as.Date")
#	tkgrid(labelRcmdr(formatFrame, text=gettext("Date format", 
#				domain="R-RcmdrPlugin.survival"), fg="blue"), sticky="nw")
#	tkgrid(dateFormat, sticky="nw")
#	tkgrid(getFrame(variableBox), formatFrame, sticky="nw")
#	tkgrid(oldVariableFrame, sticky="nw")
#	tkgrid(labelRcmdr(newVariableFrame, text=gettext("Name for date variable", 
#				domain="R-RcmdrPlugin.survival"), fg="blue"), 
#		labelRcmdr(newVariableFrame, text="   "),
#		labelRcmdr(newVariableFrame, text=gettext("Class of date variable",
#				domain="R-RcmdrPlugin.survival"), fg="blue"), sticky="nw")
#	tkgrid(labelRcmdr(radioButtonsFrame, text=paste("'Date' ", gettext("object", 
#				domain="R-RcmdrPlugin.survival"), sep="")), DateButton, sticky="nw")
#	tkgrid(labelRcmdr(radioButtonsFrame, text=paste("'date' ", gettext("object", 
#				domain="R-RcmdrPlugin.survival"), sep="")), dateButton, sticky="nw")
#	tkgrid(newVariable, labelRcmdr(newVariableFrame, text="   "), radioButtonsFrame, sticky="nw")
#	tkgrid(newVariableFrame, sticky="nw")
#	tkgrid(buttonsFrame, sticky="w")
#	dialogSuffix(rows=3, columns=1)
#}

# The following version of the toDate dialog converts only to "date" (not "Date") objects
# since "Date" objects are not currently supported in the survival package;
# the remaining infrastucture for "Date" objects is left intact.

toDate <- function(){
	dataSet <- activeDataSet()
	initializeDialog(title=gettext("Date Conversion", domain="R-RcmdrPlugin.survival"))
	oldVariableFrame <- tkframe(top)
	variableBox <- variableListBox(oldVariableFrame, Factors(), 
		title=gettext("Variable to convert", domain="R-RcmdrPlugin.survival"))
	newVariableFrame <- tkframe(top)
	newVariableName <- tclVar(gettext("date", domain="R-RcmdrPlugin.survival"))
	newVariable <- ttkentry(newVariableFrame, width="20", textvariable=newVariableName)
	dateFormatVar <- tclVar("%Y-%m-%d")
	formatFrame <- tkframe(oldVariableFrame)
	dateFormat <- ttkentry(formatFrame, width="20", textvariable=dateFormatVar)
	radioButtonsFrame <- tkframe(newVariableFrame)
	dateButton <- tkradiobutton(radioButtonsFrame)
	DateButton <- tkradiobutton(radioButtonsFrame)
	tkbind(dateButton, "<Button-1>", function() tclvalue(dateFormatVar) <- "mdy")
	tkbind(DateButton, "<Button-1>", function() tclvalue(dateFormatVar) <- "%Y-%m-%d")
	dateValue <- tclVar("Date")
	tkconfigure(dateButton, variable=dateValue, value="date")
	tkconfigure(DateButton, variable=dateValue, value="Date")
	onOK <- function(){
		x <- getSelection(variableBox)
		if (length(x) == 0){
			errorCondition(recall=toDate, message=gettextRcmdr("You must select a variable."))
			return()
		}
		newVar <- trim.blanks(tclvalue(newVariableName))
		if (!is.valid.name(newVar)){
			errorCondition(recall=toDate,
				message=paste('"', newVar, '" ', gettextRcmdr("is not a valid name."), sep=""))
			return()
		}
		if (is.element(newVar, Variables())) {
			if ("no" == tclvalue(checkReplace(newVar, gettextRcmdr("Variable")))){
				toDate()
				return()
			}
		}
		fmt <- trim.blanks(tclvalue(dateFormatVar))
		whichDate <- tclvalue(dateValue)
		closeDialog()
		command <-  if (whichDate == "Date") 
				paste(dataSet,"$",newVar, " <- as.Date(", dataSet, "$", x, ', format="', fmt, '")', sep="")
			else paste(dataSet,"$",newVar, " <- as.date(as.character(", dataSet, "$", x, '), order="', fmt, '")', sep="")
		logger(command)
		result <- justDoIt(command)
		if (class(result)[1] !=  "try-error") activeDataSet(dataSet, flushModel=FALSE)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="as.Date")
	tkgrid(labelRcmdr(formatFrame, text=gettext("Date format", 
				domain="R-RcmdrPlugin.survival"), fg="blue"), sticky="nw")
	tkgrid(dateFormat, sticky="nw")
	tkgrid(getFrame(variableBox), formatFrame, sticky="nw")
	tkgrid(oldVariableFrame, sticky="nw")
	tkgrid(labelRcmdr(newVariableFrame, text=gettext("Name for date variable", 
				domain="R-RcmdrPlugin.survival"), fg="blue"), 
		labelRcmdr(newVariableFrame, text="   "),
		labelRcmdr(newVariableFrame, text=gettext("Class of date variable",
				domain="R-RcmdrPlugin.survival"), fg="blue"), 
		sticky="nw")
	tkgrid(labelRcmdr(radioButtonsFrame, text=paste("'Date' ", gettext("object", 
					domain="R-RcmdrPlugin.survival"), sep="")), DateButton, sticky="nw")
	tkgrid(labelRcmdr(radioButtonsFrame, text=paste("'date' ", gettext("object", 
					domain="R-RcmdrPlugin.survival"), sep="")), dateButton, sticky="nw")
	tkgrid(newVariable, labelRcmdr(newVariableFrame, text="   "), 
		radioButtonsFrame, 
		sticky="nw")
	tkgrid(newVariableFrame, sticky="nw")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=3, columns=1)
}

NumericOrDate <- function(dataSet=ActiveDataSet()) {
	setdiff(Variables(), Factors())
}

# # the following function masks anova.coxph() in the survival package, to change the default
# #   test to "Chisq"
# 
# anova.coxph <- function(object, ..., test="Chisq"){
# 	survival:::anova.coxph(object, ..., test=test)
# }

mfrow <- function (n, max.plots = 0) {
	if (max.plots != 0 & n > max.plots) 
		stop(paste("number of plots =", n, " exceeds maximum =", 
						max.plots))
	rows <- round(sqrt(n))
	cols <- ceiling(n/rows)
	c(rows, cols)
}

