# last modified 2011-08-03 by J. Fox

Survdiff <- function(){
	# require(survival)
	defaults <- list(time1=NULL, event=NULL, strata=NULL, rho="0", subset=NULL)
	dialog.values <- getDialog("Survdiff", defaults)
	if (!activeDataSetP()) return()
	currentModel <- FALSE
	initializeDialog(title=gettext("Compare Survival Functions", domain="R-RcmdrPlugin.survival"))
	onOK <- function(){
		time <- getSelection(timeBox)
		if (length(time) == 1) time1 <- time
		else {
			errorCondition(recall=Survdiff, message=gettext("You must select a time-to-event variable.", 
							domain="R-RcmdrPlugin.survival"))
			return()
		}
		event <- getSelection(eventBox)
		if (length(event) == 0) {
			errorCondition(recall=Survdiff, message=gettext("You must select an event indicator.", 
							domain="R-RcmdrPlugin.survival"))
			return()
		}
		strata <- getSelection(strataBox) 
		if (length(strata) == 0) {
			errorCondition(recall=Survdiff, message=gettext("You must select strata.", 
							domain="R-RcmdrPlugin.survival"))
			return()
		}
		rho <- tclvalue(rhoValue)
		subset <- tclvalue(subsetVariable)
		putDialog("Survdiff", list(
						time1=time1,
						event=event, strata=strata, rho=rho, subset=subset
				))
		closeDialog()
		if (trim.blanks(subset) == gettext("<all valid cases>", domain="R-RcmdrPlugin.survival") 
				|| trim.blanks(subset) == ""){
			subset <- ""
		}
		else{
			subset <- paste(", subset=", subset, sep="")
		}
		formula <- paste("Surv(", time1, ",", event, ")", sep="")
		formula <- paste(formula, " ~ ", paste(strata, collapse=" + "), sep="")
		command <- paste("survdiff(", formula, ", rho=", rho,
				', data=', ActiveDataSet(), subset, ")", sep="")
		doItAndPrint(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="survdiff", reset="Survdiff")
	survFrame <- tkframe(top)
	.activeDataSet <- ActiveDataSet()
	.numeric <- NumericOrDate()
	.factors <- Factors()
	time1 <- if(!is.null(dialog.values$time1)) dialog.values$time1 else eval(parse(text=paste('attr(', .activeDataSet, ', "time1")', sep="")))
	time1 <- if (!is.null(time1)) which(time1 == .numeric) - 1 
	event <- if(!is.null(dialog.values$event)) dialog.values$event else eval(parse(text=paste('attr(', .activeDataSet, ', "event")', sep="")))
	event <- if (!is.null(event)) which(event == Numeric()) - 1 
	strata <- if(!is.null(dialog.values$strata)) dialog.values$strata else eval(parse(text=paste('attr(', .activeDataSet, ', "strata")', sep="")))
	strata <- if (!is.null(strata)) which(is.element(.factors, strata)) - 1 else -1
	timeBox <- variableListBox(survFrame, NumericOrDate(), 
			title=gettext("Time to event\n(select one)", domain="R-RcmdrPlugin.survival"),
			initialSelection=time1)
	eventBox <- variableListBox(survFrame, Numeric(), 
			title=gettext("Event indicator\n(select one)", domain="R-RcmdrPlugin.survival"),
			initialSelection=event)
	strataBox <- variableListBox(survFrame, Factors(), 
			title=gettext("Strata\n(select one or more)", domain="R-RcmdrPlugin.survival"), 
			selectmode="multiple", initialSelection=strata)
	rhoFrame <- tkframe(top)
	rhoValue <- tclVar(dialog.values$rho)
	rhoSlider <- tkscale(rhoFrame, from=0, to=1, showvalue=TRUE, variable=rhoValue,
			resolution=0.1, orient="horizontal")
#	modelFormula(hasLhs=FALSE)
	subsetBox(subset.expression=dialog.values$subset)
	tkgrid(getFrame(timeBox), labelRcmdr(survFrame, text="  "), getFrame(eventBox), sticky="sw")
	tkgrid(labelRcmdr(survFrame, text=""))
	tkgrid(getFrame(strataBox), sticky="nw")
	tkgrid(survFrame, sticky="nw")
	tkgrid(labelRcmdr(rhoFrame, text="rho", foreground="blue"), rhoSlider, sticky="sw")
	tkgrid(rhoFrame, sticky="nw")
	tkgrid(labelRcmdr(top, text=""))
	tkgrid(subsetFrame, sticky="w")
	tkgrid(labelRcmdr(top, text=""))
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=9, columns=1)
}

