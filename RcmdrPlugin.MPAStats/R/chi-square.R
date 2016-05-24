# Modified on November 26, 2013 by Jordan Gressel

# Interpretation Functions
chiSquareWords <- function(x,var1,var2){
    wrapper <- function(text){
        text2 <- strwrap(text)
        for(i in 1:length(text2)){
            cat(text2[i],"\n",sep="")
        }
    }
    # text is the test assumption
    text <- paste("Test Information: This test determines whether there is an association between categorical variables ",var1," and ",var2,".
                   \n The test assumes that observations are independent, and all expected counts exceed 10 ((row_total x column_total)/table_total >= 10).
                    \r ****************************************************************
                   \n \n",sep="")
    wrapper(text)

    pval <- x$p.value
    alpha <- .05 

    # text1 is the test results
    if(pval >= alpha){
        text1 <- paste("Test Results: There is no statistical association between ",var1," and ",var2,". \n \n",sep="")
        wrapper(text1)
    }
    else if(pval < alpha){
        text1 <- paste("Test Results: There is a statistically significant association between ",var1," and ",var2,". (chi-square = ",round(x$statistic,3)," p=",round(pval,3),").",sep="")
        wrapper(text1)
    }
}

# modified twoWayTable function from Rcmdr: R Commander
twoWayTable2 <- function(){ # dialog memory 2011-06-27 J. Fox
	Library("abind")
	defaults <- list(initial.row=NULL, initial.column=NULL, 
			initial.percents="none", initial.chisq=1, initial.chisqComp=0, initial.expected=0, 
			initial.fisher=0, initial.subset=gettextRcmdr("<all valid cases>"))
	dialog.values <- getDialog("twoWayTable2", defaults)
	initializeDialog(title=gettextRcmdr("Two-Way Table"))
	variablesFrame <- tkframe(top)
	.factors <- Factors()
	rowBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Row variable (pick one)"),
			initialSelection=varPosn(dialog.values$initial.row, "factor"))
	columnBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Column variable (pick one)"),
			initialSelection=varPosn(dialog.values$initial.column, "factor"))
	subsetBox(subset.expression=dialog.values$initial.subset)
	onOK <- function(){
		row <- getSelection(rowBox)
		column <- getSelection(columnBox)
		percents <- as.character(tclvalue(percentsVariable))
		chisq <- tclvalue(chisqTestVariable)
		chisqComp <- tclvalue(chisqComponentsVariable)
		expected <- tclvalue(expFreqVariable)
		fisher <- tclvalue(fisherTestVariable)
		initial.subset <- subset <- tclvalue(subsetVariable)
		subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) ""
				else paste(", subset=", subset, sep="")
		putDialog("twoWayTable2", list(
						initial.row=row, 
						initial.column=column, 
						initial.percents=percents, initial.chisq=chisq, initial.chisqComp=chisqComp, 
						initial.expected=expected, initial.fisher=fisher, initial.subset=initial.subset
				))
		if (length(row) == 0 || length(column) == 0){
			errorCondition(recall=twoWayTable2, message=gettextRcmdr("You must select two variables."))
			return()
		}
		if (row == column) {
			errorCondition(recall=twoWayTable2, message=gettextRcmdr("Row and column variables are the same."))
			return()
		}
		closeDialog()
		command <- paste("xtabs(~", row, "+", column, ", data=", ActiveDataSet(),
				subset, ")", sep="")
		logger(paste(".Table <- ", command, sep=""))
		doItAndPrint(paste(".Table","<-",command))
		doItAndPrint(".Table")
		if (percents == "row") doItAndPrint("rowPercents(.Table) # Row Percentages")
		if (percents == "column") doItAndPrint("colPercents(.Table) # Column Percentages")
		if (percents == "total") doItAndPrint("totPercents(.Table) # Percentage of Total")
		if (chisq == 1) {
			command <- "chisq.test(.Table, correct=FALSE)"
			logger(paste(".Test <- ", command, sep=""))
			doItAndPrint(paste(".Test","<-",command))
			doItAndPrint(".Test")
            # Inserted Code:
            doItAndPrint(paste("chiSquareWords(.Test,",'"',row,'"',",",'"',column,'"',")",sep=""))
            # End Insertion
			if (expected == 1) doItAndPrint(".Test$expected # Expected Counts")
			warnText <- NULL
			if (0 < (nlt1 <- sum(.Test$expected < 1))) warnText <- paste(nlt1,
						gettextRcmdr("expected frequencies are less than 1"))
			if (0 < (nlt5 <- sum(.Test$expected < 5))) warnText <- paste(warnText, "\n", nlt5,
						gettextRcmdr(" expected frequencies are less than 5"), sep="")
			if (!is.null(warnText)) Message(message=warnText,
						type="warning")
			if (chisqComp == 1) {
				command <- "round(.Test$residuals^2, 2) # Chi-square Components"
				doItAndPrint(command)
			}
			logger("remove(.Test)")
			# remove(.Test)   ######### used to specify envir=.GlobalEnv
		}
		if (fisher == 1) doItAndPrint("fisher.test(.Table)")
                # Inserted Code:
		doItAndPrint(paste("chiSquareWords(fisher.test(.Table),",'"',row,'"',",",'"',column,'"',")",sep=""))
                # End Insertion
		logger("remove(.Table)")
		# remove(.Table) # used to specify envir=.GlobalEnv) 
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="xtabs", reset="twoWayTable2")
	radioButtons(name="percents",
			buttons=c("rowPercents", "columnPercents", "totalPercents", "nonePercents"),
			values=c("row", "column", "total", "none"), initialValue=dialog.values$initial.percents,
			labels=gettextRcmdr(c("Row percentages", "Column percentages", "Percentages of total", "No percentages")), title=gettextRcmdr("Compute Percentages"))
	checkBoxes(frame="testsFrame", boxes=c("chisqTest", "chisqComponents", "expFreq", "fisherTest"), 
			initialValues=c(dialog.values$initial.chisq, dialog.values$initial.chisqComp, 
					dialog.values$initial.expected, dialog.values$initial.fisher),
			labels=gettextRcmdr(c("Chi-square test of independence", "Components of chi-square statistic",
							"Print expected frequencies", "Fisher's exact test")))
	tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text="    "), getFrame(columnBox), sticky="nw")
	tkgrid(variablesFrame, sticky="w")
	tkgrid(percentsFrame, sticky="w")
	tkgrid(labelRcmdr(top, text=gettextRcmdr("Hypothesis Tests"), fg="blue"), sticky="w")
	tkgrid(testsFrame, sticky="w")
	tkgrid(subsetFrame, sticky="w")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=6, columns=1)
}
