# Modified on November 27, 2013 by Jordan Gressel

# Interpretation function
wilcoxonMannWhitneyWords <- function(x,group,response){
    wrapper <- function(text){
        text2 <- strwrap(text)
        for(i in 1:length(text2)){
            cat(text2[i],"\n",sep="")
        }
    }

    pval <- x$p.value
    alpha <- .05
    statistic <- x$statistic
    
    # text is the test assumption
    text <- paste("Test Information: This test determines whether level of ",group," affects level of ",response,"; that is, whether the population median of ",response," differs between the two levels of ",group,".
                    \n The test assumes that all responses are ordinal and sampled independently.
                  \r ****************************************************************
                   \n \n",sep="")
    wrapper(text)

    if(pval >= alpha){
        text1 <- paste("Test Results: There is no significant difference in the median ",response," between the two levels of ",group,". (W=",round(statistic,3),", p=",round(pval,3),").",sep="")
        wrapper(text1)
    }
    else if(pval < alpha){
        text1 <- paste("Test Results: There is a significant difference in the median ", response," between the two levels of ",group,". (W=",round(statistic,3),", p=",round(pval,3),").",sep="")
        wrapper(text1)
    }
}

# Modified from twoSampleWilcoxonTest from Rcmdr: R Commander
twoSampleWilcoxonTest2 <- function () {
	defaults <- list(initial.group = NULL, initial.response = NULL, initial.alternative = "two.sided", 
			initial.test = "default", initial.label=NULL)
	dialog.values <- getDialog("twoSampleWilcoxonTest2", defaults)
	initializeDialog(title = gettextRcmdr("Two-Sample Wilcoxon Test"))
	groupBox <- variableListBox(top, TwoLevelFactors(), title = gettextRcmdr("Groups (pick one)"),
			initialSelection = varPosn(dialog.values$initial.group, "twoLevelFactor"))
	responseBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Response Variable (pick one)"),
			initialSelection = varPosn(dialog.values$initial.response, "numeric"))
	onOK <- function() {
		group <- getSelection(groupBox)
		if (length(group) == 0) {
			errorCondition(recall = twoSampleWilcoxonTest2, message = gettextRcmdr("You must select a groups variable."))
			return()
		}
		response <- getSelection(responseBox)
		if (length(response) == 0) {
			errorCondition(recall = twoSampleWilcoxonTest2, message = gettextRcmdr("You must select a response variable."))
			return()
		}
		alternative <- as.character(tclvalue(alternativeVariable))
		test <- as.character(tclvalue(testVariable))
		closeDialog()
		putDialog("twoSampleWilcoxonTest2", list(initial.group = group, initial.response = response, 
						initial.test = test, initial.alternative = alternative, initial.label=.groupsLabel))
		.activeDataSet <- ActiveDataSet()
		doItAndPrint(paste("tapply(", paste(.activeDataSet, "$", 
								response, sep = ""), ", ", paste(.activeDataSet, 
								"$", group, sep = ""), ", median, na.rm=TRUE)", sep = ""))
                # adding "wilcox <-"
		if (test == "default") {
			doItAndPrint(paste("wilcox <- wilcox.test(", response, " ~ ", 
							group, ", alternative=\"", alternative, "\", data=", 
							.activeDataSet, ")", sep = ""))
            doItAndPrint("wilcox")
		}
		else{doItAndPrint(paste("wilcox <- wilcox.test(", response, " ~ ", 
							group, ", alternative='", alternative, "', exact=", 
							test == "exact", ", correct=", test == "correct", 
							", data=", .activeDataSet, ")", sep = ""))
            doItAndPrint("wilcox")
        }
        # Inserted Code:
        doItAndPrint(paste("wilcoxonMannWhitneyWords(wilcox,",'"',group,'"',",",'"',response,'"',")",sep=""))
        # End Insertion
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "wilcox.test", reset = "twoSampleWilcoxonTest2")
	radioButtons(name = "alternative", buttons = c("twosided", 
					"less", "greater"), values = c("two.sided", "less", "greater"), 
			labels = gettextRcmdr(c("Two-sided", "Difference < 0", 
							"Difference > 0")), initialValue = dialog.values$initial.alternative,
			title = gettextRcmdr("Alternative Hypothesis"))
	radioButtons(name = "test", buttons = c("default", "exact", 
					"normal", "correct"), labels = gettextRcmdr(c("Default", 
							"Exact", "Normal approximation", "Normal approximation with\ncontinuity correction")), 
			initialValue = dialog.values$initial.test,
			title = gettextRcmdr("Type of Test"))
	tkgrid(getFrame(groupBox), getFrame(responseBox), sticky = "nw")
	groupsLabel(groupsBox = groupBox, columnspan = 2, initialText=dialog.values$initial.label)
	tkgrid(alternativeFrame, testFrame, sticky = "nw")
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	dialogSuffix(rows = 4, columns = 2)
}

