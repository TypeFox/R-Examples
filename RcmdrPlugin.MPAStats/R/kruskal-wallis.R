#Last modified on November 26, 2013 by Jordan Gressel

# Interpretation Function
kruskalWallisWords <- function(x,group,response){
    wrapper <- function(text){
        text2 <- strwrap(text)
        for(i in 1:length(text2)){
            cat(text2[i],"\n",sep="")
        }
    }


    pval <- round(x$p.value,3)
    alpha <- .05
    statistic <- round(x$statistic,3)
    
    # text is the test assumption
    text <- paste("Test Information: This test determines whether any level of ",group," affects ",response,"; that is, whether the true median ",response," levels of ",group," are all equivalent to each other, or if the median ",response," of at least one level of ",group," is significantly different from the others.
                   \n The test assumes that observations are independent and all samples come from populations with the same-shaped distribution.
                    \r ****************************************************************
                   \n \n",sep="")
    wrapper(text)
    
    # text1 is the test results
    if(pval >= alpha){
        text1 <- paste("Test Results: There is no significant difference in the median ",response," between levels of ",group,". (chi-square=",statistic,", p=",pval,").",sep="")
        wrapper(text1)
    }
    else if(pval < alpha){
        text1 <- paste("Test Results: At least one median ", response," among the levels of ", group, " differs from the rest. (chi-square=",statistic,", p=",pval,").",sep="")
        wrapper(text1)
    }
}

# Modified KruskalWallisTest from Rcmdr: R Commander
KruskalWallisTest2 <- function () {
	defaults <- list(initial.group = NULL, initial.response = NULL)
	dialog.values <- getDialog("KruskalWallisTest2", defaults)
	initializeDialog(title = gettextRcmdr("Kruskal-Wallis Rank Sum Test"))
	groupBox <- variableListBox(top, Factors(), title = gettextRcmdr("Groups (pick one)"),
			initialSelection = varPosn(dialog.values$initial.group, "factor"))
	responseBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Response Variable (pick one)"),
			initialSelection = varPosn(dialog.values$initial.response, "numeric"))
	onOK <- function() {
		group <- getSelection(groupBox)
		if (length(group) == 0) {
			errorCondition(recall = KruskalWallisTest2, message = gettextRcmdr("You must select a groups variable."))
			return()
		}
		response <- getSelection(responseBox)
		closeDialog()
		putDialog("KruskalWallisTest2", list(initial.group = group, initial.response = response))
		if (length(response) == 0) {
			errorCondition(recall = KruskalWallisTest2, message = gettextRcmdr("You must select a response variable."))
			return()
		}
		.activeDataSet <- ActiveDataSet()
		doItAndPrint(paste("tapply(", paste(.activeDataSet, "$", 
								response, sep = ""), ", ", paste(.activeDataSet, 
								"$", group, sep = ""), ", median, na.rm=TRUE)", sep = ""))
                # Edited Code (added "kruskal <- ")
		doItAndPrint(paste("kruskal <- kruskal.test(", response, " ~ ", 
						group, ", data=", .activeDataSet, ")", sep = ""))
        # Inserted Code:
        doItAndPrint("kruskal")
        doItAndPrint(paste("kruskalWallisWords(kruskal,",'"',group,'"',",",'"',response,'"',")",sep=""))
        # End Insertion
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "kruskal.test", reset = "KruskalWallisTest2")
	tkgrid(getFrame(groupBox), getFrame(responseBox), sticky = "nw")
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	dialogSuffix(rows = 2, columns = 2)
}
