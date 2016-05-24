# Modified on November 27, 2013 by Jordan Gressel

# Interpretation function
wilcoxonWords <- function(x,var1,var2){
    wrapper <- function(text){
        text2 <- strwrap(text)
        for(i in 1:length(text2)){
            cat(text2[i],"\n",sep="")
        }
    }

    pval <- x$p.value
    alpha <- .05
    alternative <- paste(x$alternative," than ",sep="")
    if(alternative == "two.sided than "){
        alternative <- "different from "
    }
    null <- x$null.value
    Wstat <- x$statistic
    
    # text is the test assumption
    text <- paste("Test Information: This test determines the median distance between ",var1," and ",var2," in the population, where ",var1," and ",var2," are organized into connected pairs; that is, whether the true mean population difference ",var1," - ",var2," is ",alternative,null,".
                  \n The test assumes that the data are paired and from same population, that each pair is chosen randomly and independently, and that the data are measured at least on an ordinal scale. Differences are distributed symmetrically.
                   \r ****************************************************************
                  \n \n",sep="")
    wrapper(text)
    

    if(pval >= alpha){
        text1 <- paste("Test Results: There is insufficient evidence to conclude that the true median paired distance between ",var1," and ",var2," is ",alternative,null,". (W=",round(Wstat,3),", p=",round(pval,3),").",sep="")
        wrapper(text1)
    }

    else if(pval < alpha){
        text1 <- paste("Test Results: The true median paired distance between ",var1," and ",var2," is ",alternative,null,". (W=",round(Wstat,3),", p=",round(pval,3),").",sep="")
        wrapper(text1)
    }    
}

# Modified from pairedWilcoxonTest from Rcmdr: R Commander
pairedWilcoxonTest2 <- function () {
	defaults <- list(initial.x = NULL, initial.y = NULL, initial.alternative = "two.sided", 
			initial.test = "default")
	dialog.values <- getDialog("pairedWilcoxonTest2", defaults)
	initializeDialog(title = gettextRcmdr("Paired Wilcoxon Test"))
	.numeric <- Numeric()
	xBox <- variableListBox(top, .numeric, title = gettextRcmdr("First variable (pick one)"), 
			initialSelection = varPosn(dialog.values$initial.x, "numeric"))
	yBox <- variableListBox(top, .numeric, title = gettextRcmdr("Second variable (pick one)"), 
			initialSelection = varPosn(dialog.values$initial.y, "numeric"))
	onOK <- function() {
		x <- getSelection(xBox)
		y <- getSelection(yBox)
		closeDialog()
		alternative <- as.character(tclvalue(alternativeVariable))
		test <- as.character(tclvalue(testVariable))
		putDialog("pairedWilcoxonTest2", list(initial.x = x, initial.y = y, 
						initial.test = test, initial.alternative = alternative))
		if (length(x) == 0 | length(y) == 0) {
			errorCondition(recall = pairedWilcoxonTest2, message = gettextRcmdr("You must select two variables."))
			return()
		}
		if (x == y) {
			errorCondition(recall = pairedWilcoxonTest2, message = gettextRcmdr("The two variables must be different."))
			return()
		}
		.activeDataSet <- ActiveDataSet()
		doItAndPrint(paste("median(", .activeDataSet, "$", x, 
						" - ", .activeDataSet, "$", y, ", na.rm=TRUE) # median difference", 
						sep = ""))


                # Inserting: "wilcox <- "
		if (test == "default") {
            doItAndPrint(paste("wilcox <- wilcox.test(", .activeDataSet, 
							"$", x, ", ", .activeDataSet, "$", y, ", alternative='", 
							alternative, "', paired=TRUE)", sep = ""))
            doItAndPrint("wilcox")
		}
		else if (test == "exact") {
			doItAndPrint(paste("wilcox <- wilcox.test(", .activeDataSet, 
							"$", x, ", ", .activeDataSet, "$", y, ", alternative='", 
							alternative, "', exact=TRUE, paired=TRUE)", sep = ""))
            doItAndPrint("wilcox")
		}
		else {
			doItAndPrint(paste("wilcox <- wilcox.test(", .activeDataSet, 
							"$", x, ", ", .activeDataSet, "$", y, ", alternative='", 
							alternative, "', correct=", test == "correct", 
							", exact=FALSE, paired=TRUE)", sep = ""))
            doItAndPrint("wilcox")
		}
        # Insertion:
        doItAndPrint(paste("wilcoxonWords(wilcox,",'"',x,'"',",",'"',y,'"',")",sep=""))
        # End Insertion
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject = "wilcox.test", reset = "pairedWilcoxonTest2")
	radioButtons(name = "alternative", buttons = c("twosided", 
					"less", "greater"), values = c("two.sided", "less", "greater"), 
			labels = gettextRcmdr(c("Two-sided", "Difference < 0", 
							"Difference > 0")), title = gettextRcmdr("Alternative Hypothesis"), 
			initialValue = dialog.values$initial.alternative)
	radioButtons(name = "test", buttons = c("default", "exact", 
					"normal", "correct"), labels = gettextRcmdr(c("Default", 
							"Exact", "Normal approximation", "Normal approximation with\ncontinuity correction")), 
			title = gettextRcmdr("Type of Test"), initialValue = dialog.values$initial.test)
	tkgrid(getFrame(xBox), getFrame(yBox), sticky = "nw")
	tkgrid(alternativeFrame, testFrame, sticky = "nw")
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	dialogSuffix(rows = 3, columns = 2)
}
