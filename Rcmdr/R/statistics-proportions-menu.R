# Statistics Menu dialogs

# last modified 2015-12-15 by J. Fox

    # Proportions menu

singleProportionTest <- function () {
    defaults <- list (initial.x = NULL, initial.alternative = "two.sided", initial.level = ".95", 
        initial.test = "normal" , initial.p = ".5", initial.tab=0)
    dialog.values <- getDialog ("singleProportionTest", defaults)
    initializeDialog(title = gettextRcmdr("Single-Sample Proportion Test"), use.tabs=TRUE)
    xBox <- variableListBox(dataTab, TwoLevelFactors(), title = gettextRcmdr("Variable (pick one)"),
        initialSelection = varPosn(dialog.values$initial.x,"twoLevelFactor"))
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        x <- getSelection(xBox)
        if (length(x) == 0) {
            errorCondition(recall = singleProportionTest, message = gettextRcmdr("You must select a variable."))
            return()
        }
        alternative <- as.character(tclvalue(alternativeVariable))
        level <- tclvalue(confidenceLevel)
        test <- as.character(tclvalue(testVariable))
        p <- tclvalue(pVariable)
        putDialog ("singleProportionTest", list (initial.x = x, initial.alternative = alternative, 
            initial.level = level, initial.test = test ,initial.p = p, initial.tab=tab))
        closeDialog()
        command <- paste("local({\n  .Table <- xtabs(~", x, ", data=", ActiveDataSet(), 
            ')\n  cat("\\nFrequency counts (test is for first level):\\n")\n  print(.Table)')
        if (test == "normal") 
            doItAndPrint(paste(command, "\n  prop.test(rbind(.Table), alternative='", 
                alternative, "', p=", p, ", conf.level=", level, 
                ", correct=FALSE)\n})", sep = ""))
        else if (test == "corrected") 
            doItAndPrint(paste(command, "\n  prop.test(rbind(.Table), alternative='", 
                alternative, "', p=", p, ", conf.level=", level, 
                ", correct=TRUE)\n})", sep = ""))
        else doItAndPrint(paste(command, "\n  binom.test(rbind(.Table), alternative='", 
            alternative, "', p=", p, ", conf.level=", level, 
            ")\n})", sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "prop.test", reset = "singleProportionTest", apply = "singleProportionTest")
    radioButtons(optionsTab, name = "alternative", buttons = c("twosided", 
        "less", "greater"), values = c("two.sided", "less", "greater"), 
        labels = gettextRcmdr(c("Population proportion != p0", 
            "Population proportion < p0", "Population proportion > p0")), 
        title = gettextRcmdr("Alternative Hypothesis"), initialValue = dialog.values$initial.alternative)
    rightFrame <- tkframe(optionsTab)
    confidenceFrame <- tkframe(rightFrame)
    confidenceLevel <- tclVar(dialog.values$initial.level)
    confidenceField <- ttkentry(confidenceFrame, width = "6", 
        textvariable = confidenceLevel)
    pFrame <- tkframe(rightFrame)
    pVariable <- tclVar(dialog.values$initial.p)
    pField <- ttkentry(pFrame, width = "6", textvariable = pVariable)
    radioButtons(optionsTab, name = "test", buttons = c("normal", "corrected", 
        "exact"), labels = gettextRcmdr(c("Normal approximation", 
            "Normal approximation with\ncontinuity correction", "Exact binomial")), 
        title = gettextRcmdr("Type of Test"), initialValue = dialog.values$initial.test)
    tkgrid(getFrame(xBox), sticky = "nw")
    tkgrid(labelRcmdr(pFrame, text = gettextRcmdr("Null hypothesis: p = "), 
        fg = getRcmdr("title.color"), font="RcmdrTitleFont"), pField, sticky = "w")
    tkgrid(pFrame, sticky = "w")
    tkgrid(labelRcmdr(rightFrame, text = ""))
    tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level: "), 
        fg = getRcmdr("title.color"), font="RcmdrTitleFont"), confidenceField, sticky = "w")
    tkgrid(confidenceFrame, sticky = "w")
    tkgrid(alternativeFrame, labelRcmdr(optionsTab, text="   "), rightFrame, sticky = "nw")
    tkgrid(testFrame, sticky = "w")
    tkgrid.configure(confidenceField, sticky = "e")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

twoSampleProportionsTest <- function () {
    Library("abind")
    defaults <- list(initial.groups = NULL, initial.response = NULL, initial.alternative = "two.sided", 
        initial.confidenceLevel = ".95", initial.test = "normal", initial.label=NULL, initial.tab=0)
    dialog.values <- getDialog("twoSampleProportionsTest", defaults)
    initializeDialog(title = gettextRcmdr("Two-Sample Proportions Test"), use.tabs=TRUE)
    .twoLevelFactors <- TwoLevelFactors()
    groupsBox <- variableListBox(dataTab, .twoLevelFactors, title = gettextRcmdr("Groups (pick one)"), 
        initialSelection = varPosn(dialog.values$initial.groups, "twoLevelFactor"))
    xBox <- variableListBox(dataTab, .twoLevelFactors, title = gettextRcmdr("Response Variable (pick one)"), 
        initialSelection = varPosn(dialog.values$initial.response, "twoLevelFactor"))
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        groups <- getSelection(groupsBox)
        if (length(groups) == 0) {
            errorCondition(recall = twoSampleProportionsTest, 
                message = gettextRcmdr("You must select a groups variable."))
            return()
        }
        x <- getSelection(xBox)
        if (length(x) == 0) {
            errorCondition(recall = twoSampleProportionsTest, 
                message = gettextRcmdr("You must select a response variable."))
            return()
        }
        if (x == groups) {
            errorCondition(recall = twoSampleProportionsTest, 
                message = gettextRcmdr("Groups and response variables must be different."))
            return()
        }
        alternative <- as.character(tclvalue(alternativeVariable))
        level <- tclvalue(confidenceLevel)
        test <- as.character(tclvalue(testVariable))
        closeDialog()
        putDialog("twoSampleProportionsTest", list(initial.groups = groups, initial.response = x, 
            initial.test = test, initial.alternative = alternative, initial.confidenceLevel = level,
            initial.label=.groupsLabel, initial.tab=tab))
        command <- paste("local({  .Table <- xtabs(~", groups, "+", x, ", data=", 
            ActiveDataSet(), ")", sep = "")
        command <- paste(command, '\n  cat("\\nPercentage table:\\n")', sep="")
        command <- paste(command, "\n  print(rowPercents(.Table))", sep="")
        if (test == "normal") 
            doItAndPrint(paste(command, "\n  prop.test(.Table, alternative='", 
                alternative, "', conf.level=", level, ", correct=FALSE)\n})", 
                sep = ""))
        else doItAndPrint(paste(command, "\n  prop.test(.Table, alternative='", 
            alternative, "', conf.level=", level, ", correct=TRUE)\n})", 
            sep = ""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "prop.test", reset = "twoSampleProportionsTest", apply = "twoSampleProportionsTest")
    radioButtons(optionsTab, name = "alternative", buttons = c("twosided", 
        "less", "greater"), values = c("two.sided", "less", "greater"), 
        labels = gettextRcmdr(c("Two-sided", "Difference < 0", 
            "Difference > 0")), initialValue = dialog.values$initial.alternative, 
        title = gettextRcmdr("Alternative Hypothesis"))
    rightFrame <- tkframe(optionsTab)
    confidenceFrame <- tkframe(rightFrame)
    confidenceLevel <- tclVar(dialog.values$initial.confidenceLevel)
    confidenceField <- ttkentry(confidenceFrame, width = "6", 
        textvariable = confidenceLevel)
    radioButtons(optionsTab, name = "test", buttons = c("normal", "corrected"), 
        labels = gettextRcmdr(c("Normal approximation", "Normal approximation with\ncontinuity correction")), 
        initialValue = dialog.values$initial.test, 
        title = gettextRcmdr("Type of Test"))
    tkgrid(getFrame(groupsBox), labelRcmdr(dataTab, text="  "), getFrame(xBox), sticky = "nw")
    groupsLabel(optionsTab, initialText=dialog.values$initial.label)
    tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Confidence Level: "), 
        fg = getRcmdr("title.color"), font="RcmdrTitleFont"), confidenceField, sticky = "w")
    tkgrid(confidenceFrame, sticky = "w")
    tkgrid(alternativeFrame, rightFrame, sticky = "nw")
    tkgrid(testFrame, sticky = "w")
    tkgrid.configure(confidenceField, sticky = "e")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

