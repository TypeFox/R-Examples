# Last modified Feb 16, 2008


`powerAnovatest` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Power for balanced ANOVA"))
    optionsFrame <- tkframe(top)
    groups <- tclVar(gettextRcmdr(""))
    entryGroups <- tkentry(optionsFrame, width = "10", textvariable = groups)
    ssize <- tclVar(gettextRcmdr(""))
    entrySsize <- tkentry(optionsFrame, width = "10", textvariable = ssize)
    between.var <- tclVar(gettextRcmdr(""))
    entryBetween.var <- tkentry(optionsFrame, width = "10", textvariable = between.var)
    within.var <- tclVar(gettextRcmdr(""))
    entryWithin.var <- tkentry(optionsFrame, width = "10", textvariable = within.var)
    sig.level <- tclVar(gettextRcmdr("0.05"))
    entrySig.level <- tkentry(optionsFrame, width = "10", textvariable = sig.level)
    power <- tclVar(gettextRcmdr(""))
    entryPower <- tkentry(optionsFrame, width = "10", textvariable = power)
    onOK <- function() {
        closeDialog()
        groupsValue <- trim.blanks(tclvalue(groups))
        ssizeValue <- trim.blanks(tclvalue(ssize))
        between.varValue <- trim.blanks(tclvalue(between.var))
        within.varValue <- trim.blanks(tclvalue(within.var))
        sig.levelValue <- trim.blanks(tclvalue(sig.level))
        powerValue <- trim.blanks(tclvalue(power))
        if (ssizeValue == "") {
            ssizeValue = "NULL"
        }
        else if (groupsValue == "") {
            groupsValue = "NULL"
        }
        else if (between.varValue == "") {
            between.varValue = "NULL"
        }
        else if (within.varValue == "") {
            within.varValue = "NULL"
        }
        else if (sig.levelValue == "") {
            sig.levelValue = "NULL"
        }
        else if (powerValue == "") {
            powerValue = "NULL"
        }
        else {
            errorCondition(recall = powerAnovatest, message = gettextRcmdr("Exactly one field must be left empty."))
            return()
        }
        command <- paste("power.anova.test(groups=", groupsValue, 
            ", n=", ssizeValue, ", between.var=", between.varValue, 
            ", within.var=", within.varValue, ", sig.level=", 
            sig.levelValue, ", power=", powerValue, ")", sep = "")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "power.anova.test")
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters: (leave exactly one field blank)"), 
        fg = "blue"), columnspan = 2, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Number of groups:")), 
        entryGroups, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Number of observations (per group):")), 
        entrySsize, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Between group variance:")), 
        entryBetween.var, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Within group variance:")), 
        entryWithin.var, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Significance level (Type I error probability):")), 
        entrySig.level, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Power of test (1 minus Type II error probability):")), 
        entryPower, sticky = "w")
    tkgrid(optionsFrame, sticky = "w")
    dialogSuffix(rows = 1, columns = 1)
}


`powerProptest` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Power for two Proportions"))
    optionsFrame <- tkframe(top)
    ssize <- tclVar(gettextRcmdr(""))
    entrySsize <- tkentry(optionsFrame, width = "10", textvariable = ssize)
    p1 <- tclVar(gettextRcmdr(""))
    entryP1 <- tkentry(optionsFrame, width = "10", textvariable = p1)
    p2 <- tclVar(gettextRcmdr(""))
    entryP2 <- tkentry(optionsFrame, width = "10", textvariable = p2)
    sig.level <- tclVar(gettextRcmdr("0.05"))
    entrySig.level <- tkentry(optionsFrame, width = "10", textvariable = sig.level)
    power <- tclVar(gettextRcmdr(""))
    entryPower <- tkentry(optionsFrame, width = "10", textvariable = power)
    strictVariable <- tclVar("1")
    strictCheckBox <- tkcheckbutton(optionsFrame, variable = strictVariable)
    radioButtons(optionsFrame, "alternative", buttons = c("two.sided", 
        "one.sided"), labels = gettextRcmdr(c("Two sided", "One sided")), 
        title = gettextRcmdr("Alternative Hypothesis"))
    onOK <- function() {
        closeDialog()
        ssizeValue <- trim.blanks(tclvalue(ssize))
        p1Value <- trim.blanks(tclvalue(p1))
        p2Value <- trim.blanks(tclvalue(p2))
        sig.levelValue <- trim.blanks(tclvalue(sig.level))
        powerValue <- trim.blanks(tclvalue(power))
        strictValue <- tclvalue(strictVariable) == "1"
        alternativeValue <- tclvalue(alternativeVariable)
        if (ssizeValue == "") {
            ssizeValue = "NULL"
        }
        else if (p1Value == "") {
            p1Value = "NULL"
        }
        else if (p2Value == "") {
            p2Value = "NULL"
        }
        else if (sig.levelValue == "") {
            sig.levelValue = "NULL"
        }
        else if (powerValue == "") {
            powerValue = "NULL"
        }
        else {
            errorCondition(recall = powerProptest, message = gettextRcmdr("Exactly one field must be left empty."))
            return()
        }
        command <- paste("power.prop.test(n=", ssizeValue, ", p1=", 
            p1Value, ", p2=", p2Value, ", sig.level=", sig.levelValue, 
            ", power=", powerValue, ", alternative=\"", alternativeValue, 
            "\", strict=", strictValue, ")", sep = "")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "power.prop.test")
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters: (leave exactly one field blank)"), 
        fg = "blue"), columnspan = 2, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Number of observations (per group):")), 
        entrySsize, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("probability in one group:")), 
        entryP1, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("probability in other group:")), 
        entryP2, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Significance level (Type I error probability):")), 
        entrySig.level, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Power of test (1 minus Type II error probability):")), 
        entryPower, sticky = "w")
    tkgrid(optionsFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    tkgrid(alternativeFrame, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Use strict interpretation in two-sided case:")), 
        strictCheckBox, sticky = "w")
    dialogSuffix(rows = 3, columns = 1)
}



`powerTtest` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Power Calculations for t-Tests"))
    optionsFrame <- tkframe(top)
    ssize <- tclVar(gettextRcmdr(""))
    entrySsize <- tkentry(optionsFrame, width = "10", textvariable = ssize)
    delta <- tclVar(gettextRcmdr(""))
    entryDelta <- tkentry(optionsFrame, width = "10", textvariable = delta)
    sd <- tclVar(gettextRcmdr("1"))
    entrySd <- tkentry(optionsFrame, width = "10", textvariable = sd)
    sig.level <- tclVar(gettextRcmdr("0.05"))
    entrySig.level <- tkentry(optionsFrame, width = "10", textvariable = sig.level)
    power <- tclVar(gettextRcmdr(""))
    entryPower <- tkentry(optionsFrame, width = "10", textvariable = power)
    strictVariable <- tclVar("1")
    strictCheckBox <- tkcheckbutton(optionsFrame, variable = strictVariable)
    radioButtons(optionsFrame, "type", buttons = c("two.sample", 
        "one.sample", "paired"), labels = gettextRcmdr(c("Two sample", 
        "One sample", "Paired")), title = gettextRcmdr("Type of t Test"))
    radioButtons(optionsFrame, "alternative", buttons = c("two.sided", 
        "one.sided"), labels = gettextRcmdr(c("Two sided", "One sided")), 
        title = gettextRcmdr("Alternative Hypothesis"))
    onOK <- function() {
        closeDialog()
        ssizeValue <- trim.blanks(tclvalue(ssize))
        deltaValue <- trim.blanks(tclvalue(delta))
        sdValue <- trim.blanks(tclvalue(sd))
        sig.levelValue <- trim.blanks(tclvalue(sig.level))
        powerValue <- trim.blanks(tclvalue(power))
        strictValue <- tclvalue(strictVariable) == "1"
        typeValue <- tclvalue(typeVariable)
        alternativeValue <- tclvalue(alternativeVariable)
        if (ssizeValue == "") {
            ssizeValue = "NULL"
        }
        else if (deltaValue == "") {
            deltaValue = "NULL"
        }
        else if (sdValue == "") {
            sdValue = "NULL"
        }
        else if (sig.levelValue == "") {
            sig.levelValue = "NULL"
        }
        else if (powerValue == "") {
            powerValue = "NULL"
        }
        else {
            errorCondition(recall = powerTtest, message = gettextRcmdr("Exactly one field must be left empty."))
            return()
        }
        command <- paste("power.t.test(n=", ssizeValue, ", delta=", 
            deltaValue, ", sd=", sdValue, ", sig.level=", sig.levelValue, 
            ", power=", powerValue, ", type=\"", typeValue, "\", alternative=\"", 
            alternativeValue, "\", strict=", strictValue, ")", 
            sep = "")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "power.t.test")
    tkgrid(tklabel(top, text = gettextRcmdr("Parameters: (leave exactly one field blank)"), 
        fg = "blue"), columnspan = 2, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Number of observations (per group):")), 
        entrySsize, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("True difference in means:")), 
        entryDelta, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Standard deviation:")), 
        entrySd, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Significance level (Type I error probability):")), 
        entrySig.level, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Power of test (1 minus Type II error probability):")), 
        entryPower, sticky = "w")
    tkgrid(typeFrame, sticky = "w", columnspan = 2)
    tkgrid(optionsFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    tkgrid(alternativeFrame, sticky = "w")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Use strict interpretation in two-sided case:")), 
        strictCheckBox, sticky = "w")
    dialogSuffix(rows = 3, columns = 1)
}
