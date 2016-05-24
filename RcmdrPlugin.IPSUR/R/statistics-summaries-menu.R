# Last modified Feb 16, 2008


`frequencyDistribution.ipsur` <-
function () 
{
    initializeDialog(title = gettextRcmdr("Frequency Distribution"))
    xBox <- variableListBox(top, Variables(), title = gettextRcmdr("Variable (pick one)"))
    optionsFrame <- tkframe(top)
    goodnessOfFitVariable <- tclVar("0")
    goodnessOfFitCheckBox <- tkcheckbutton(optionsFrame, variable = goodnessOfFitVariable)
    onOK <- function() {
        x <- getSelection(xBox)
        if (length(x) == 0) {
            errorCondition(recall = frequencyDistribution.ipsur, message = gettextRcmdr("You must select a variable."))
            return()
        }
        goodnessOfFit <- tclvalue(goodnessOfFitVariable)
        closeDialog()
        .activeDataSet <- ActiveDataSet()
        command <- paste(".Table <- table(", .activeDataSet, "$", x, ")", 
            sep = "")
        logger(paste(".Table <-", command))
        justDoIt(command)
        doItAndPrint(".Table  # frequencies")
        doItAndPrint(".Table/sum(.Table)  # relative frequencies")
        env <- environment()
        if (goodnessOfFit == 1) {
            initializeDialog(subwin, title = gettextRcmdr("Goodness-of-Fit Test"))
            hypothesisFrame <- tkframe(subwin)
            levs <- eval(parse(text = paste("levels(as.factor(", 
                .activeDataSet, "$", x, "))", sep = "")))
            n.levs <- length(levs)
            assign(".entry.1", tclVar(paste("1/", n.levs, sep = "")), 
                envir = env)
            make.entries <- "tklabel(hypothesisFrame, text='Hypothesized probabilities:   ')"
            make.lev.names <- "tklabel(hypothesisFrame, text='Factor levels:')"
            for (i in 1:n.levs) {
                entry.varname <- paste(".entry.", i, sep = "")
                assign(entry.varname, tclVar(paste("1/", n.levs, 
                  sep = "")), envir = env)
                make.entries <- paste(make.entries, ", ", "tkentry(hypothesisFrame, width='5', textvariable=", 
                  entry.varname, ")", sep = "")
                make.lev.names <- paste(make.lev.names, ", tklabel(hypothesisFrame, text='", 
                  levs[i], "')", sep = "")
            }
            eval(parse(text = paste("tkgrid(", make.lev.names, 
                ", sticky='w')", sep = "")), envir = env)
            eval(parse(text = paste("tkgrid(", make.entries, 
                ", stick='w')", sep = "")), envir = env)
            tkgrid(hypothesisFrame, sticky = "w")
            onOKsub <- function() {
                probs <- rep(NA, n.levs)
                for (i in 1:n.levs) {
                  entry.varname <- paste(".entry.", i, sep = "")
                  res <- try(entry <- eval(parse(text = eval(parse(text = paste("tclvalue(", 
                    entry.varname, ")", sep = "")), envir = env))), 
                    silent = TRUE)
                  if (class(res) == "try-error") {
                    errorCondition(subwin, message = gettextRcmdr("Invalid entry."))
                    return()
                  }
                  if (length(entry) == 0) {
                    errorCondition(subwin, message = gettextRcmdr("Missing entry."))
                    return()
                  }
                  opts <- options(warn = -1)
                  probs[i] <- as.numeric(entry)
                  options(opts)
                }
                probs <- na.omit(probs)
                if (length(probs) != n.levs) {
                  errorCondition(subwin, message = sprintf(gettextRcmdr("Number of valid entries (%d)\nnot equal to number levels (%d)."), 
                    length(probs), n.levs))
                  return()
                }
                if (any(probs < 0)) {
                  errorCondition(subwin, message = gettextRcmdr("Negative probabilities not allowed."))
                  return()
                }
                if (abs(sum(probs) - 1) > 0.001) {
                  Message(message = gettextRcmdr("Probabilities rescaled to sum to 1."), 
                    type = "warning")
                  probs <- probs/sum(probs)
                }
                closeDialog(subwin)
                command <- paste(".Probs <- c(", paste(probs, collapse = ","), 
                  ")", sep = "")
                logger(paste(".Probs <-", command))
                justDoIt(command)
                doItAndPrint("chisq.test(.Table, p=.Probs)")
                logger("remove(.Probs)")
                remove(.Probs, envir = .GlobalEnv)
            }
            subOKCancelHelp(subwin)
            tkgrid(subButtonsFrame, sticky = "w")
            dialogSuffix(subwin, rows = 2, columns = 1, onOK = onOKsub, 
                focus = subwin)
        }
        logger("remove(.Table)")
        remove(.Table, envir = .GlobalEnv)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "table")
    tkgrid(getFrame(xBox), sticky = "nw")
    tkgrid(tklabel(optionsFrame, text = gettextRcmdr("Chi-square goodness-of-fit test")), 
        goodnessOfFitCheckBox, sticky = "w")
    tkgrid(optionsFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix(rows = 3, columns = 2)
}



`numericalSummaries.ipsur` <-
function () 
{
    require(abind)
    require(e1071)
    initializeDialog(title = gettextRcmdr("Numerical Summaries"))
    xBox <- variableListBox(top, Numeric(), selectmode = "multiple", 
        title = gettextRcmdr("Variables (pick one or more)"))
    checkBoxes(frame = "checkBoxFrame", boxes = c("mean", "sd", 
        "skewness", "kurtosis"), initialValues = c("1", "1", 
        "1", "1"), labels = gettextRcmdr(c("Mean", "Standard Deviation", 
        "Skewness", "Kurtosis")))
    quantilesVariable <- tclVar("0")
    quantilesFrame <- tkframe(top)
    quantilesCheckBox <- tkcheckbutton(quantilesFrame, variable = quantilesVariable)
    quantiles <- tclVar("0,.25,.5,.75,1")
    quantilesEntry <- tkentry(quantilesFrame, width = "20", textvariable = quantiles)
    groupsBox(recall = numericalSummaries.ipsur, label = gettextRcmdr("Summarize by:"), 
        initialLabel = gettextRcmdr("Summarize by groups"))
    onOK <- function() {
        x <- getSelection(xBox)
        if (length(x) == 0) {
            errorCondition(recall = numericalSummaries.ipsur, message = gettextRcmdr("You must select a variable."))
            return()
        }
        closeDialog()
        quants <- paste("c(", gsub(" ", ",", tclvalue(quantiles)), 
            ")")
        .activeDataSet <- ActiveDataSet()
        vars <- if (length(x) == 1) 
            paste("\"", x, "\"", sep = "")
        else paste("c(", paste("\"", x, "\"", collapse = ", ", 
            sep = ""), ")", sep = "")
        vars <- paste(.activeDataSet, "[,", vars, "]", sep = "")
        stats <- paste("c(", paste(c("\"mean\"", "\"sd\"", "\"skewness\"", 
            "\"kurtosis\"", "\"quantiles\"")[c(tclvalue(meanVariable), 
            tclvalue(sdVariable), tclvalue(skewnessVariable), 
            tclvalue(kurtosisVariable), tclvalue(quantilesVariable)) == 
            1], collapse = ", "), ")", sep = "")
        if (stats == "c()") {
            errorCondition(recall = numericalSummaries.ipsur, message = gettextRcmdr("No statistics selected."))
            return()
        }
        command <- if (.groups != FALSE) {
            grps <- paste(.activeDataSet, "$", .groups, sep = "")
            paste("numSummaryIPSUR(", vars, ", groups=", grps, 
                ", statistics=", stats, ")", sep = "")
        }
        else paste("numSummaryIPSUR(", vars, ", statistics=", 
            stats, ")", sep = "")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "numSummaryIPSUR")
    tkgrid(getFrame(xBox), sticky = "nw")
    tkgrid(checkBoxFrame, sticky = "w")
    tkgrid(tklabel(quantilesFrame, text = gettextRcmdr("Quantiles")), 
        quantilesCheckBox, tklabel(quantilesFrame, text = gettextRcmdr(" quantiles:")), 
        quantilesEntry, sticky = "w")
    tkgrid(quantilesFrame, sticky = "w")
    tkgrid(groupsFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky = "w")
    dialogSuffix(rows = 6, columns = 1)
}
