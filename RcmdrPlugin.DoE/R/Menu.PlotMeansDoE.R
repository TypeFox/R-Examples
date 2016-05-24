Menu.PlotMeansDoE <- function ()
{
    initializeDialog(title = gettextRcmdr("Main Effects and Interaction Plots (general)"))
    tkgrid(tklabel(top, text=gettextRcmdr("CAUTION: This dialog plots means only.\nFor unbalanced designs, hidden factors may be responsible for mean differences.\nEffects plots from linear models may be more useful.")),columnspan=2)
    groupBox <- variableListBox(top, Factors(), title = gettextRcmdr("Factors (select one or two)"),
        selectmode = "multiple")
    responseBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Response Variable (select one)"))
    onOK <- function() {
        groups <- getSelection(groupBox)
        response <- getSelection(responseBox)
        closeDialog()
        if (0 == length(groups)) {
            errorCondition(recall = PlotMeansDoE.menu, message = gettextRcmdr("No factors selected."))
            return()
        }
        if (2 < length(groups)) {
            errorCondition(recall = PlotMeansDoE.menu, message = gettextRcmdr("More than two factors selected."))
            return()
        }
        if (0 == length(response)) {
            errorCondition(recall = PlotMeansDoE.menu, message = gettextRcmdr("No response variable selected."))
            return()
        }
        .activeDataSet <- ActiveDataSet()
        error.bars <- tclvalue(errorBarsVariable)
        level <- if (error.bars == "conf.int")
            paste(", level=", tclvalue(levelVariable), sep = "")
        else ""
        if (length(groups) == 1)
            doItAndPrint(paste("plotMeans(", .activeDataSet,
                "$", response, ", ", .activeDataSet, "$", groups[1],
                ", error.bars=\"", error.bars, "\"", level, ", main=\"Means from ", .activeDataSet,
                "\", ylab=\"", response,"\", xlab=\"",groups[1],"\")",
                sep = ""))
        else {
            if (eval(parse(text = paste("length(levels(", .activeDataSet,
                "$", groups[1], ")) < length(levels(", .activeDataSet,
                "$", groups[2], "))", sep = ""))))
                groups <- rev(groups)
            doItAndPrint(paste("plotMeans(", .activeDataSet,
                "$", response, ", ", .activeDataSet, "$", groups[1],
                ", ", .activeDataSet, "$", groups[2], ", error.bars=\"",
                error.bars, "\"", level, ", main=\"Means from ", .activeDataSet,
                "\", ylab=\"", response, "\", xlab=\"",groups[1],"\", legend.lab=\"", groups[2],"\")", sep = ""))
        }
        activateMenus()
        tkfocus(CommanderWindow())
    }
  #  optionsFrame <- tkframe(top)
    errorBarsVariable <- tclVar("none")
  #  seButton <- ttkradiobutton(optionsFrame, variable = errorBarsVariable,
  #      value = "se")
  #  sdButton <- ttkradiobutton(optionsFrame, variable = errorBarsVariable,
  #      value = "sd")
  #  confIntButton <- ttkradiobutton(optionsFrame, variable = errorBarsVariable,
  #      value = "conf.int")
  #  noneButton <- ttkradiobutton(optionsFrame, variable = errorBarsVariable,
  #      value = "none")
  #  levelVariable <- tclVar("0.95")
  #  levelEntry <- ttkentry(optionsFrame, width = "6", textvariable = levelVariable)
  #  buttonsFrame <- tkframe(top)
    OKCancelHelp(helpSubject = "PlotMeansDoE.menu")
    tkgrid(getFrame(groupBox), getFrame(responseBox), sticky = "nw")
  #  tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Error Bars"),
  #      fg = "blue"), sticky = "w")
  #  tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Standard errors")),
  #      seButton, sticky = "w")
  #  tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Standard deviations")),
  #      sdButton, sticky = "w")
  #  tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Confidence intervals")),
  #      confIntButton, labelRcmdr(optionsFrame, text = gettextRcmdr("   Level of confidence:")),
  #      levelEntry, sticky = "w")
  #  tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("No error bars")),
  #      noneButton, sticky = "w")
  #  tkgrid(optionsFrame, columnspan = 2, sticky = "w")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 3, columns = 2)
}
