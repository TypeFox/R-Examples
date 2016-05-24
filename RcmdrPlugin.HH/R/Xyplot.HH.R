## The following function by Richard Heiberger, with small modifications by J. Fox
## with more modifications by Richard Heiberger.
## 2008-01-03 added conditions, layout, and multiple colors
## 2012-08-19 rmh added memory to the dialogs, using John Fox's getDialog and putDialog functions

Xyplot.HH <- function() {
  Library("lattice")
  defaults <- list(initial.predictor = NULL, initial.response = NULL,
                   initial.auto.key = 1, initial.outer = 0,
                   initial.x.relation = "same", initial.y.relation = "same",
                   initial.layoutColumns = "", initial.layoutRows = "",
                   initial.conditions = FALSE,
                   initial.groups = FALSE,
                   initial.points = 1, initial.lines = 0)
  dialog.values <- getDialog("Xyplot", defaults)
  initializeDialog(title=gettextRcmdr("XY Conditioning Plot"))
  predictorFrame <- tkframe(top)
  predictorBox <-
    variableListBox(predictorFrame, Numeric(),
                    title=gettextRcmdr("Explanatory variables (pick one or more)"),
                    selectmode="multiple",
                    initialSelection = varPosn (dialog.values$initial.predictor, "numeric"))
  responseBox <- variableListBox(predictorFrame, Numeric(),
                                 title=gettextRcmdr("Response variables (pick one or more)"),
                                 selectmode="multiple",
                                 initialSelection = varPosn (dialog.values$initial.response, "numeric"))
  cgFrame <- tkframe(top)
  conditions.if <-
    length(dialog.values$initial.conditions) == 1 &&
      dialog.values$initial.conditions == FALSE
  conditionsBox <- variableListBox(cgFrame, Factors(),
                                   title=gettextRcmdr("Conditions '|' (pick zero or more)"),
                                   selectmode="multiple",
                                   initialSelection=if (conditions.if) FALSE else
                                   varPosn(dialog.values$initial.conditions, "factor"))
  groups.if <-
    length(dialog.values$initial.groups) == 1 &&
      dialog.values$initial.groups == FALSE
  groupsBox <- variableListBox(cgFrame, Factors(),
                               title=gettextRcmdr("Groups 'groups=' (pick zero or more)"),
                               selectmode="multiple",
                               initialSelection=if (groups.if) FALSE else
                               varPosn(dialog.values$initial.groups, "factor"))
  checkBoxes(frame="optionsFrame",
             boxes=c("auto.key", "outer"),
             initialValues=c(dialog.values$initial.auto.key, dialog.values$initial.outer),
             labels=gettextRcmdr(c("Automatically draw key",
               "Different panels for different y ~ x combinations")))
  relationFrame <- tkframe(top)
  radioButtons(window=relationFrame,
               name="x.relation",
               buttons=c("same", "free", "sliced"),
               labels=gettextRcmdr(c("Identical", "Free", "Same range")),
               title=gettextRcmdr("X-Axis Scales in Different Panels"),
               initialValue = dialog.values$initial.x.relation)
  radioButtons(window=relationFrame,
               name="y.relation",
               buttons=c("same", "free", "sliced"),
               labels=gettextRcmdr(c("Identical", "Free", "Same range")),
               title=gettextRcmdr("Y-Axis Scales in Different Panels"),
               initialValue = dialog.values$initial.y.relation)

  scalarsFrame <- tkframe(top)
  layoutColumnsVar <- tclVar(dialog.values$initial.layoutColumns)
  layoutColumnsEntry <- tkentry(scalarsFrame, width="6", textvariable=layoutColumnsVar)
  layoutRowsVar <- tclVar(dialog.values$initial.layoutRows)
  layoutRowsEntry <- tkentry(scalarsFrame, width="6", textvariable=layoutRowsVar)

  checkBoxes(frame="typeFrame",
             boxes=c("points", "lines"),
             initialValues=c(dialog.values$initial.points, dialog.values$initial.lines),
             labels=gettextRcmdr(c("Points", "Lines")))

  onOK <- function() {
    predictor <- getSelection(predictorBox)
    response <- getSelection(responseBox)
    conditions <- getSelection(conditionsBox)
    groups <- getSelection(groupsBox)
    closeDialog()

    if (0 == length(response)) {
      errorCondition(recall=Xyplot.HH,
                     message=gettextRcmdr("At least one response variable must be selected."))
      return()
    }
    if (0 == length(predictor)) {
      errorCondition(recall=Xyplot.HH,
                     message=gettextRcmdr("At least one explanatory variable must be selected."))
      return()
    }
    auto.key <- ("1" == tclvalue(auto.keyVariable))
    outer    <- ("1" == tclvalue(outerVariable))
    x.relation <- as.character(tclvalue(x.relationVariable))
    y.relation <- as.character(tclvalue(y.relationVariable))

    layoutColumns  <- as.numeric(tclvalue(layoutColumnsVar))
    layoutRows     <- as.numeric(tclvalue(layoutRowsVar))

    points <- ("1" == tclvalue(pointsVariable))
    lines  <- ("1" == tclvalue(linesVariable))

    putDialog ("Xyplot", list(initial.predictor = predictor, initial.response = response,
                              initial.auto.key = auto.key, initial.outer = outer,
                              initial.x.relation = x.relation,
                              initial.y.relation = y.relation,
                              initial.layoutColumns = tclvalue(layoutColumnsVar),
                              initial.layoutRows = tclvalue(layoutRowsVar),
                              initial.conditions = if (length(conditions) != 0) conditions else FALSE,
                              initial.groups = if (length(groups) != 0) groups else FALSE,
                              initial.points = points,
                              initial.lines = lines))

    layout.command <- ""
    number.na <- is.na(layoutColumns) + is.na(layoutRows)

    if (number.na==1) {
      errorCondition(recall=Xyplot.HH,
                     message=gettextRcmdr("Both or neither layout values must be numbers."))
      return()
    }
    if (number.na==0) layout.command <- deparse(c(layoutColumns, layoutRows))

    .activeDataSet <- ActiveDataSet()



    conditions.command <-
      if (length(conditions) == 0) {
        if (outer) {
          if (layout.command == "")
            paste(", layout=c(",
                  length(predictor),
                  ",",
                  length(response),
                  ")")
          else
            paste(", layout=", layout.command, sep="")
        }
        else
          if (layout.command != "")
            paste(", layout=", layout.command, sep="")
      }
      else {  ## (length(conditions) > 0)
        if (outer) {
          condition.levels <- prod(sapply(conditions, d.f=get(.activeDataSet),
                                          function(g, d.f) length(levels(d.f[[g]]))))
          if (layout.command != "")
            paste(", layout=", layout.command, sep="")
          else
            paste(", layout=c(",
                  condition.levels,
                  "*",
                  length(predictor),
                  ",",
                  length(response),
                  ")",
                  ## ", between=list(x=c(0,0, 1, 0,0), y=1)",
                  ", between=list(x=c(",
                  paste(rep(c(rep(0, condition.levels-1), 1),
                            length=condition.levels*length(predictor)-1),
                        collapse=","),
                  "), y=1)")
        }
        else
          if (layout.command != "")
            paste(", layout=", layout.command, sep="")
      }


    groups.command <- switch(as.character(length(groups)),
                             "0"="",
                             "1"=paste(", groups=", groups, sep=""),
                             paste(", groups=interaction(",
                                   paste(groups, collapse=","),
                                   ")", sep=""))

    if(!(points || lines)) {
      errorCondition(recall=Xyplot.HH,
                     message=gettextRcmdr("Choose at least one of points or lines."))
      return()
    }

    type.command <- paste(", type=",
                          deparse(c("p"[points], "l"[lines])),
                          sep="")

    xyplot.command <- paste("xyplot(",
                            paste(response, collapse=" + "),
                            " ~ ",
                            paste(predictor, collapse=" + "),
                            if (length(conditions) > 0)
                            paste(" |",
                                  paste(conditions, collapse=" + ")
                                  ) else "",
                            if (outer) ", outer=TRUE",
                            conditions.command,
                            groups.command,
                            type.command,
                            ", pch=16",
                            if (auto.key) ", auto.key=list(border=TRUE), par.settings=simpleTheme(pch=16)" else "",
                            paste(", scales=list(x=list(relation='",
                                  x.relation,
                                  "'), y=list(relation='",
                                  y.relation,
                                  "'))", sep=""),
                            ", data=", .activeDataSet, ")", sep="")
    doItAndPrint(xyplot.command)
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="xyplot", reset = "Xyplot")
  tkgrid(getFrame(predictorBox), getFrame(responseBox),
         columnspan=1, sticky="w")
  tkgrid(predictorFrame, sticky="w")
  tkgrid(getFrame(conditionsBox),
         tklabel(cgFrame, text=gettextRcmdr("           ")),
         getFrame(groupsBox),
         columnspan=1, sticky="w")
  tkgrid(cgFrame, sticky="w")

  tkgrid(tklabel(top, text=gettextRcmdr("Options"), fg="blue"), sticky="w")
  tkgrid(optionsFrame, sticky="w")
  tkgrid(tklabel(top, text=gettextRcmdr("Plot Type (one or both)"), fg="blue"), sticky="w")
  tkgrid(typeFrame, sticky="w")

  tkgrid(x.relationFrame, y.relationFrame, columnspan=2, sticky="w")
  tkgrid(relationFrame, sticky="w")

  tkgrid(tklabel(top, text=gettextRcmdr("Layout"), fg="blue"),
         sticky="w")
  tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("number of columns:")), layoutColumnsEntry, sticky="w")
  tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("number of rows:")), layoutRowsEntry, sticky="w")
  tkgrid(scalarsFrame, sticky="w")

  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  dialogSuffix(rows=6, columns=2)
}
