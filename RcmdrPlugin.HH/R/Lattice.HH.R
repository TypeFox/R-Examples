## The following function by Richard Heiberger, with small modifications by J. Fox
## with more modifications by Richard Heiberger.
## 2008-01-03 added conditions, layout, and multiple colors
## Extension to more lattice functions by Richard Heiberger 2009-06-18
## Extensions to barchart function by Richard Heiberger 2009-11-29
## 2012-08-19 rmh added memory to the dialogs, using John Fox's getDialog and putDialog functions


Xyplot.HH2 <- function() {
  Library("lattice")
  defaults <- list(initial.predictor = NULL, initial.response = NULL,
                   initial.auto.key = 1, initial.outer = 0,
                   initial.x.relation = "same", initial.y.relation = "same",
                   initial.layoutColumns = "", initial.layoutRows = "",
                   initial.conditions = FALSE,
                   initial.groups = FALSE,
                   initial.points = 1, initial.lines = 0,
                   initial.functionName="xyplot",
                   initial.panelName="DEFAULT",
                   initial.x.log="Linear",
                   initial.y.log="Linear",
                   initial.horizontal="NULL")
  dialog.values <- getDialog("Xyplot.HH2", defaults)

  initializeDialog(title=gettextRcmdr("Lattice Plot"))
  predictorFrame <- tkframe(top)
  predictorBox <- variableListBox(predictorFrame, Variables(),
                                  title=gettextRcmdr("Explanatory variables (pick one or more)"),
                                  selectmode="multiple",
                                  initialSelection = varPosn (dialog.values$initial.predictor, "all"))

  response.if <-
    (length(dialog.values$initial.response) == 1 &&
     dialog.values$initial.response == FALSE) ||
       length(dialog.values$initial.response) == 0
  responseBox <- variableListBox(predictorFrame, Variables(),
                                 title=gettextRcmdr("Response variables (pick zero or more)"),
                                 selectmode="multiple",
                                 initialSelection = if (response.if) FALSE else
                                 varPosn (dialog.values$initial.response, "all"))
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
                               varPosn (dialog.values$initial.groups, "factor"))

  functionFrame <- tkframe(top)
  functionNameBox <- variableListBox(functionFrame, latticeFunctions(),
                                     title=gettextRcmdr("Pick lattice function"),
                                     selectmode="single",
                                     initialSelection=match(dialog.values$initial.functionName, latticeFunctions(), 0) - 1,
                                     listHeight=5)
  panelNameBox <- variableListBox(functionFrame, latticePanelFunctions(),
                                  title=gettextRcmdr("Pick panel function (or accept DEFAULT)"),
                                  selectmode="single",
                                  initialSelection=match(dialog.values$initial.panelName, latticePanelFunctions(), 0) - 1,
                                  listHeight=5)

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

  logFrame <- tkframe(top)
  radioButtons(window=logFrame,
               name="x.log",
               buttons=c("Linear", "Ten", "Two", "e"),
               labels=gettextRcmdr(c("Linear","10","2","e")),
               title=gettextRcmdr("X-Axis Linear or Log Scale"),
               initialValue = dialog.values$initial.x.log)
  radioButtons(window=logFrame,
               name="y.log",
               buttons=c("Linear", "Ten", "Two", "e"),
               labels=gettextRcmdr(c("Linear","10","2","e")),
               title=gettextRcmdr("Y-Axis Linear or Log Scale"),
               initialValue = dialog.values$initial.y.log)

  scalarsFrame <- tkframe(top)
  layoutColumnsVar <- tclVar(dialog.values$initial.layoutColumns)
  layoutColumnsEntry <- tkentry(scalarsFrame, width="6", textvariable=layoutColumnsVar)
  layoutRowsVar <- tclVar(dialog.values$initial.layoutRows)
  layoutRowsEntry <- tkentry(scalarsFrame, width="6", textvariable=layoutRowsVar)

  radioButtons(name="horizontal",
               buttons=gettextRcmdr(c("NULL","TRUE","FALSE")),
               labels=c("DEFAULT", "horizontal", "vertical"),
               title=gettextRcmdr("Horizontal"),
               initialValue = dialog.values$initial.horizontal)

  checkBoxes(frame="typeFrame",
             boxes=c("points", "lines"),
             initialValues=c(dialog.values$initial.points, dialog.values$initial.lines),
             labels=gettextRcmdr(c("Points", "Lines")))

  onOK <- function() {
    predictor <- getSelection(predictorBox)
    response <- getSelection(responseBox)
    conditions <- getSelection(conditionsBox)
    groups <- getSelection(groupsBox)
    functionName <- getSelection(functionNameBox)
    panelName <- getSelection(panelNameBox)
    closeDialog()

    if (0 == length(response) && functionName == "xyplot") {
      errorCondition(recall=Xyplot.HH2, message=gettextRcmdr("At least one response variable must be selected."))
      return()
    }
    if (0 == length(predictor)) {
      errorCondition(recall=Xyplot.HH2, message=gettextRcmdr("At least one explanatory variable must be selected."))
      return()
    }
    auto.key <- ("1" == tclvalue(auto.keyVariable))
    outer    <- ("1" == tclvalue(outerVariable))
    x.relation <- as.character(tclvalue(x.relationVariable))
    y.relation <- as.character(tclvalue(y.relationVariable))

    x.log <- switch(as.character(tclvalue(x.logVariable)), Linear=FALSE, Ten=10, Two=2, e="'e'")
    y.log <- switch(as.character(tclvalue(y.logVariable)), Linear=FALSE, Ten=10, Two=2, e="'e'")

    horizontal <- tclvalue(horizontalVariable)

    layoutColumns  <- as.numeric(tclvalue(layoutColumnsVar))
    layoutRows     <- as.numeric(tclvalue(layoutRowsVar))

    points <- ("1" == tclvalue(pointsVariable))
    lines  <- ("1" == tclvalue(linesVariable))

    putDialog ("Xyplot.HH2", list(initial.predictor = predictor, initial.response = response,
                                  initial.auto.key = auto.key, initial.outer = outer,
                                  initial.x.relation = x.relation,
                                  initial.y.relation = y.relation,
                                  initial.layoutColumns = tclvalue(layoutColumnsVar),
                                  initial.layoutRows = tclvalue(layoutRowsVar),
                                  initial.conditions = if (length(conditions) != 0) conditions else FALSE,
                                  initial.groups = if (length(groups) != 0) groups else FALSE,
                                  initial.points = points, initial.lines = lines,
                                  initial.functionName=functionName,
                                  initial.panelName=panelName,
                                  initial.x.log=as.character(tclvalue(x.logVariable)),
                                  initial.y.log=as.character(tclvalue(y.logVariable)),
                                  initial.horizontal=horizontal)
               )

    layout.command <- ""
    number.na <- is.na(layoutColumns) + is.na(layoutRows)

    if (number.na==1) {
      errorCondition(recall=Xyplot.HH2,
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
      errorCondition(recall=Xyplot.HH2,
                     message=gettextRcmdr("Choose at least one of points or lines."))
      return()
    }

    type.command <-  paste(", type=",
                           deparse(c("p"[points], "l"[lines])),
                           sep="")

    panel.command <-
      if (panelName == "DEFAULT")
        ""
      else
        paste(", panel=", panelName, sep="")

    horizontal.command <- if (horizontal == "NULL") "" else paste(", horizontal=", horizontal, sep="")

    functionFormula <- if (functionName == "splom")
      splomFormula(predictor, .activeDataSet)
    else
      usualFormula(functionName, response, predictor, .activeDataSet)

    data.command <- if ((functionName == "splom") && (length(groups) == 0 && length(conditions) == 0))
      ""
    else
      paste(", data=", .activeDataSet)
    xyplot.command <- paste(functionFormula,
                            if (length(conditions) > 0)
                            paste(" |",
                                  paste(conditions, collapse=" + ")
                                  ) else "",
                            if (outer) ", outer=TRUE",
                            conditions.command,
                            groups.command,
                            data.command,
                            type.command,
                            horizontal.command,
                            ", par.settings=simpleTheme(pch=16)",
                            panel.command,
                            if (functionName == "barchart"  || panelName == "panel.barchart") ", origin=0" else "",
                            if (auto.key) ", auto.key=list(border=TRUE)" else "",
                            paste(", scales=list(x=list(relation='",
                                  x.relation,
                                  "', log=", x.log, "), y=list(relation='",
                                  y.relation,
                                  "', log=", y.log, "))", sep=""),
                            ')', sep="")

    doItAndPrint(xyplot.command)

    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="xyplot", reset = "Xyplot.HH2")
  tkgrid(getFrame(predictorBox), getFrame(responseBox),
         columnspan=1, sticky="w")
  tkgrid(predictorFrame, sticky="w")
  tkgrid(getFrame(conditionsBox),
         tklabel(cgFrame, text=gettextRcmdr("           ")),
         getFrame(groupsBox),
         columnspan=1, sticky="w")
  tkgrid(cgFrame, sticky="w")

  tkgrid(tklabel(top, text=gettextRcmdr("Function"), fg="blue"), sticky="w")
  tkgrid(getFrame(functionNameBox), getFrame(panelNameBox), columnspan=1, sticky="w")
  tkgrid(functionFrame, sticky="w")

  tkgrid(tklabel(top, text=gettextRcmdr("Options"), fg="blue"), sticky="w")
  tkgrid(optionsFrame, sticky="w")
  tkgrid(tklabel(top, text=gettextRcmdr("Plot Type (one or both)"), fg="blue"), sticky="w")
  tkgrid(typeFrame, sticky="w")

  tkgrid(x.relationFrame, y.relationFrame, columnspan=2, sticky="w")
  tkgrid(relationFrame, sticky="w")

  tkgrid(x.logFrame, y.logFrame, columnspan=2, sticky="w")
  tkgrid(logFrame, sticky="w")

  tkgrid(horizontalFrame, sticky="w")

  tkgrid(tklabel(top, text=gettextRcmdr("Layout"), fg="blue"),
         sticky="w")
  tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("number of columns:")), layoutColumnsEntry, sticky="w")
  tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("number of rows:")), layoutRowsEntry, sticky="w")
  tkgrid(scalarsFrame, sticky="w")

  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  dialogSuffix(rows=7, columns=2)
}


latticeFunctions <- function() {
  putRcmdr("functions",
           c("xyplot", "bwplot", "splom", "barchart", "dotplot", "likert")
           )
  getRcmdr("functions")
}



latticePanelFunctions <- function() {
  putRcmdr("panels",
           c("DEFAULT",
             "panel.xyplot",
             "panel.bwplot",
             "panel.bwplot.intermediate.hh",
             "panel.pairs",
             "panel.splom",
             "panel.barchart",
             "panel.dotplot"))
  getRcmdr("panels")
}

splomFormula <- function(predictor, data.frame.name) {
  expl.subscr <- paste('"', paste(predictor, collapse='","'), '"', sep="")
  paste("splom( ~ ", data.frame.name, "[, c(", expl.subscr, ")]", sep="")
}

usualFormula <- function(functionName, response, predictor, data.frame.name) {
  paste(functionName, "(",
        paste(response, collapse=" + "),
        " ~ ",
        paste(predictor, collapse=" + "),
        sep="")
}

## source("~/HH-R.package/RcmdrPlugin.HH2/R/Lattice.HH.R")
## source("~/WindowsC/HOME/rmh/HH-R.package/RcmdrPlugin.HH2/R/Lattice.HH.R")

