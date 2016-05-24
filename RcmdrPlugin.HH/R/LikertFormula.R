## The following function by Richard Heiberger, with small modifications by J. Fox
## with more modifications by Richard Heiberger.
## 2008-01-03 added conditions, layout, and multiple colors
## Extension to more lattice functions by Richard Heiberger 2009-06-18
## Extensions to barchart function by Richard Heiberger 2009-11-29
## 2012-08-19 rmh added memory to the dialogs, using John Fox's getDialog and putDialog functions
## 2012-12-27 rmh LikertFormula based on Xyplot.HH2 and PlotLikertDialog

LikertFormula <- function() {
  Library("lattice")

  defaults <- list(initial.predictor = NULL,
                   initial.response = NULL,
                   initial.layoutColumns = "",
                   initial.layoutRows = "",
                   initial.main = "",
                   initial.LikertPlotName = "LikertPlot",
                   initial.boxWidthValue = "",
                   initial.boxWidthUnitValue = "mm",
                   initial.ReferenceZeroValue = "",
##                   initial.BrewerPaletteName = "",
                   initial.xlabName = "",
                   initial.ylabName = "",
                   initial.subName = "",
                   initial.resizewidthName = "",
                   initial.resizeheightName = "",
                   initial.conditions = 0,
                   initial.horizontal = 1,
                   initial.as.percent = 0,
                   initial.rowsort = "levels",
                   initial.reverse = 1,
                   initial.relationx = "same",
                   initial.relationy = "same",
                   initial.strip = 1,
                   initial.strip.left = 0,
                   initial.useOuterStrips=0,
                   initial.combineLimits=0)
  dialog.values <- getDialog("LikertFormula", defaults)

  initializeDialog(title=gettextRcmdr("Likert Formula"))

  ngad <- names(get(ActiveDataSet()))  ## intermediate variable used to maintain variables in original order

  predictorFrame <- tkframe(top)

  response.if <-
    (length(dialog.values$initial.response) == 1 &&
     dialog.values$initial.response == 0) ||
       length(dialog.values$initial.response) == 0
  responseBox <- variableListBox(predictorFrame, ngad[ngad %in% Factors()],
                                 title=gettextRcmdr("Response variables (pick zero or one) ~ "),
                                 selectmode="single",
                                 initialSelection = if (response.if) FALSE else
                                 varPosnOriginal(dialog.values$initial.response, "factor"))
  predictorBox <- variableListBox(predictorFrame, ngad[ngad %in% Numeric()],
                                  title=gettextRcmdr("Explanatory variables (pick one or more) | "),
                                  selectmode="multiple",
                                  initialSelection = varPosnOriginal(dialog.values$initial.predictor, "numeric"))

  ## cgFrame <- tkframe(top)
  conditions.if <-
    length(dialog.values$initial.conditions) == 1 &&
      dialog.values$initial.conditions == FALSE
  conditionsBox <- variableListBox(predictorFrame, ngad[ngad %in% Factors()],   ## cgFrame
                                   title=gettextRcmdr("Conditions (pick zero or more)"),
                                   selectmode="multiple",
                                   initialSelection=if (conditions.if) FALSE else
                                   varPosnOriginal(dialog.values$initial.conditions, "factor"))

  checkBoxes(frame="optionsFrame",
             boxes=c("horizontal","as.percent", "reverse","strip","strip.left","useOuterStrips","combineLimits"),
             initialValues=c(
               dialog.values$initial.horizontal,
               dialog.values$initial.as.percent,
               dialog.values$initial.reverse, ##c(1, 0, 1),
               dialog.values$initial.strip, ## 1
               dialog.values$initial.strip.left, ## 0
               dialog.values$initial.useOuterStrips, ## 0
               dialog.values$initial.combineLimits), ## 0
             labels=gettextRcmdr(c("Horizontal Bars", "Plot Percents", "Reverse Row Order",
               "Strip","Strip Left", "useOuterStrips", "combineLimits")))

  rowsortFrame <- tkframe(top)
  relationxFrame <- tkframe(top)
  relationyFrame <- tkframe(top)
  radioButtons(name="rowsort",
               buttons=c("levels", "positive.order", "data.order"),
               labels=gettextRcmdr(c("levels", "positive.order", "data.order")),
               title=gettextRcmdr("Row Sort Sequence"),
               initialValue = dialog.values$initial.rowsort)

  radioButtons(name="relationx",
               buttons=c("same", "free", "sliced"),
               labels=gettextRcmdr(c("same", "free", "sliced")),
               title=gettextRcmdr("x relation"),
               initialValue = dialog.values$initial.relationx)
  radioButtons(name="relationy",
               buttons=c("same", "free", "sliced"),
               labels=gettextRcmdr(c("same", "free", "sliced")),
               title=gettextRcmdr("y relation"),
               initialValue = dialog.values$initial.relationy)

  scalarsFrame <- tkframe(top)
  layoutColumnsVar <- tclVar(dialog.values$initial.layoutColumns)
  layoutColumnsEntry <- tkentry(scalarsFrame, width="6", textvariable=layoutColumnsVar)
  layoutRowsVar <- tclVar(dialog.values$initial.layoutRows)
  layoutRowsEntry <- tkentry(scalarsFrame, width="6", textvariable=layoutRowsVar)


  likertFrame <- tkframe(top)
  mainVar   <- tclVar(dialog.values$initial.main)
  mainEntry <- tkentry(likertFrame, width="30", textvariable=mainVar)
  LikertPlotNameVar   <- tclVar(dialog.values$initial.LikertPlotName)
  LikertPlotNameEntry <- tkentry(likertFrame, width="16", textvariable=LikertPlotNameVar)
  boxWidthNumberVar   <- tclVar(dialog.values$initial.boxWidthValue)
  boxWidthNumberEntry <- tkentry(likertFrame, width="16", textvariable=boxWidthNumberVar)
  boxWidthUnitVar   <- tclVar(dialog.values$initial.boxWidthUnitValue)
  boxWidthUnitEntry <- tkentry(likertFrame, width="16", textvariable=boxWidthUnitVar)

  ReferenceZeroVar   <- tclVar(dialog.values$initial.ReferenceZeroValue)
  ReferenceZeroEntry <- tkentry(likertFrame, width="16", textvariable=ReferenceZeroVar)
##  BrewerPaletteVar   <- tclVar(dialog.values$initial.BrewerPaletteName)
##  BrewerPaletteEntry <- tkentry(likertFrame, width="16", textvariable=BrewerPaletteVar)
  xlabVar   <- tclVar(dialog.values$initial.xlabName)
  xlabEntry <- tkentry(likertFrame, width="16", textvariable=xlabVar)
  ylabVar   <- tclVar(dialog.values$initial.ylabName)
  ylabEntry <- tkentry(likertFrame, width="16", textvariable=ylabVar)
  subVar   <- tclVar(dialog.values$initial.subName)
  subEntry <- tkentry(likertFrame, width="16", textvariable=subVar)
  resizewidthVar   <- tclVar(dialog.values$initial.resizewidthName)
  resizewidthEntry <- tkentry(likertFrame, width="16", textvariable=resizewidthVar)
  resizeheightVar   <- tclVar(dialog.values$initial.resizeheightName)
  resizeheightEntry <- tkentry(likertFrame, width="16", textvariable=resizeheightVar)

  ## relationFrame <- tkframe(top)
  ## relationxVar   <- tclVar(dialog.values$initial.relationxName)
  ## relationxEntry <- tkentry(relationFrame, width="16", textvariable=relationxVar)
  ## relationyVar   <- tclVar(dialog.values$initial.relationyName)
  ## relationyEntry <- tkentry(relationFrame, width="16", textvariable=relationyVar)
  ## tkgrid(getFrame(relationxEntry), getFrame(relationyEntry),
  ##        columnspan=1, sticky="w")
  ## tkgrid(relationFrame, sticky="w")

  onOK <- function() {
    predictor <- getSelection(predictorBox)
    response <- getSelection(responseBox)
    conditions <- getSelection(conditionsBox)
    closeDialog()

    ## if (0 == length(response)) {
    ##   errorCondition(recall=LikertFormula, message=gettextRcmdr("One response variable must be selected."))
    ##   return()
    ## }
    if (0 == length(predictor)) {
      errorCondition(recall=LikertFormula, message=gettextRcmdr("At least one explanatory variable must be selected."))
      return()
    }
    horizontal <-     ("1" == tclvalue(horizontalVariable))
    as.percent <-     ("1" == tclvalue(as.percentVariable))
    reverse    <-     ("1" == tclvalue(reverseVariable))
    strip      <-     ("1" == tclvalue(stripVariable))
    strip.left <-     ("1" == tclvalue(strip.leftVariable))
    useOuterStrips <-     ("1" == tclvalue(useOuterStripsVariable))
    combineLimits <-     ("1" == tclvalue(combineLimitsVariable))
    rowsort <- as.character(tclvalue(rowsortVariable))

    layoutColumns  <- as.numeric(tclvalue(layoutColumnsVar))
    layoutRows     <- as.numeric(tclvalue(layoutRowsVar))

    main  <- tclvalue(mainVar)
    LikertPlotName  <- tclvalue(LikertPlotNameVar)
    if (nchar(LikertPlotName)==0) LikertPlotName <- "LikertPlot"
    boxWidthValue  <- tclvalue(boxWidthNumberVar)
    boxWidthUnitValue  <- tclvalue(boxWidthUnitVar)
    ReferenceZeroValue  <- tclvalue(ReferenceZeroVar)
##    BrewerPaletteName  <- tclvalue(BrewerPaletteVar)
    xlabName  <- tclvalue(xlabVar)
    ylabName  <- tclvalue(ylabVar)
    subName   <- tclvalue(subVar)
    resizewidthName   <- tclvalue(resizewidthVar)
    resizeheightName   <- tclvalue(resizeheightVar)
    relationx <- as.character(tclvalue(relationxVariable))
    relationy <- as.character(tclvalue(relationyVariable))

    putDialog ("LikertFormula", list(initial.predictor = predictor,
                                     initial.response = response,
                                     initial.layoutColumns = tclvalue(layoutColumnsVar),
                                     initial.layoutRows = tclvalue(layoutRowsVar),
                                     initial.main = main,
                                     initial.LikertPlotName = LikertPlotName,
                                     initial.boxWidthValue = boxWidthValue,
                                     initial.boxWidthUnitValue = boxWidthUnitValue,
                                     initial.ReferenceZeroValue = ReferenceZeroValue,
##                                     initial.BrewerPaletteName = BrewerPaletteName,
                                     initial.xlabName = xlabName,
                                     initial.ylabName = ylabName,
                                     initial.subName = subName,
                                     initial.resizewidthName = resizewidthName,
                                     initial.resizeheightName = resizeheightName,
                                     initial.conditions = if (length(conditions) != 0) conditions else FALSE,
                                     initial.horizontal = horizontal,
                                     initial.as.percent = as.percent,
                                     initial.rowsort=rowsort,
                                     initial.reverse = reverse,
                                     initial.relationx = relationx,
                                     initial.relationy = relationy,
                                     initial.strip = strip,
                                     initial.strip.left = strip.left,
                                     initial.useOuterStrips = useOuterStrips,
                                     initial.combineLimits = combineLimits
                                     )
               )

    layout.command <- ""
    number.na <- is.na(layoutColumns) + is.na(layoutRows)

    if (number.na==1) {
      errorCondition(recall=LikertFormula,
                     message=gettextRcmdr("Both or neither layout values must be numbers."))
      return()
    }
    if (number.na==0) layout.command <- paste(", layout=", deparse(c(layoutColumns, layoutRows)), sep="")

    .activeDataSet <- ActiveDataSet()

    functionFormula <- LikertFormulaConstruct("likert", response, predictor)

    data.command <- paste(", data=", .activeDataSet)

    likert.command <- paste(
      LikertPlotName, " <- ",
      functionFormula,
      if (length(conditions) > 0)
      paste(" |",
            paste(conditions, collapse=" + ")
            ) else "",
      layout.command,
      data.command,
      if (!horizontal)
      ", horizontal=FALSE, auto.key=list(reverse=TRUE, space='right', columns=1, padding.text=2)",
      if (as.percent) ", as.percent=TRUE",
      if (as.percent && length(conditions) > 1) ", rightAxis=FALSE",
      if (rowsort=="positive.order") ", positive.order=TRUE",
      if (rowsort=="data.order") ", data.order=TRUE",
      if (!reverse) ", reverse=FALSE",
      if (nchar(main) != 0) paste(", main='", main, "'", sep=""),
      if (nchar(boxWidthValue) != 0) paste(
                 ", box.width=unit(", boxWidthValue,
                 ",'", boxWidthUnitValue, "')", sep=""),
      ", par.settings=simpleTheme(pch=16)",
      if (nchar(ReferenceZeroValue) != 0) paste(
                 ", ReferenceZero=", ReferenceZeroValue, sep=""),
##      if (nchar(BrewerPaletteName) != 0) paste(
##                 ", BrewerPaletteName='", BrewerPaletteName, "'", sep=""),
      if (nchar(xlabName) != 0) {
        if (xlabName == "NULL")
          ", xlab=NULL"
        else
          paste(", xlab='", xlabName, "'", sep="")
      },
      if (nchar(ylabName) != 0) {
        if (ylabName == "NULL")
          ", ylab=NULL"
        else
          paste(", ylab='", ylabName, "'", sep="")
      },
      if (nchar(subName) != 0) paste(
                 ", sub='", subName, "'", sep=""),
      ", scales=list(x=list(relation='", relationx, "'), y=list(relation='", relationy, "'))",
      if (nchar(resizewidthName)  != 0) {
        if (resizewidthName == "rowSums")
          resizewidthName <-  paste("rowSums(", .activeDataSet, "[,", deparse(predictor, width.cutoff=500), "])", sep="")
        if (resizewidthName != "nrow") paste(", w.resizePanels=", resizewidthName, sep="")
      },
      if (nchar(resizeheightName) != 0) {
        if (resizeheightName == "rowSums")
          resizeheightName <- paste("rowSums(", .activeDataSet, "[,", deparse(predictor, width.cutoff=500), "])", sep="")
        if (resizeheightName != "nrow") paste(", h.resizePanels=", resizeheightName, sep="")
      },

      paste(", strip=", strip, sep=""),
      paste(", strip.left=", strip.left, sep=""),
      ")",
      sep="")

    doItAndPrint(likert.command)
    if (combineLimits || useOuterStrips) {
      StripsLimits <- LikertPlotName
      if (useOuterStrips) StripsLimits <- paste("useOuterStrips(", StripsLimits, ")", sep="")
      if (combineLimits) StripsLimits <- paste("combineLimits(", StripsLimits, ")", sep="")
      StripsLimits <- paste(LikertPlotName, "<-", StripsLimits)
      doItAndPrint(StripsLimits)
    }
    doItAndPrint(LikertPlotName)
    activateMenus()

    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="likert", reset = "LikertFormula")
  tkgrid(getFrame(responseBox), getFrame(predictorBox), getFrame(conditionsBox),
         columnspan=1, sticky="w")
  tkgrid(predictorFrame, sticky="w")
  ## tkgrid(getFrame(conditionsBox),
  ##        columnspan=1, sticky="w")
  ## tkgrid(cgFrame, sticky="w")

  tkgrid(tklabel(top, text=gettextRcmdr("Options"), fg="blue"), sticky="w")
  tkgrid(optionsFrame, sticky="w")

  tkgrid(rowsortFrame, sticky="w")
  tkgrid(relationxFrame, sticky="w")
  tkgrid(relationyFrame, sticky="w")

  tkgrid(tklabel(top, text=gettextRcmdr("Labels"), fg="blue"), sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("Main title:")), mainEntry, sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("LikertPlot Name:")), LikertPlotNameEntry, sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("box.width:")), boxWidthNumberEntry, sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("box.width unit:")), boxWidthUnitEntry, sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("ReferenceZero:")),  ReferenceZeroEntry   , sticky="w")
##  tkgrid(tklabel(likertFrame, text=gettextRcmdr("RColorBrewerPalette:")),  BrewerPaletteEntry   , sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("x label:")),  xlabEntry   , sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("y label:")),  ylabEntry   , sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("subtitle:")),  subEntry   , sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("Resize Width:")),  resizewidthEntry   , sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("Resize Height:")),  resizeheightEntry   , sticky="w")
  tkgrid(likertFrame, sticky="w")

  tkgrid(tklabel(top, text=gettextRcmdr("Layout"), fg="blue"), sticky="w")
  tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("number of columns:")), layoutColumnsEntry, sticky="w")
  tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("number of rows:")), layoutRowsEntry, sticky="w")
  tkgrid(scalarsFrame, sticky="w")

  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  dialogSuffix(rows=7, columns=2)
}


LikertFormulaConstruct <- function(functionName, response, predictor) {
  paste(functionName, "(",
        paste(response, collapse=" + "),
        " ~ ",
        paste(predictor, collapse=" + "),
        sep="")
}

varPosnOriginal <-
function (variables, type = c("all", "factor", "numeric", "nonfactor",
    "twoLevelFactor"))
{
    if (is.null(variables))
        return(NULL)
    ngad <- names(get(ActiveDataSet()))  ## intermediate variable used to maintain variables in original order
    type <- match.arg(type)
    vars <- switch(type, all = ngad[ngad %in% Variables()], factor = ngad[ngad %in% Factors()],
        numeric = ngad[ngad %in% Numeric()], nonfactor = setdiff(ngad[ngad %in% Variables()],
            ngad[ngad %in% Factors()]), twoLevelFactor = ngad[ngad %in% TwoLevelFactors()])
    if (any(!variables %in% vars))
        NULL
    else apply(outer(variables, vars, "=="), 1, which) - 1
}
## source("~/HH-R.package/RcmdrPlugin.HH/R/LikertFormula.R")
## source("~/WindowsC/HOME/rmh/HH-R.package/RcmdrPlugin.HH/R/LikertFormula.R")

