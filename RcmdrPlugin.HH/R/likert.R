listAllLikertCapable <- function (envir = .GlobalEnv, ...)
{
    objects <- ls(envir = envir, ...)
    if (length(objects) == 0)
        NULL
    else objects[sapply(objects, function(x) is.likertCapable(get(x)))]
}


PlotLikertDialog <- function() {
  initializeDialog(title=gettextRcmdr("Likert Plots"))
  allLikertCapables <- listAllLikertCapable()
  activeDatasetPosition <-  match(ActiveDataSet(), allLikertCapables, 0)
  if (activeDatasetPosition > getRcmdr("variable.list.height")) {
    allLikertCapables <- c(allLikertCapables[activeDatasetPosition],
                           allLikertCapables[-activeDatasetPosition])
    activeDatasetPosition <- 1
  }
  likertFrame <- tkframe(top)
  likertBox <- variableListBox(likertFrame, allLikertCapables,
                               title=gettextRcmdr("Likert Capable Objects (pick one)"),
                               selectmode="single",
                               initialSelection=activeDatasetPosition-1)

  mainVar   <- tclVar("")
  mainEntry <- tkentry(likertFrame, width="30", textvariable=mainVar)
  LikertPlotNameVar   <- tclVar("LikertPlot")
  LikertPlotNameEntry <- tkentry(likertFrame, width="16", textvariable=LikertPlotNameVar)
  boxWidthNumberVar   <- tclVar("")
  boxWidthNumberEntry <- tkentry(likertFrame, width="16", textvariable=boxWidthNumberVar)
  boxWidthUnitVar   <- tclVar("mm")
  boxWidthUnitEntry <- tkentry(likertFrame, width="16", textvariable=boxWidthUnitVar)

  ## horizontalFrame <- tkframe(top)
  ## radioButtons(name="horizontal",
  ##              buttons=c("TRUE", "FALSE"),
  ##              labels=gettextRcmdr(c("Horizontal", "Vertical")),
  ##              title=gettextRcmdr("Horizontal of Bars"))

  checkBoxes(frame="optionsFrame",
             boxes=c("horizontal","as.percent","positive.order"),
             initialValues=c(1,0,0),
             labels=gettextRcmdr(c("Horizontal Bars","Plot Percents","Sort by Total Positive")))

  layoutXVar        <- tclVar("")
  layoutXEntry      <- tkentry(likertFrame, width="8", textvariable=layoutXVar)
  layoutYVar        <- tclVar("")
  layoutYEntry      <- tkentry(likertFrame, width="8", textvariable=layoutYVar)
  ReferenceZeroVar   <- tclVar("")
  ReferenceZeroEntry <- tkentry(likertFrame, width="16", textvariable=ReferenceZeroVar)
  BrewerPaletteVar   <- tclVar("")
  BrewerPaletteEntry <- tkentry(likertFrame, width="16", textvariable=BrewerPaletteVar)

  onOK <- function() {
                                        #on.exit(recover())
    table <- getSelection(likertBox)
    closeDialog()
    if (0 == length(table)) {
      errorCondition(recall=PlotLikertDialog,
                     message=gettextRcmdr("Exactly one Likert Capable Object must be selected."))
      return()
    }

    main  <- tclvalue(mainVar)
    LikertPlotName  <- tclvalue(LikertPlotNameVar)
    if (nchar(LikertPlotName)==0) LikertPlotName <- "LikertPlot"
    boxWidthValue  <- tclvalue(boxWidthNumberVar)
    boxWidthUnitValue  <- tclvalue(boxWidthUnitVar)
##    horizontal <- as.logical(tclvalue(horizontalVariable))
    horizontal <- ("1" == tclvalue(horizontalVariable))
    as.percent <- ("1" == tclvalue(as.percentVariable))
    positive.order <- ("1" == tclvalue(positive.orderVariable))
    layoutXValue  <- tclvalue(layoutXVar)
    layoutYValue  <- tclvalue(layoutYVar)
    ReferenceZeroValue  <- tclvalue(ReferenceZeroVar)
    BrewerPaletteName  <- tclvalue(BrewerPaletteVar)

    command1 <- paste(LikertPlotName, " <- plot.likert(", table,
                      if (nchar(main) != 0) paste(", main='", main, "'", sep=""),
                      if (nchar(boxWidthValue) != 0) paste(
                                 ", box.width=unit(", boxWidthValue,
                                 ",'", boxWidthUnitValue, "')", sep=""),
                      if (!horizontal)
                         paste(", horizontal=FALSE, auto.key=list(reverse=TRUE, space='right', columns=1, padding.text=2)"),
                      if (as.percent) paste(", as.percent=TRUE"),
                      if (positive.order) paste(", positive.order=TRUE"),
                      if (nchar(layoutYValue) != 0 &&
                          nchar(layoutYValue) != 0 )
                             paste(", layout=c(", layoutXValue, ",", layoutYValue, ")", sep=""),
                      if (nchar(ReferenceZeroValue) != 0) paste(
                                 ", ReferenceZero=", ReferenceZeroValue, sep=""),
                      if (nchar(BrewerPaletteName) != 0) paste(
                                 ", BrewerPaletteName='", BrewerPaletteName, "'", sep=""),
                      ")", sep="")
    doItAndPrint(command1)
    doItAndPrint(LikertPlotName)
    activateMenus()
    tkfocus(CommanderWindow())
    if (version$os == "mingw32") justDoIt("bringToTop()")
  }
  OKCancelHelp(helpSubject="plot.likert")
  tkgrid(getFrame(likertBox),
         columnspan=1, sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("Main title:")), mainEntry, sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("LikertPlot Name:")), LikertPlotNameEntry, sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("box.width:")), boxWidthNumberEntry, sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("box.width unit:")), boxWidthUnitEntry, sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("layout x:")),         layoutXEntry          , sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("layout y:")),         layoutYEntry          , sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("ReferenceZero:")),  ReferenceZeroEntry   , sticky="w")
  tkgrid(tklabel(likertFrame, text=gettextRcmdr("RColorBrewerPalette:")),  BrewerPaletteEntry   , sticky="w")
  tkgrid(likertFrame, sticky="w")
##  tkgrid(horizontalFrame, sticky="w")
  tkgrid(optionsFrame, sticky="w")


  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  dialogSuffix(rows=7, columns=2)
}

## source("c:/HOME/rmh/HH-R.package/RcmdrPlugin.HH/R/likert.R")
## source("x:/HOME/rmh/HH-R.package/RcmdrPlugin.HH/R/likert.R")
## source("/Users/rmh/WindowsC/HOME/rmh/HH-R.package/RcmdrPlugin.HH/R/likert.R")
