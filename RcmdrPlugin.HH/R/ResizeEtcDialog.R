ResizeEtcDialog <- function() {
    initializeDialog(title=gettextRcmdr("Resize Panels"))
    resizeFrame <- tkframe(top)

    c.listVar            <- tclVar("")
    condlevelsNameVar    <- tclVar("")
    x.sameVar            <- tclVar("")
    y.sameVar            <- tclVar("")
    layoutVar            <- tclVar("")
    strip.valuesVar      <- tclVar("")
    strip.left.valuesVar <- tclVar("")
    strip.parVar         <- tclVar("")
    strip.left.parVar    <- tclVar("")
    resize.heightVar     <- tclVar("")
    resize.widthVar      <- tclVar("")
    mainVar              <- tclVar("")
    mainMiddleVar        <- tclVar("")

    c.listEntry            <- tkentry(resizeFrame, width="48", textvariable=c.listVar)
    condlevelsNameEntry    <- tkentry(resizeFrame, width="48", textvariable=condlevelsNameVar)
    x.sameEntry            <- tkentry(resizeFrame, width="48", textvariable=x.sameVar)
    y.sameEntry            <- tkentry(resizeFrame, width="48", textvariable=y.sameVar)
    layoutEntry            <- tkentry(resizeFrame, width="48", textvariable=layoutVar)
    strip.valuesEntry      <- tkentry(resizeFrame, width="48", textvariable=strip.valuesVar)
    strip.left.valuesEntry <- tkentry(resizeFrame, width="48", textvariable=strip.left.valuesVar)
    strip.parEntry         <- tkentry(resizeFrame, width="48", textvariable=strip.parVar)
    strip.left.parEntry    <- tkentry(resizeFrame, width="48", textvariable=strip.left.parVar)
    resize.heightEntry     <- tkentry(resizeFrame, width="48", textvariable=resize.heightVar)
    resize.widthEntry      <- tkentry(resizeFrame, width="48", textvariable=resize.widthVar)
    mainEntry              <- tkentry(resizeFrame, width="48", textvariable=mainVar)
    mainMiddleEntry        <- tkentry(resizeFrame, width="48", textvariable=mainMiddleVar)


    onOK <- function() {
      #on.exit(recover())
      c.listValue            <- tclvalue(c.listVar)
      condlevelsNameValue    <- tclvalue(condlevelsNameVar)
      x.sameValue            <- tclvalue(x.sameVar)
      y.sameValue            <- tclvalue(y.sameVar)
      layoutValue            <- tclvalue(layoutVar)
      strip.valuesValue      <- tclvalue(strip.valuesVar)
      strip.left.valuesValue <- tclvalue(strip.left.valuesVar)
      strip.parValue         <- tclvalue(strip.parVar)
      strip.left.parValue    <- tclvalue(strip.left.parVar)
      resize.heightValue     <- tclvalue(resize.heightVar)
      resize.widthValue      <- tclvalue(resize.widthVar)
      mainValue              <- tclvalue(mainVar)
      mainMiddleValue        <- tclvalue(mainMiddleVar)

      closeDialog()

      if (nchar(c.listValue) == 0) {
        errorCondition(recall=ResizeEtcDialog,
                       message=gettextRcmdr("c.list must be specified."))
        return()
      }

      ##command.xmiddle <- "x.middle <- diff(current.panel.limits()$xlim)/2"
      ##doItAndPrint(command.xmiddle)
      if (nchar(mainMiddleValue)==0) mainMiddleValue <- ".5"
      command <- paste("ResizeEtc(",
                                                               paste(  "c.list=",            c.listValue, sep=""),
                       if (nchar(condlevelsNameValue   ) != 0) paste(", condlevelsName='",   condlevelsNameValue   , "'", sep=""),
                       if (nchar(x.sameValue           ) != 0) paste(", x.same=",            x.sameValue           , sep=""),
                       if (nchar(y.sameValue           ) != 0) paste(", y.same=",            y.sameValue           , sep=""),
                       if (nchar(layoutValue           ) != 0) paste(", layout=",            layoutValue           , sep=""),
                       if (nchar(strip.valuesValue     ) != 0) paste(", strip.values=",      strip.valuesValue     , sep=""),
                       if (nchar(strip.left.valuesValue) != 0) paste(", strip.left.values=", strip.left.valuesValue, sep=""),
                       if (nchar(strip.parValue        ) != 0) paste(", strip.par=",         strip.parValue        , sep=""),
                       if (nchar(strip.left.parValue   ) != 0) paste(", strip.left.par=",    strip.left.parValue   , sep=""),
                       if (nchar(resize.heightValue    ) != 0) paste(", resize.height=",     resize.heightValue    , sep=""),
                       if (nchar(resize.widthValue     ) != 0) paste(", resize.width=",      resize.widthValue     , sep=""),
                       if (nchar(mainValue             ) != 0) paste(", main='",             mainValue             , "'", sep=""),
                       if (nchar(mainMiddleValue       ) != 0) paste(", main.middle=",       mainMiddleValue       , sep=""),
                        ")", sep="")

      doItAndPrint(command)
      activateMenus()
      tkfocus(CommanderWindow())
      if (version$os == "mingw32") justDoIt("bringToTop()")
    }
    OKCancelHelp(helpSubject="ResizeEtcDialog")
    tkgrid(tklabel(resizeFrame, text=gettextRcmdr("c.list:")),            c.listEntry           , sticky="w")
    tkgrid(tklabel(resizeFrame, text=gettextRcmdr("condlevelsName:")),    condlevelsNameEntry   , sticky="w")
    tkgrid(tklabel(resizeFrame, text=gettextRcmdr("x.same:")),            x.sameEntry           , sticky="w")
    tkgrid(tklabel(resizeFrame, text=gettextRcmdr("y.same:")),            y.sameEntry           , sticky="w")
    tkgrid(tklabel(resizeFrame, text=gettextRcmdr("layout:")),            layoutEntry           , sticky="w")
    tkgrid(tklabel(resizeFrame, text=gettextRcmdr("strip.values:")),      strip.valuesEntry     , sticky="w")
    tkgrid(tklabel(resizeFrame, text=gettextRcmdr("strip.left.values:")), strip.left.valuesEntry, sticky="w")
    tkgrid(tklabel(resizeFrame, text=gettextRcmdr("strip.par:")),         strip.parEntry        , sticky="w")
    tkgrid(tklabel(resizeFrame, text=gettextRcmdr("strip.left.par:")),    strip.left.parEntry   , sticky="w")
    tkgrid(tklabel(resizeFrame, text=gettextRcmdr("resize.height:")),     resize.heightEntry    , sticky="w")
    tkgrid(tklabel(resizeFrame, text=gettextRcmdr("resize.width:")),      resize.widthEntry     , sticky="w")
    tkgrid(tklabel(resizeFrame, text=gettextRcmdr("main:")),              mainEntry             , sticky="w")
    tkgrid(tklabel(resizeFrame, text=gettextRcmdr("mainMiddle:")),        mainMiddleEntry       , sticky="w")
    tkgrid(resizeFrame, sticky="w")


    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=13, columns=2)
    }


listAllTrellisObjects <- function (envir = .GlobalEnv, ...) {
    objects <- ls(envir = envir, ...)
    if (length(objects) == 0)
        return(NULL)
    objects[sapply(objects, function(.x) {
        "trellis" %in% class(get(.x, envir = envir))
    })]
}

## source("c:/HOME/rmh/HH-R.package/RcmdrPlugin.HH2/R/ResizeEtcDialog.R")
