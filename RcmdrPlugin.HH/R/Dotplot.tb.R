"DotplottbRcmdr" <-
function() {
    initializeDialog(title=gettextRcmdr("dotplot with tiebreakers"))
    factorBox <- variableListBox(top, Factors(), title=gettextRcmdr("Factors (pick zero or one)"), selectmode="single", initialSelection=-1)
    responseBox <- variableListBox(top, Numeric(), title=gettextRcmdr("~ Response Variable (pick one)"))
    groupBox <- variableListBox(top, Factors(), title=gettextRcmdr(" | Groups (pick zero or more)"), selectmode="multiple", initialSelection=-1)

    scalarsFrame <- tkframe(top)

    factor.jitterVar <- tclVar(".5")
    factor.jitterEntry <- tkentry(scalarsFrame, width="6", textvariable=factor.jitterVar)

    cexVar <- tclVar("1")
    cexEntry <- tkentry(scalarsFrame, width="6", textvariable=cexVar)

    pchVar <- tclVar("16")
    pchEntry <- tkentry(scalarsFrame, width="6", textvariable=pchVar)
    
    layxVar <- tclVar("")
    layxEntry <- tkentry(scalarsFrame, width="6", textvariable=layxVar)
    
    layyVar <- tclVar("")
    layyEntry <- tkentry(scalarsFrame, width="6", textvariable=layyVar)
    
    betxVar <- tclVar("0")
    betxEntry <- tkentry(scalarsFrame, width="6", textvariable=betxVar)
    
    betyVar <- tclVar("1")
    betyEntry <- tkentry(scalarsFrame, width="6", textvariable=betyVar)

    altxVar <- tclVar("1")
    altxEntry <- tkentry(scalarsFrame, width="6", textvariable=altxVar)
    
    altyVar <- tclVar("1")
    altyEntry <- tkentry(scalarsFrame, width="6", textvariable=altyVar)

    onOK <- function() {
        factor <- getSelection(factorBox)
        groups <- getSelection(groupBox)
        response <- getSelection(responseBox)

        factor.jitter  <- as.numeric(tclvalue(factor.jitterVar))
        cex            <- as.numeric(tclvalue(cexVar))
        old.warn <- options(warn=-1)
        pch            <- as.numeric(tclvalue(pchVar))
        if (is.na(pch)) pch <- paste("\'", tclvalue(pchVar), "\'", sep="")
        options(old.warn)
        layx           <- as.numeric(tclvalue(layxVar))
        layy           <- as.numeric(tclvalue(layyVar))
        betx           <- as.numeric(tclvalue(betxVar))
        bety           <- as.numeric(tclvalue(betyVar))
        altx           <- as.numeric(tclvalue(altxVar))
        alty           <- as.numeric(tclvalue(altyVar))
        
        closeDialog()
        if (0 == length(response)) {
            errorCondition(recall=DotplottbRcmdr,
                           message=gettextRcmdr("No response variable selected."))
            return()
            }
        .activeDataSet <- ActiveDataSet()
        groups.formula <- paste(groups, collapse=' + ')
        if (nchar(groups.formula) > 0) {
          groups.formula <- paste(" | ", groups.formula, sep="")
          freq.formula <- paste("with(get(.activeDataSet),",
                                "tapply(",
                                "get(response)",
                                ", data.frame(factor(get(response))",
                                if (length(factor)==1) ", get(factor)",
                                if (length(groups)>0) ", get(paste(groups, collapse=','))",
                                "), length))")
          freq <- eval(parse(text=freq.formula))
          max.freq <- max(freq, na.rm=TRUE)
          groups.formula <- paste(groups.formula,
                                  ", max.freq=", max.freq, sep="")
        }
        if (length(factor)==0) {
          n.response <- eval(parse(text='length(get(.activeDataSet)[[response]])'))
          factor <- paste("rep('', nrow(", .activeDataSet, "))", sep="")
        }
        dptb.command <- paste("dotplot(", factor, " ~ ",  response,
                              groups.formula,
                              ", data=", .activeDataSet, ",", sep="")

        if (is.na(layx) && !is.na(layy)) layx <- 1
        if (!is.na(layx) && is.na(layy)) layy <- 1
        if (is.na(layx) && is.na(layy))
          layout <- ""
        else
          layout <- paste(", layout=c(", layx, ",", layy, ")", sep="")

        between <- paste(", between=list(x=", betx, ",y=", bety, ")", sep="")
        alternating <- paste("scales=list(",
                             "x=list(alternating=", altx, "), ",
                             "y=list(alternating=", alty, "))", sep="")
        dptb.command <- paste(dptb.command,
                              " ", ## "\n        ",
                              "panel=panel.dotplot.tb",
                              ", factor=", factor.jitter,
                              ", cex=", cex,
                              ", pch=", pch,
                              layout,
                              between,
                              ", ", ## ",\n        ",
                              alternating,
                              ')',
                              sep="")
        doItAndPrint(dptb.command)
        activateMenus()
        tkfocus(CommanderWindow())
        }
    optionsFrame <- tkframe(top)
    buttonsFrame <- tkframe(top)
    OKCancelHelp(helpSubject="panel.dotplot.tb")
    tkgrid(getFrame(factorBox),
           getFrame(responseBox),
           getFrame(groupBox),
           sticky="nw")

    tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("Change defaults if needed:")))
    tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("factor.jitter:")), factor.jitterEntry, sticky="w")
    tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("cex:")), cexEntry, sticky="w")
    tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("pch:")), pchEntry, sticky="w")
    tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("layx:")), layxEntry, sticky="w")
    tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("layy:")), layyEntry, sticky="w")
    tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("betx:")), betxEntry, sticky="w")
    tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("bety:")), betyEntry, sticky="w")
    tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("altx:")), altxEntry, sticky="w")
    tkgrid(tklabel(scalarsFrame, text=gettextRcmdr("alty:")), altyEntry, sticky="w")
    tkgrid(scalarsFrame, sticky="w")
    
    tkgrid(optionsFrame, columnspan=2, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=3, columns=2)
  }
##source("~/HH-R.package/RcmdrPlugin.HH/R/Dotplot.tb.R")
##source("~/HH-R.package/HH/R/panel.dotplot.tb.R")
