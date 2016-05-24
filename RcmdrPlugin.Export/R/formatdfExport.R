## adapted from Rcmdr utilities.R

##FIXME CRAN checks: 
##listMatrixObjects: possible error in ls(envir = envir, ...): ... used
##in a situation where it does not exist
listMatrixObjects <- function(envir=.GlobalEnv, ...){
    objects <- ls(envir=envir, ...)
    if (length(objects) == 0) NULL
    else objects[sapply(objects,
      function(.x) "matrix" == (class(get(.x, envir=envir))[1]))]
}

## adapted from Rcmdr
formatdfExport <- function(){
    dataFrames <- listDataSets()
    matrices <- listMatrixObjects()
    framesMatrices <- sort(append(dataFrames, matrices))
    if (length(framesMatrices) == 0){
        Message(message=gettextRcmdr("There are no available data frames or matrices."),
          type="error")
        tkfocus(CommanderWindow())
        return()
    }
    initializeDialog(title=gettextRcmdr("Format for LaTeX or HTML export"))
    varsFrame <- tkframe(top)
    framesMatricesBox <- variableListBox(varsFrame, framesMatrices, 
      title=gettextRcmdr("Select data frame or matrix (pick one)"))
    optionsFrame <- tkframe(top)
    decInput <- tclVar("2")
    decField <- tkentry(optionsFrame, width="15", textvariable=decInput)
    cdecInput <- tclVar("")
    cdecField <- tkentry(optionsFrame, width="15", textvariable=cdecInput)
    dollarVariable <- tclVar("1")
    dollarCheckBox <- tkcheckbutton(optionsFrame, variable=dollarVariable)
    na.blankVariable <- tclVar("0")
    na.blankCheckBox <- tkcheckbutton(optionsFrame, variable=na.blankVariable)
    na.dotVariable <- tclVar("0")
    na.dotCheckBox <- tkcheckbutton(optionsFrame, variable=na.dotVariable)
    blank.dotVariable <- tclVar("0")
    blank.dotCheckBox <- tkcheckbutton(optionsFrame, variable=blank.dotVariable)
    col.justInput <- tclVar("")
    col.justField <- tkentry(optionsFrame, width="15", textvariable=col.justInput)
    cdotVariable <- tclVar("0")
    cdotCheckBox <- tkcheckbutton(optionsFrame, variable=cdotVariable)
    scientificInput <- tclVar("c(-4,4)")
    scientificField <- tkentry(optionsFrame, width="15", textvariable=scientificInput)
    math.rowVariable <- tclVar("0")
    math.rowCheckBox <- tkcheckbutton(optionsFrame, variable=math.rowVariable)
    math.colVariable <- tclVar("0")
    math.colCheckBox <- tkcheckbutton(optionsFrame, variable=math.colVariable)
    slashVariable <- tclVar("0")
    slashCheckBox <- tkcheckbutton(optionsFrame, variable=slashVariable)
    onOK <- function(){
        objectName <- getSelection(framesMatricesBox)
        dec <- paste(tclvalue(decInput))
        cdec <- paste(tclvalue(cdecInput))
        dollar <- paste(tclvalue(dollarVariable))
        na.blank <- paste(tclvalue(na.blankVariable))
        na.dot <- paste(tclvalue(na.dotVariable))
        blank.dot <- paste(tclvalue(blank.dotVariable))
        col.just <- paste(tclvalue(col.justInput))
        cdot <- paste(tclvalue(cdotVariable))
        scientific <- paste(tclvalue(scientificInput))
        math.row <- paste(tclvalue(math.rowVariable))
        math.col <- paste(tclvalue(math.colVariable))
        slash <- paste(tclvalue(slashVariable))
        closeDialog()
        if (length(objectName) == 0) {
            tkfocus(CommanderWindow())
            Message(message=gettextRcmdr("Please select an object."),
              type="error")
            return(formatdfExport())
        }
        if (dec != ""){
            dec <- paste(", dec=", paste(tclvalue(decInput)), sep="")
        }
        if (cdec != ""){
            cdec <- paste(", cdec=", paste(tclvalue(cdecInput)), sep="")
        }
        if (dollar == "0"){
            dollar <- paste(", numeric.dollar=FALSE", sep="")
        } else {
            dollar <- paste("", sep="")
        }

        if (na.blank == "1"){
            na.blank <- paste(", na.blank=TRUE", sep="")
        } else {
            na.blank <- paste("", sep="")
        }
        if (na.dot == "1"){
            na.dot <- paste(", na.dot=TRUE", sep="")
        } else {
            na.dot <- paste("", sep="")
        }
        if (blank.dot == "1"){
            blank.dot <- paste(", blank.dot=TRUE", sep="")
        } else {
            blank.dot <- paste("", sep="")
        }
        if (col.just != ""){
            col.just <- paste(", col.just=", paste(tclvalue(col.justInput)), sep="")
        }
        if (cdot == "1"){
            cdot <- paste(", cdot=TRUE", sep="")
        } else {
            cdot <- paste("", sep="")
        }
        if(scientific == "c(-4,4)"){
            scientific <- paste("", sep="")
        } else if(scientific != ""){
            scientific <- paste(", scientific=", paste(tclvalue(scientificInput)), sep="")
        }
        if (math.row == "1"){
            math.row <- paste(", math.row.names=TRUE", sep="")
        } else {
            math.row <- paste("", sep="")
        }
        if (math.col == "1"){
            math.col <- paste(", math.col.names=TRUE", sep="")
        } else {
            math.col <- paste("", sep="")
        }
        if (slash == "1"){
            slash <- paste(", double.slash=TRUE", sep="")
        } else {
            slash <- paste("", sep="")
        }
        doItAndPrint(paste("format.df(", objectName, dec, cdec, dollar, 
          na.blank, na.dot, blank.dot, col.just, cdot, scientific, 
          math.row, math.col, slash, ")", sep=""))
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="format.df")
    tkgrid(getFrame(framesMatricesBox), sticky="nw")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Arguments"), fg="blue"))
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Decimals:")), decField, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Col. decimals:")), cdecField, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Numeric dollar")), dollarCheckBox, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Blank NA")), na.blankCheckBox, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Dot NA")), na.dotCheckBox, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Dot blank")), blank.dotCheckBox, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Col. justify:")), col.justField, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Centered dot")), cdotCheckBox, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Scientific notation:")), scientificField, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Math. row-names")), math.rowCheckBox, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Math. col.-names")), math.colCheckBox, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Escape backslashes")), slashCheckBox, sticky="w")
    tkgrid(varsFrame, tklabel(top, text=" "), optionsFrame, sticky="nw")
    tkgrid(buttonsFrame, columnspan=3, sticky="w")
    dialogSuffix(rows=2, columns=2)
}
