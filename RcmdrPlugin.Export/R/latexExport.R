# This file is part of the RcmdrPlugin.Export package.
# The current code is adapted from Rcmdr.
# file created: 02 Feb 2008
# last modified: 22 Oct 2009

latexExport <- function(){
    .tmpObject <- popOutput(keep=TRUE)
    objectClass <- class(.tmpObject)
    objectClass <- objectClass[[1]]
    objectName <- paste(".", objectClass, sep="")
    if (is.null(objectClass) == "TRUE" | objectClass == "logical"){
        .tmpObject <- popOutput()
        rm(.tmpObject)
        Message(message=paste("latex() cannot export objects of class '", objectClass,
          "'", sep=""), type="note")
        Message(message=paste("the stack is probably empty", sep=""), type="warning")
        return()
    } else if (objectClass == "function" | objectClass == "help_files_with_topic" |
        objectClass == "packageIQR" | objectClass == "trellis" | objectClass == "xtable" |
        objectClass == "latex" | objectClass == "stem.leaf" | objectClass == "multinom"){
        .tmpObject <- popOutput()
        rm(.tmpObject)
        Message(message=paste("latex() cannot export objects of class '", objectClass,
          "'", sep=""), type="note")
        return(latexExport())
    } else {
        initializeDialog(title=gettextRcmdr("Export objects using latex()"))
    }
    dataFrame <- tkframe(top)
    xBox <- variableListBox(dataFrame, paste(objectClass), title=gettextRcmdr("Object class"))
    optionsFrame <- tkframe(top)
    captionInput <- tclVar("")
    captionField <- tkentry(optionsFrame, width="15", textvariable=captionInput)
    radioButtons(window=optionsFrame, name="caption.loc", buttons=c("top", "bottom"),
      values=c("top", "bottom"), labels=gettextRcmdr(c("Top", "Bottom")),
      title=gettextRcmdr("Caption loc."))
    labelInput <- tclVar("")
    labelField <- tkentry(optionsFrame, width="15", textvariable=labelInput)
    if (objectClass == "numSummary") {
        objectCol <- ncol(.tmpObject$table)
    } else if (objectClass == "summary.lm") {
        objectCol <- ncol(.tmpObject$coefficients)
    } else if (objectClass == "rcorr"){
        objectCol <- ncol(.tmpObject$r)
    } else if (objectClass == "rcorr.adjust"){
        objectCol <- ncol(.tmpObject$R$r)
    } else if (objectClass == "confint.glht"){
        objectCol <- ncol(.tmpObject$confint)
    } else if (objectClass == "factanal"){
        objectCol <- ncol(.tmpObject$loadings)
    } else if (objectClass == "reliability"){
        objectCol <- ncol(.tmpObject$rel.matrix)
    } else if (objectClass == "summary.aov" & length(.tmpObject)==1) {
        objectCol <- ncol(.tmpObject[[1]])
    } else if (objectClass == "summary.princomp") {
        objectCol <- length(.tmpObject$sdev)
    } else {
        objectCol <- ncol(.tmpObject)
    }
    if (is.null(objectCol)) {
        digitsVector <- paste("2")
    } else if (is.na(objectCol)) {
        digitsVector <- paste("2")
    } else {
        digitsVector <- paste("c(", paste(rep("2,", (objectCol-1)), collapse=""), "2)", sep="")
    }
    digitsInput <- tclVar(paste(digitsVector))
    digitsField <- tkentry(optionsFrame, width="15", textvariable=digitsInput)
    additionalFrame <- tkframe(top)
    sizeInput <- tclVar("normal")
    sizeField <- tkentry(additionalFrame, width="15", textvariable=sizeInput)
    centerInput <- tclVar("center")
    centerField <- tkentry(additionalFrame, width="15", textvariable=centerInput)
    naVariable <- tclVar("1")
    naCheckBox <- tkcheckbutton(additionalFrame, variable=naVariable)
    fileInput <- tclVar("")
    fileField <- tkentry(additionalFrame, width="15", textvariable=fileInput)
    checkboxesFrame <- tkframe(top)
    appendVariable <- tclVar("1")
    appendCheckBox <- tkcheckbutton(checkboxesFrame, variable=appendVariable)
    visualiseVariable <- tclVar("0")
    visualiseCheckBox <- tkcheckbutton(checkboxesFrame, variable=visualiseVariable)
    longtableVariable <- tclVar("0")
    longtableCheckBox <- tkcheckbutton(checkboxesFrame, variable=longtableVariable)
    tab.envVariable <- tclVar("1")
    tab.envCheckBox <- tkcheckbutton(checkboxesFrame, variable=tab.envVariable)
    landscapeVariable <- tclVar("0")
    landscapeCheckBox <- tkcheckbutton(checkboxesFrame, variable=landscapeVariable)
    booktabsVariable <- tclVar("0")
    booktabsCheckBox <- tkcheckbutton(checkboxesFrame, variable=booktabsVariable)
    ctableVariable <- tclVar("0")
    ctableCheckBox <- tkcheckbutton(checkboxesFrame, variable=ctableVariable)
    vbarVariable <- tclVar("0")
    vbarCheckBox <- tkcheckbutton(checkboxesFrame, variable=vbarVariable)
    nomarginsVariable <- tclVar("1")
    nomarginsCheckBox <- tkcheckbutton(checkboxesFrame, variable=nomarginsVariable)
    onOK <- function(){
        caption <- paste(tclvalue(captionInput))
        label <- paste(tclvalue(labelInput))
        digits <- paste(tclvalue(digitsInput))
        caption.loc <- tclvalue(caption.locVariable)
        size <- paste(tclvalue(sizeInput))
        center <- paste(tclvalue(centerInput))
        na <- paste(tclvalue(naVariable))
        file <- paste(tclvalue(fileInput))
        append <- paste(tclvalue(appendVariable))
        visualise <- paste(tclvalue(visualiseVariable))
        longtable <- paste(tclvalue(longtableVariable))
        tab.env <- paste(tclvalue(tab.envVariable))
        landscape <- paste(tclvalue(landscapeVariable))
        booktabs <- paste(tclvalue(booktabsVariable))
        ctable <- paste(tclvalue(ctableVariable))
        closeDialog()
        vbar <- paste(tclvalue(vbarVariable))
        nomargins <- paste(tclvalue(nomarginsVariable))
        closeDialog()
        if (caption != ""){
            caption <- paste(", caption=", '"', paste(tclvalue(captionInput)), '"', sep="")
            }
        if (label != ""){
            label <- paste(", label=", '"', paste(tclvalue(labelInput)), '"', sep="")
            }
        if (digits != ""){
            digits <- paste(", cdec=", paste(tclvalue(digitsInput)), sep="")
            }
        if (caption != ""){
            if (caption.loc == "top"){
                caption.loc <- paste("", sep="")
            } else if (caption.loc == "bottom"){
                caption.loc <- paste(", caption.loc=", '"', paste(tclvalue(caption.locVariable)), '"', sep="")
            }
        } else {
            caption.loc <- paste("", sep="")
        }
        if (size != ""){
            if (size == "normal"){
                size <- paste("", sep="")
                }
            else size <- paste(", size=", '"', paste(tclvalue(sizeInput)), '"', sep="")
            }
        if (center != ""){
            if (center == "center"){
                center <- paste("", sep="")
                }
            else center <- paste(", center=", '"', paste(tclvalue(centerInput)), '"', sep="")
            }
        if (na == "1"){
            na <- paste("", sep="")
            }
        else if (na == "0"){
            na <- paste(', na.blank=FALSE', sep="")
            }
        if (file == ""){
            inObject <- paste("", sep="")
            if (visualise == "1"){
                secondTime <- TRUE
            } else {
                secondTime <- FALSE
            }
        } else if (file != ""){
            secondTime <- FALSE
            if (visualise == "1"){
                inObject <- paste("", sep="")
            } else {
                inObject <- paste(objectName, " <- ", sep="")
            }
        }
        if (file != ""){
            file <- paste(', file="', file, '.tex"', sep="")
            if (append == "1"){
                append <- paste(", append=TRUE", sep="")
            } else if (append == "0"){
                append <- paste("", sep="")
            }
        } else if (file == ""){
            file <- paste(', file=""', sep="")
            append <- paste("", sep="")
        }
        if (longtable == "1"){
            longtable <- paste(', longtable=TRUE', sep="")
        } else {
            longtable <- paste("", sep="")
        }
        if (tab.env == "0"){
            tab.env <- paste(', table.env=FALSE', sep="")
        } else {
            tab.env <- paste("", sep="")
        }
        if (landscape == "1"){
            landscape <- paste(', landscape=TRUE', sep="")
        } else if (landscape == "0"){
            landscape <- paste("", sep="")
        }
        if (booktabs == "1"){
            booktabs <- paste(', booktabs=TRUE', sep="")
            }
        else if (booktabs == "0"){
            booktabs <- paste("", sep="")
            }
        if (ctable == "1"){
            ctable <- paste(', ctable=TRUE', sep="")
            }
        else if (ctable == "0"){
            ctable <- paste("", sep="")
            }
        if (vbar == "1"){
            vbar <- paste(', vbar=TRUE', sep="")
            }
        else if (vbar == "0"){
            vbar <- paste("", sep="")
            }
        if (nomargins == "1"){
            nomargins <- paste("", sep="")
            }
        else if (nomargins == "0"){
            nomargins <- paste(', nomargins=FALSE', sep="")
            }
        functionName <- "latex"
        objectCommandName <- NULL
        commandRepeat <- 1
        if (objectClass == "numSummary"){
            objectCommandName <- paste(objectName, "$table", sep="")
        } else if (objectClass == "summary.lm") {
            objectCommandName <- paste(objectName, "$coefficients", sep="")
### use *[i], like in xtableExport()
        } else if (objectClass == "summary.multinom"){
            objectCommandName1 <- paste("as.data.frame(", objectName, "$coefficients)", sep="")
            objectCommandName2 <- paste("as.data.frame(", objectName, "$standard.errors)", sep="")
            objectCommandName <- c(objectCommandName1, objectCommandName2)
            commandRepeat <- 2
        } else if (objectClass == "polr"){
            objectCommandName1 <- paste("as.data.frame(", objectName, "$coefficients)", sep="")
            objectCommandName2 <- paste("as.data.frame(", objectName, "$zeta)", sep="")
            objectCommandName <- c(objectCommandName1, objectCommandName2)
            commandRepeat <- 2
        } else if (objectClass == "summary.polr"){
            objectCommandName <- paste(objectName, "$coefficients", sep="")
        } else if (objectClass == "reliability"){
            objectCommandName <- paste(objectName, "$rel.matrix", sep="")
        } else if (objectClass == "confint.glht"){
            objectCommandName <- paste(objectName, "$", "confint", sep="")
        } else if (objectClass == "factanal"){
            objectCommandName <- paste("as.table(", objectName, "$loadings)", sep="")
        } else if (objectClass == "outlier.test"){
            objectCommandName <- paste("as.matrix(", objectName, "$test)", sep="")
        } else if (objectClass == "array" | objectClass == "integer" |
          objectClass == "character" | objectClass == "numeric"){
            objectCommandName <- paste("as.data.frame(", objectName, ")", sep="")
        ###FIXME support for `rcorr' possibly buggy
        } else if (objectClass == "rcorr"){
            objectCommandName <- paste(objectName, sep="")
            functionName <- "latex.list"
        } else if (objectClass == "rcorr.adjust"){
            commandRepeat <- 4
            objectCommandList <- c("$R$r", "$R$n", "$R$P", "$P")
            for (i in 1:commandRepeat){
                objectCommandName[i] <- paste(objectName, objectCommandList[i], sep="")
            }
        } else if (objectClass == "by" & is.list(.tmpObject)==TRUE){
            commandRepeat <- length(.tmpObject)
#            objectCommandName <- NULL
            for (i in 1:commandRepeat){
                objectCommandName[i] <- paste("as.matrix(", objectName,
                  "[[", i, "]])", sep="")
            }
        } else if (objectClass == "table"){
            objectCommandName <- paste("as.matrix(", objectName, ")", sep="")
        } else if (objectClass == "summary.aov" & length(.tmpObject)==1) {
            objectCommandName <- paste(objectName, "[[1]]", sep="")
        } else if (objectClass == "summary.princomp"){
            objectCommandName <- paste(objectName, "$sdev", sep="")
        } else {
            objectName <- paste(".object", sep="")
            objectCommandName <- paste(objectName)
        }
        assign(objectName, .tmpObject)
        if (inObject != ""){
            inObject <- paste(objectName, " <- ", sep="")
        }
        
        cmds <- character(3)
        cmds[1] <- paste0("local({\n  ", 
                         "## retrieve the last printed object\n  ",
                         objectName, " <- popOutput()")
        #justDoIt(paste(objectName, " <- .tmpObject", sep=""))
        #eval(parse(text=paste(objectName, " <- .tmpObject", sep="")))
        #logger(paste(objectName, " <- popOutput()   ## retrieve the last printed object", sep=""))
        .matPercentage <- FALSE
        if (objectClass == "matrix"){
            eval(parse(text=paste('.matPercentage <- !(nrow(as.matrix(grep("%",
                colnames(', objectCommandName, "), fixed=TRUE))) == 0)",
                sep="")))
            }
        need.sanitize <- (objectClass == "numSummary" || .matPercentage == TRUE)
        if (need.sanitize){
            run.sanitize <- paste("  ## escape strings for LaTeX\n  ",
                "colnames(", objectCommandName,
                ") <- \n    latexTranslate(", "colnames(", objectCommandName,
                "))", sep="")
            cmds[2] <- run.sanitize
            #logger(run.sanitize)
            #eval(parse(text=run.sanitize))
        } else {
            cmds[2] <- ""
        }
        run.command <- character(commandRepeat)
        for (i in 1:commandRepeat){
            run.command[i] <- paste("  ", inObject, functionName, "(", objectCommandName[i],
              caption, caption.loc, label, digits, size, na, file, append,
              longtable, tab.env, landscape, booktabs, ctable, vbar, nomargins,
              center, ', title="")', sep="")
            #logger(run.command)
            #eval(parse(text=run.command))
        }
        if (secondTime){
            file <- paste0("")  ##default file runs DVI preview
            ##without print() preview works only for last call
            if(commandRepeat > 1){
                usePrint <- paste0("print(")
                endPrint <- paste0(")")
            } else {
                usePrint <- paste0("")
                endPrint <- paste0("")
            }
            run.preview <- character(commandRepeat)
            for (i in 1:commandRepeat){
                run.preview[i] <- paste0(if(i==1) "  ## DVI preview\n  " else "  ", 
                                        usePrint, inObject, functionName, "(", objectCommandName[i],
                                        caption, caption.loc, label, digits, size, na, file, append,
                                        longtable, tab.env, landscape, booktabs, ctable, vbar, 
                                        nomargins, center, ', title="")', endPrint)
                #logger(run.preview)
                #eval(parse(text=run.preview))
            }
        }
        commands <- paste(c(cmds[1], if(need.sanitize) run.sanitize, 
            run.command, if(secondTime) run.preview), collapse="\n")
        doItAndPrint(paste(commands, "\n})", sep=""))
        #eval(parse(text=paste('rm(list=c(', '"',  objectName, '"))', sep="")))
        #logger(paste("remove(", objectName, ")", sep=""))
        
        tkdestroy(top)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="latex")
    tkgrid(getFrame(xBox), sticky="nw")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Arguments"), fg="blue"))
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Caption:")), captionField, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Label:")), labelField, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Digits:")), digitsField, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr(" ")), sticky="w")
    tkgrid(caption.locFrame, sticky="sw")
    tkgrid(tklabel(additionalFrame, text=gettextRcmdr("Printing    "), fg="blue"))
    tkgrid(tklabel(additionalFrame, text=gettextRcmdr("Size:")), sizeField, sticky="w")
    tkgrid(tklabel(additionalFrame, text=gettextRcmdr("Center:")), centerField, sticky="w")
    tkgrid(tklabel(additionalFrame, text=gettextRcmdr("Blank NA")), naCheckBox, sticky="w")
    tkgrid(tklabel(additionalFrame, text=gettextRcmdr("File:")), fileField, sticky="w")
    tkgrid(tklabel(checkboxesFrame, text=gettextRcmdr("Append")), appendCheckBox, sticky="w")
    tkgrid(tklabel(checkboxesFrame, text=gettextRcmdr("Preview")), visualiseCheckBox, sticky="w")
    tkgrid(tklabel(checkboxesFrame, text=gettextRcmdr("Table env.")), tab.envCheckBox, sticky="w")
    tkgrid(tklabel(checkboxesFrame, text=gettextRcmdr("Use 'longtable'")), longtableCheckBox, sticky="w")
    tkgrid(tklabel(checkboxesFrame, text=gettextRcmdr("Landscape")), landscapeCheckBox, sticky="w")
    tkgrid(tklabel(checkboxesFrame, text=gettextRcmdr("Use 'booktabs'")), booktabsCheckBox, sticky="w")
    tkgrid(tklabel(checkboxesFrame, text=gettextRcmdr("Use 'ctable'")), ctableCheckBox, sticky="w")
    tkgrid(tklabel(checkboxesFrame, text=gettextRcmdr("Vertical bar")), vbarCheckBox, sticky="w")
    tkgrid(tklabel(checkboxesFrame, text=gettextRcmdr("No margins")), nomarginsCheckBox, sticky="w")
    tkgrid(dataFrame, tklabel(top, text=" "), additionalFrame, sticky="nw")
    tkgrid(optionsFrame, tklabel(top, text=" "), checkboxesFrame, sticky="nw")
    tkgrid(buttonsFrame, columnspan=3, sticky="w")
    dialogSuffix(rows=4, columns=3)
}
