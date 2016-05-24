# This file is part of the RcmdrPlugin.Export package.
# The current code is based on code taken from Rcmdr.
# file created: 02 Feb 2008
# last modified: 07 Oct 2008

xtableExport <- function(){
    .tmpObject <- popOutput(keep=TRUE)
    objectClass <- class(.tmpObject)
    objectClass <- objectClass[[1]]
    if (is.null(objectClass) == "TRUE" | objectClass == "logical"){
        .tmpObject <- popOutput()
        rm(.tmpObject)
        Message(message=paste("xtable() cannot export objects of class '", objectClass,
          "'", sep=""), type="note")
        Message(message=paste("the stack is probably empty", sep=""), type="warning")
        return()
        }
    else if (objectClass == "function" | objectClass == "help_files_with_topic" |
      objectClass == "packageIQR" | objectClass == "trellis" | objectClass == "xtable" |
      objectClass == "latex" | objectClass == "htest" | objectClass == "stem.leaf" |
      objectClass == "multinom" | objectClass == "try-error"){
        .tmpObject <- popOutput()
        rm(.tmpObject)
        Message(message=paste("xtable() cannot export objects of class '", objectClass,
          "'", sep=""), type="note")
        return(xtableExport())
        }
    else {
        initializeDialog(title=gettextRcmdr("Export objects using xtable()"))
        }
    dataFrame <- tkframe(top)
    xBox <- variableListBox(dataFrame, paste(objectClass), title=gettextRcmdr("Object class"))
    radioButtons(window=dataFrame, name="type", buttons=c("latex", "html"),
      values=c("LaTeX", "HTML"), labels=gettextRcmdr(c("LaTeX", "HTML")),
      title=gettextRcmdr("Export format"))
    optionsFrame <- tkframe(top)
    captionInput <- tclVar("")
    captionField <- tkentry(optionsFrame, width="15", textvariable=captionInput)
    radioButtons(window=optionsFrame, name="caption.place", 
      buttons=c("top", "bottom"),
      values=c("top", "bottom"), initialValue="bottom", 
      labels=gettextRcmdr(c("Top", "Bottom")), 
      title=gettextRcmdr("Caption place."))
    labelInput <- tclVar("")
    labelField <- tkentry(optionsFrame, width="15", textvariable=labelInput)
    alignInput <- tclVar("")
    alignField <- tkentry(optionsFrame, width="15", textvariable=alignInput)
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
        }
    else if (is.na(objectCol)) {
        digitsVector <- paste("2")
        }
    else digitsVector <- paste("c(", paste(rep("2,", objectCol), collapse=""), "2)", sep="")
    digitsInput <- tclVar(paste(digitsVector))
    digitsField <- tkentry(optionsFrame, width="15", textvariable=digitsInput)
    displayInput <- tclVar("")
    displayField <- tkentry(optionsFrame, width="15", textvariable=displayInput)
    additionalFrame <- tkframe(top)
    sizeInput <- tclVar("normal")
    sizeField <- tkentry(additionalFrame, width="15", textvariable=sizeInput)
    naInput <- tclVar("")
    naField <- tkentry(additionalFrame, width="15", textvariable=naInput)
    fileInput <- tclVar("")
    fileField <- tkentry(additionalFrame, width="15", textvariable=fileInput)
    appendVariable <- tclVar("1")
    appendCheckBox <- tkcheckbutton(additionalFrame, variable=appendVariable)
    floatVariable <- tclVar("1")
    floatCheckBox <- tkcheckbutton(additionalFrame, variable=floatVariable)
    radioButtons(window=additionalFrame, name="tab.env", 
      buttons=c("tabular", "longtable"),
      values=c("tabular", "longtable"), #initialValue="tabular", 
      labels=gettextRcmdr(c("tabular", "longtable")), 
      title=gettextRcmdr("Tabular env."))
    onOK <- function(){
        type <- tclvalue(typeVariable)
        caption <- paste(tclvalue(captionInput))
        caption.place <- tclvalue(caption.placeVariable)
        label <- paste(tclvalue(labelInput))
        align <- paste(tclvalue(alignInput))
        digits <- paste(tclvalue(digitsInput))
        if (digits == "2"){
            digits <- paste("")
            }
        else if (digits == digitsVector){
            digits <- paste("")
            }
        display <- paste(tclvalue(displayInput))
        size <- paste(tclvalue(sizeInput))
        na <- paste(tclvalue(naInput))
        file <- paste(tclvalue(fileInput))
        append <- paste(tclvalue(appendVariable))
        float <- paste(tclvalue(floatVariable))
        tab.env <- tclvalue(tab.envVariable)
        closeDialog()
        if (caption != ""){
            caption <- paste(", caption=", '"', paste(tclvalue(captionInput)), '"', sep="")
            if (caption.place == "top"){
                caption.place <- paste(", caption.placement=", '"', paste(tclvalue(caption.placeVariable)), '"', sep="")
            } else 
            caption.place <- paste("", sep="")
        } else 
        caption.place <- paste("", sep="")
        if (label != ""){
            label <- paste(", label=", '"', paste(tclvalue(labelInput)), '"', sep="")
            }
        if (align != ""){
            align <- paste(", align=", '"', paste(tclvalue(alignInput)), '"', sep="")
            }
        if (digits != ""){
            digits <- paste(", digits=", paste(tclvalue(digitsInput)), sep="")
            }
        if (display != ""){
            display <- paste(", display=", '"', paste(tclvalue(displayInput)), '"', sep="")
            }
        if (size != ""){
            if (size == "normal"){
                size <- paste("", sep="")
                }
            else size <- paste(", size=", '"', paste(tclvalue(sizeInput)), '"', sep="")
            }
        if (na != ""){
            na <- paste(", NA.string=", '"', paste(tclvalue(naInput)), '"', sep="")
            }
        if (file != ""){
            if (type == "LaTeX"){
                file <- paste(", file=", '"', paste(tclvalue(fileInput)), '.tex"', sep="")
                }
            else if (type == "HTML"){
                file <- paste(", file=", '"', paste(tclvalue(fileInput)), '.html"', sep="")
                }
            if (append == "1"){
                append <- paste(", append=TRUE", sep="")
                }
            else if (append == "0"){
                append <- paste("", sep="")
                }
        }
        else if (file == ""){
           append <- paste("", sep="")
           }
        if (float == "1"){
            float <- paste("", sep="")
        } else 
        if (float == "0"){
            float <- paste(", floating=FALSE", sep="")
        }
        if (tab.env == "longtable"){
            tab.env <- paste(", tabular.environment=", '"', paste(tclvalue(tab.envVariable)), '"', sep="")
        } else 
        tab.env <- paste("", sep="")
        objectCommandName <- NULL
        commandRepeat <- 1
        objectName <- paste(".", objectClass, sep="")
        if (objectClass == "numSummary"){
            objectCommandName <- paste("as.table(", objectName, "$table)", sep="")
        } else if (objectClass == "summary.multinom"){
            objectCommandName[1] <- paste("as.data.frame(", objectName, "$coefficients)", sep="")
            objectCommandName[2] <- paste("as.data.frame(", objectName, "$standard.errors)", sep="")
            commandRepeat <- 2
        } else if (objectClass == "polr"){
            objectCommandName[1] <- paste("as.data.frame(", objectName, "$coefficients)", sep="")
            objectCommandName[2] <- paste("as.data.frame(", objectName, "$zeta)", sep="")
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
        } else if (objectClass == "rcorr"){
            commandRepeat <- 3
            objectCommandList <- c("$r", "$n", "$P")
            for (i in 1:commandRepeat){
                objectCommandName[i] <- paste(objectName, objectCommandList[i], sep="")
            }
        } else if (objectClass == "rcorr.adjust"){
            commandRepeat <- 4
            objectCommandList <- c("$R$r", "$R$n", "$R$P", "$P")
            for (i in 1:commandRepeat){
                objectCommandName[i] <- paste(objectName, objectCommandList[i], sep="")
            }
        } else if (objectClass == "by" & is.list(.tmpObject)==TRUE){
            commandRepeat <- length(.tmpObject)
            for (i in 1:commandRepeat){
                objectCommandName[i] <- paste("as.matrix(", objectName,
                  "[[", i, "]])", sep="")
            }
        } else if (objectClass == "summary.princomp"){
            objectCommandName <- paste("as.matrix(", objectName, "$sdev)", sep="")
        } else {
            objectName <- paste(".object", sep="")
            objectCommandName <- paste(objectName)
        }
        assign(objectName, .tmpObject)
        cmds <- character(1)
        cmds[1] <- paste0("local({\n  ", 
                          "## retrieve the last printed object\n  ",
                          objectName, " <- popOutput()")
        # justDoIt(paste(objectName, " <- .tmpObject", sep=""))
        # logger(paste(objectName, " <- popOutput()   ## retrieve the last printed object", sep=""))
        ##always use print for length(commandRepeat)>1
        if(type == "LaTeX"){
            if(size!="" | na!="" | file!="" | caption.place!="" | 
              float!="" | tab.env!="" | commandRepeat>1){
                usePrint <- paste("print(", sep="")
                printType <- paste(", type=", '"', "latex", '"', ")", sep="")
            } else {
                usePrint <- paste("", sep="")
                printType <- paste("", sep="")
            }
        }
        else if (type == "HTML"){
            usePrint <- paste("print(", sep="")
            printType <- paste(", type=", '"', "html", '"', ")", sep="")
        }
        run.command <- character(commandRepeat)
        for (i in 1:commandRepeat){
            run.command[i] <- paste0("  ", usePrint, "xtable(", objectCommandName[i], 
              caption, label, align, digits, display, ")", caption.place, 
              size, na, file, append, float, tab.env, printType)
        }
        commands <- paste(c(cmds[1], run.command), collapse="\n")
        doItAndPrint(paste(commands, "\n})", sep=""))
        
        # justDoIt(paste('rm(list=c(".tmpObject", "',  objectName, '"))', sep=""))
        #    logger(paste("remove(", objectName, ")", sep=""))
        tkdestroy(top)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="xtable")
    tkgrid(getFrame(xBox), sticky="nw")
    tkgrid(typeFrame, sticky="sw")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Arguments"), fg="blue"))
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Caption:")), captionField, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Label:")), labelField, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Align:")), alignField, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Digits:")), digitsField, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Display:")), displayField, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr(" ")), sticky="w")
    tkgrid(caption.placeFrame, sticky="sw")
    tkgrid(tklabel(additionalFrame, text=gettextRcmdr("Printing"), fg="blue"))
    tkgrid(tklabel(additionalFrame, text=gettextRcmdr("Size:")), sizeField, sticky="w")
    tkgrid(tklabel(additionalFrame, text=gettextRcmdr("NA string:")), naField, sticky="w")
    tkgrid(tklabel(additionalFrame, text=gettextRcmdr("File:")), fileField, sticky="w")
    tkgrid(tklabel(additionalFrame, text=gettextRcmdr("Append")), appendCheckBox, sticky="w")
    tkgrid(tklabel(additionalFrame, text=gettextRcmdr("Floating")), floatCheckBox, sticky="w")
    tkgrid(tklabel(additionalFrame, text=gettextRcmdr(" ")), sticky="w")
    tkgrid(tab.envFrame, sticky="sw")
    tkgrid(dataFrame, tklabel(top, text=" "), optionsFrame, tklabel(top, text=" "),
      additionalFrame, sticky="nw")
    tkgrid(buttonsFrame, columnspan=3, sticky="w")
    dialogSuffix(rows=3, columns=3)
}
