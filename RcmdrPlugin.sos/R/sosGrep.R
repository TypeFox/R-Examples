### code adapted from groupsBox() in Rcmdr utilities
### file created 04 Dec 2009 

sosGrep <- defmacro(recall=NULL, bLabel=gettextRcmdr("Filter by"), eLabel=gettextRcmdr("in"), 
    initialLabel=gettextRcmdr("Filter by"), expr={
    env <- environment()
    .grep <- ""
    .column <- NULL
    .ignCase <- NULL
    .grepLabel <- tclVar(paste(initialLabel, "...", sep=""))
    onGrep <- function(){
        initializeDialog(subdialog, title=gettextRcmdr("Filter the search results"))
        subGrepFrame <- tkframe(subdialog)
        grepInput <- tclVar("")
        grepField <- tkentry(subGrepFrame, width="15", textvariable=grepInput)
        radioButtons(window=subGrepFrame, name="column", 
          buttons=c("package", "function", "description"),
          values=c("Package", "Function", "Description"), 
          initialValue=..values[2],
          labels=gettextRcmdr(c("Package", "Function", "Description"))
          , title="", right.buttons=FALSE
          )
        ignCaseVariable <- tclVar("1")
        ignCaseCheckBox <- tkcheckbutton(subGrepFrame, variable=ignCaseVariable)
        onOKsub <- function() {
            grep <- paste(tclvalue(grepInput))
            column <- tclvalue(columnVariable)
            ignCase <- paste(tclvalue(ignCaseVariable))
            if (grep == ""){
                assign(".grep", "", envir=env)
                tclvalue(.grepLabel) <- paste(initialLabel, "...", sep="")
                tkconfigure(grepButton, foreground="black")
                if (GrabFocus()) tkgrab.release(subdialog)
                tkdestroy(subdialog)
                tkwm.deiconify(top)
                if (GrabFocus()) tkgrab.set(top)
                tkfocus(top)
                tkwait.window(top)
                return()
                }
            if (column == "Package"){
                column <- paste(", column=", "'", paste(tclvalue(columnVariable)), "'", sep="")
                columnLab <- paste("'", paste(tclvalue(columnVariable)), "' name", sep="")
            } else if (column == "Function"){
                column <- paste("", sep="")
                columnLab <- paste("'", paste(tclvalue(columnVariable)), "' name", sep="")
            } else if (column == "Description"){
                column <- paste(", column=", "'", paste(tclvalue(columnVariable)), "'", sep="")
                columnLab <- paste("function '", paste(tclvalue(columnVariable)), "'", sep="")
            }
            if (ignCase == "1"){
                ignCase <- paste(', ignore.case=TRUE', sep="")
            } else {
                ignCase <- paste("", sep="")
            }
            assign(".grep", grep, envir=env)
            assign(".column", column, envir=env)
            assign(".ignCase", ignCase, envir=env)
            tclvalue(.grepLabel) <- paste(bLabel, " '", grep, "' ", 
              eLabel, " ", columnLab, sep="")
            tkconfigure(grepButton, foreground="blue")
            if (GrabFocus()) tkgrab.release(subdialog)
            tkdestroy(subdialog)
            tkwm.deiconify(top)
            if (GrabFocus()) tkgrab.set(top)
            tkfocus(top)
            tkwait.window(top)
        }
        subOKCancelHelp(helpSubject="grepFn")
        tkgrid(tklabel(subGrepFrame, text=gettextRcmdr("Filter the search"), fg="blue"), 
          tklabel(subGrepFrame, text=gettextRcmdr("results..."), fg="blue"), sticky="w")
        tkgrid(tklabel(subGrepFrame, text=gettextRcmdr("Match:")), grepField, sticky="w")
        tkgrid(tklabel(subGrepFrame, text=gettextRcmdr("Column:")), columnFrame, sticky="w")
        tkgrid(tklabel(subGrepFrame, text=gettextRcmdr("Ignore case")), ignCaseCheckBox, sticky="w")
        tkgrid(subGrepFrame, sticky="w")
        tkgrid(subButtonsFrame, sticky="w")
        dialogSuffix(subdialog, onOK=onOKsub, rows=2, columns=1, focus=subdialog)
    }
#    grepFrame <- tkframe(top)
    grepButton <- tkbutton(searchFrame, textvariable=.grepLabel, command=onGrep, borderwidth=3)
#    tkgrid(tklabel(searchFrame, text=gettextRcmdr(" ")), grepButton, sticky="w")
})
