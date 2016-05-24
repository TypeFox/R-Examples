readDataSetFacto <-
function() {
    initializeDialog(title=.Facto_gettext("Read Text Data From File, Clipboard, or URL"))
    optionsFrame <- tkframe(top)
    dsname <- tclVar(.Facto_gettext("Dataset"))
    entryDsname <- ttkentry(optionsFrame, width="20", textvariable=dsname)
    radioButtons(optionsFrame, "location", buttons=c("local", "clipboard", "url"),
        labels=.Facto_gettext(c("Local file system", "Clipboard", "Internet URL")), title=.Facto_gettext("Location of Data File"))
    headerVariable <- tclVar("1")
    nameRows <- tclVar("0")
    headerCheckBox <- tkcheckbutton(optionsFrame, variable=headerVariable)
    rowCheckBox <- tkcheckbutton(optionsFrame, variable=nameRows)
 ##   clipboardVariable <- tclVar("0")
 ##   clipboardCheckBox <- tkcheckbutton(optionsFrame, variable=clipboardVariable)
    radioButtons(optionsFrame, "delimiter", buttons=c("whitespace", "commas", "tabs"),
        labels=.Facto_gettext(c("White space", "Commas", "Tabs")), title=.Facto_gettext("Field Separator"))
    otherButton <- ttkradiobutton(delimiterFrame, variable=delimiterVariable, value="other")
    otherVariable <- tclVar("")
    otherEntry <- ttkentry(delimiterFrame, width="4", textvariable=otherVariable)
    radioButtons(optionsFrame, "decimal", buttons=c("period", "comma"),
        labels=.Facto_gettext(c("Period [.]", "Comma [,]")), title=.Facto_gettext("Decimal-Point Character"))
    missingVariable <- tclVar("NA")
    missingEntry <- ttkentry(optionsFrame, width="8", textvariable=missingVariable)
    onOK <- function(){
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == ""){
            errorCondition(recall=readDataSet,
                message=.Facto_gettext("You must enter a name for the data set."))
                return()
                }
        if (!is.valid.name(dsnameValue)){
            errorCondition(recall=readDataSet,
                message=paste('"', dsnameValue, '" ', .Facto_gettext("is not a valid name."), sep=""))
            return()
            }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, .Facto_gettext("Data set")))){
                readDataSet()
                return()
                }
            }
##        clip <- tclvalue(clipboardVariable) == "1"
        location <- tclvalue(locationVariable)
        file <- if (location == "clipboard") "clipboard"
            else if (location == "local") tclvalue(tkgetOpenFile(filetypes=
                .Facto_gettext('{"Text Files" {".txt" ".TXT" ".dat" ".DAT" ".csv" ".CSV"}} {"All Files" {"*"}}')))
            else {
                initializeDialog(subdialog, title=.Facto_gettext("Internet URL"))
                onOKsub <- function(){
                    closeDialog(subdialog)
                    }
                urlFrame <- tkframe(subdialog)
                urlVar <- tclVar("")
                url <- ttkentry(urlFrame, font=getRcmdr("logFont"), width="30", textvariable=urlVar)
                urlXscroll <- ttkscrollbar(urlFrame,
                        orient="horizontal", command=function(...) tkxview(url, ...))
                tkconfigure(url, xscrollcommand=function(...) tkset(urlXscroll, ...))
                subOKCancelHelp()
                tkgrid(url, sticky="w")
                tkgrid(urlXscroll, sticky="ew")
                tkgrid(urlFrame, sticky="nw")
                tkgrid(subButtonsFrame, sticky="w")
                dialogSuffix(subdialog, rows=2, columns=1, focus=url, onOK=onOKsub)
                tclvalue(urlVar)
                }
        if (file == "") {
            if (getRcmdr("grab.focus")) tkgrab.release(top)
            tkdestroy(top)
            return()
            }
        head <- tclvalue(headerVariable) == "1"
        row <- tclvalue(nameRows) == "1"
        delimiter <- tclvalue(delimiterVariable)
        del <- if (delimiter == "whitespace") ""
            else if (delimiter == "commas") ","
            else if (delimiter == "tabs") "\\t"
            else tclvalue(otherVariable)
        miss <- tclvalue(missingVariable)
        dec <- if (tclvalue(decimalVariable) == "period") "." else ","
        if (row) command <- paste('read.table("', file,'", header=', head,
            ', sep="', del, '", na.strings="', miss, '", dec="', dec, '", row.names=1, strip.white=TRUE)', sep="")
        else command <- paste('read.table("', file,'", header=', head,
            ', sep="', del, '", na.strings="', miss, '", dec="', dec, '", strip.white=TRUE)', sep="")
        logger(paste(dsnameValue, " <- ", command, sep=""))
        result <- justDoIt(command)
        if (class(result)[1] !=  "try-error"){
            gassign(dsnameValue, result)
            activeDataSet(dsnameValue)
            }
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="read.table")
    tkgrid(labelRcmdr(optionsFrame, text=.Facto_gettext("Enter name for data set:")), entryDsname, sticky="w")
    tkgrid(labelRcmdr(optionsFrame, text=.Facto_gettext("Variable names in file:")), headerCheckBox, sticky="w")
    tkgrid(labelRcmdr(optionsFrame, text=.Facto_gettext("Row names in the first columns:")), rowCheckBox, sticky="w")
##    tkgrid(labelRcmdr(optionsFrame, text=.Facto_gettext("Read data from clipboard:")), clipboardCheckBox, sticky="w")
    tkgrid(labelRcmdr(optionsFrame, text=.Facto_gettext("Missing data indicator:")), missingEntry, sticky="w")
    tkgrid(locationFrame, sticky="w")
    tkgrid(labelRcmdr(delimiterFrame, text=.Facto_gettext("Other")), otherButton,
        labelRcmdr(delimiterFrame, text=paste("  ",.Facto_gettext("Specify:"))), otherEntry, sticky="w")
    tkgrid(delimiterFrame, sticky="w", columnspan=2)
    tkgrid(decimalFrame, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=5, columns=1)
    }
