# Copied from R-forge Dec 1, 2015 by Jessica Reese
# all credit goes to John Fox and those he lists as contributers

## Import/Export Data

# selecting an active data set (button in Rcmdr window)
# selectActiveDataSet <- function(){
#   dataSets <- listDataSets()
#   .activeDataSet <- ActiveDataSet()
#   if ((length(dataSets) == 1) && !is.null(.activeDataSet)) {
#     Message(message=gettextRcmdr("There is only one dataset in memory."),
#             type="warning")
#     tkfocus(CommanderWindow())
#     return()
#   }
#   if (length(dataSets) == 0){
#     Message(message=gettextRcmdr("There are no data sets from which to choose."),
#             type="error")
#     tkfocus(CommanderWindow())
#     return()
#   }
#   initializeDialog(title=gettextRcmdr("Select Data Set"))
#   dataSetsBox <- variableListBox(top, dataSets, title=gettextRcmdr("Data Sets (pick one)"),
#                                  initialSelection=if (is.null(.activeDataSet)) NULL else which(.activeDataSet == dataSets) - 1)
#   onOK <- function(){
#     activeDataSet(getSelection(dataSetsBox))
#     closeDialog()
#     tkfocus(CommanderWindow())
#   }
#   OKCancelHelp()
#   tkgrid(getFrame(dataSetsBox), sticky="nw")
#   tkgrid(buttonsFrame, sticky="w")
#   dialogSuffix(rows=2, columns=1)
# }
# 
# # load data set
# loadDataSet <- function() {
#   file <- tclvalue(tkgetOpenFile(filetypes=
#                                    gettextRcmdr('{"All Files" {"*"}} {"R Data Files" {".RData" ".rda" ".Rda" ".RDA"}}')))
#   if (file == "") return()
#   command <- paste('load("', file,'")', sep="")
#   dsname <- justDoIt(command)
#   logger(command)
#   if (class(dsname)[1] !=  "try-error") activeDataSet(dsname)
#   tkfocus(CommanderWindow())
# }
# 
# # import data set from file, clipboard, or url
# readDataSet <- function() {
#   initializeDialog(title=gettextRcmdr("Read Text Data From File, Clipboard, or URL"))
#   optionsFrame <- tkframe(top)
#   dsname <- tclVar(gettextRcmdr("Dataset"))
#   entryDsname <- ttkentry(optionsFrame, width="20", textvariable=dsname)
#   radioButtons(optionsFrame, "location", buttons=c("local", "clipboard", "url"), 
#                labels=gettextRcmdr(c("Local file system", "Clipboard", "Internet URL")), title=gettextRcmdr("Location of Data File"))
#   headerVariable <- tclVar("1")
#   headerCheckBox <- tkcheckbutton(optionsFrame, variable=headerVariable)
#   radioButtons(optionsFrame, "delimiter", buttons=c("whitespace", "commas", "tabs"),
#                labels=gettextRcmdr(c("White space", "Commas", "Tabs")), title=gettextRcmdr("Field Separator"))
#   otherButton <- ttkradiobutton(delimiterFrame, variable=delimiterVariable, value="other")
#   otherVariable <- tclVar("")
#   otherEntry <- ttkentry(delimiterFrame, width="4", textvariable=otherVariable)
#   radioButtons(optionsFrame, "decimal", buttons=c("period", "comma"),
#                labels=gettextRcmdr(c("Period [.]", "Comma [,]")), title=gettextRcmdr("Decimal-Point Character"))
#   missingVariable <- tclVar("NA")
#   missingEntry <- ttkentry(optionsFrame, width="8", textvariable=missingVariable)
#   onOK <- function(){
#     closeDialog()
#     dsnameValue <- trim.blanks(tclvalue(dsname))
#     if (dsnameValue == ""){
#       errorCondition(recall=readDataSet,
#                      message=gettextRcmdr("You must enter a name for the data set."))
#       return()
#     }
#     if (!is.valid.name(dsnameValue)){
#       errorCondition(recall=readDataSet,
#                      message=paste('"', dsnameValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
#       return()
#     }
#     if (is.element(dsnameValue, listDataSets())) {
#       if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
#         readDataSet()
#         return()
#       }
#     }
#     location <- tclvalue(locationVariable)
#     file <- if (location == "clipboard") "clipboard" 
#     else if (location == "local") tclvalue(tkgetOpenFile(filetypes=
#                                                            gettextRcmdr('{"All Files" {"*"}} {"Text Files" {".txt" ".TXT" ".dat" ".DAT" ".csv" ".CSV"}}')))
#     else {
#       initializeDialog(subdialog, title=gettextRcmdr("Internet URL"))
#       onOKsub <- function(){
#         closeDialog(subdialog)
#       }
#       urlFrame <- tkframe(subdialog)
#       urlVar <- tclVar("")
#       url <- ttkentry(urlFrame, font=getRcmdr("logFont"), width="30", textvariable=urlVar)
#       urlXscroll <- ttkscrollbar(urlFrame,
#                                  orient="horizontal", command=function(...) tkxview(url, ...))
#       tkconfigure(url, xscrollcommand=function(...) tkset(urlXscroll, ...))
#       subOKCancelHelp()
#       tkgrid(url, sticky="w")
#       tkgrid(urlXscroll, sticky="ew")
#       tkgrid(urlFrame, sticky="nw")
#       tkgrid(subButtonsFrame, sticky="w")
#       dialogSuffix(subdialog, rows=2, columns=1, focus=url, onOK=onOKsub)
#       tclvalue(urlVar)
#     }
#     if (file == "") {
#       if (getRcmdr("grab.focus")) tkgrab.release(top)
#       tkdestroy(top)
#       return()
#     }
#     head <- tclvalue(headerVariable) == "1"
#     delimiter <- tclvalue(delimiterVariable)
#     del <- if (delimiter == "whitespace") ""
#     else if (delimiter == "commas") ","
#     else if (delimiter == "tabs") "\\t"
#     else tclvalue(otherVariable)
#     miss <- tclvalue(missingVariable)
#     dec <- if (tclvalue(decimalVariable) == "period") "." else ","
#     command <- paste('read.table("', file,'", header=', head,
#                      ', sep="', del, '", na.strings="', miss, '", dec="', dec, '", strip.white=TRUE)', sep="")
#     logger(paste(dsnameValue, " <- ", command, sep=""))
#     result <- justDoIt(command)
#     if (class(result)[1] !=  "try-error"){
#       gassign(dsnameValue, result)
#       activeDataSet(dsnameValue)
#     }
#     tkfocus(CommanderWindow())
#   }
#   OKCancelHelp(helpSubject="read.table")
#   tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="w")
#   tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Variable names in file:")), headerCheckBox, sticky="w")
#   tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Missing data indicator:")), missingEntry, sticky="w")
#   tkgrid(locationFrame, sticky="w")
#   tkgrid(labelRcmdr(delimiterFrame, text=gettextRcmdr("Other")), otherButton,
#          labelRcmdr(delimiterFrame, text=gettextRcmdr("  Specify:")), otherEntry, sticky="w")
#   tkgrid(delimiterFrame, sticky="w", columnspan=2)
#   tkgrid(decimalFrame, sticky="w")
#   tkgrid(optionsFrame, sticky="w")
#   tkgrid(buttonsFrame, sticky="w")
#   dialogSuffix(rows=5, columns=1)
# }
# 
# # import data set from SPSS file
# importSPSS <- function() {
#   Library("foreign")
#   initializeDialog(title=gettextRcmdr("Import SPSS Data Set"))
#   dsname <- tclVar(gettextRcmdr("Dataset"))
#   entryDsname <- ttkentry(top, width="20", textvariable=dsname)
#   asFactor <- tclVar("1")
#   asFactorCheckBox <- tkcheckbutton(top, variable=asFactor)
#   toLower <- tclVar("1")
#   toLowerCheckBox <- tkcheckbutton(top, variable=toLower)
#   maxLevels <- tclVar("Inf")
#   entryMaxLevels <- ttkentry(top, width="5", textvariable=maxLevels)
#   onOK <- function(){
#     closeDialog()
#     dsnameValue <- trim.blanks(tclvalue(dsname))
#     if (dsnameValue == ""){
#       errorCondition(recall=importSPSS,
#                      message=gettextRcmdr("You must enter the name of a data set."))
#       return()
#     }
#     if (!is.valid.name(dsnameValue)){
#       errorCondition(recall=importSPSS,
#                      message=paste('"', dsnameValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
#       return()
#     }
#     if (is.element(dsnameValue, listDataSets())) {
#       if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
#         importSPSS()
#         return()
#       }
#     }
#     file <- tclvalue(tkgetOpenFile(
#       filetypes=gettextRcmdr('{"All Files" {"*"}} {"SPSS portable files" {".por" ".POR"}} {"SPSS save files" {".sav" ".SAV"}}')))
#     if (file == "") {
#       tkfocus(CommanderWindow())
#       return()
#     }
#     factor <- tclvalue(asFactor) == "1"
#     levels <- as.numeric(tclvalue(maxLevels))
#     command <- paste('read.spss("', file,'", use.value.labels=', factor,
#                      ", max.value.labels=", levels, ", to.data.frame=TRUE)", sep="")
#     logger(paste(dsnameValue, " <- ", command, sep=""))
#     result <- justDoIt(command)
#     if (class(result)[1] !=  "try-error"){
#       gassign(dsnameValue, result)
#       if (tclvalue(toLower) == "1") 
#         doItAndPrint(paste("colnames(", dsnameValue, ") <- tolower(colnames(",
#                            dsnameValue, "))", sep=""))
#       activeDataSet(dsnameValue)
#     }
#     tkfocus(CommanderWindow())
#   }
#   OKCancelHelp(helpSubject="read.spss")
#   tkgrid(labelRcmdr(top, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="w")
#   tkgrid(labelRcmdr(top, text=gettextRcmdr("Convert value labels\nto factor levels"), justify="left"),
#          asFactorCheckBox, sticky="w")
#   tkgrid(labelRcmdr(top, text=gettextRcmdr("Convert variable names\nto lower case"), justify="left"),
#          toLowerCheckBox, sticky="w")
#   tkgrid(labelRcmdr(top, text=gettextRcmdr("Maximum number\nof value labels\nfor factor conversion"), justify="left"),
#          entryMaxLevels, sticky="w")
#   tkgrid(buttonsFrame, columnspan="2", sticky="w")
#   tkgrid.configure(entryDsname, sticky="w")
#   tkgrid.configure(asFactorCheckBox, sticky="w")
#   tkgrid.configure(toLowerCheckBox, stick="w")
#   tkgrid.configure(entryMaxLevels, sticky="w")
#   dialogSuffix(rows=5, columns=2, focus=entryDsname)
# }
# 
# # import data set from SAS file
# importSAS <- function() {
#   # the following local function is adapted from ?chartr
#   capwords <- function(s) {
#     cap <- function(s) paste(toupper(substring(s,1,1)), {s <- substring(s, 2); tolower(s)}, sep = "", collapse = " " )
#     sapply(strsplit(s, split = " "), cap)
#   }
#   Library("foreign")
#   file <- tclvalue(tkgetOpenFile(
#     filetypes=gettextRcmdr('{"All Files" {"*"}} {"SAS xport files" {".xpt" ".XPT" ".xport" ".XPORT"}}')))
#   if (file == "") {
#     tkfocus(CommanderWindow())
#     return()
#   }
#   command <- paste('read.xport("', file,'")', sep="")
#   logger(paste(".Datasets <- ", command, sep=""))
#   result <- justDoIt(command)
#   if (class(result)[1] !=  "try-error"){
#     gassign(".Datasets", result)
#     if (is.data.frame(.Datasets)){
#       getdsname <- function(){
#         initializeDialog(title=gettextRcmdr("Data Set Name"))
#         dsname <- tclVar(gettextRcmdr("Dataset"))
#         entryDsname <- ttkentry(top, width="20", textvariable=dsname)
#         onOK <- function(){
#           closeDialog()
#           dsnameValue <- trim.blanks(tclvalue(dsname))
#           if (dsnameValue == ""){
#             errorCondition(recall=getdsname,
#                            message=gettextRcmdr("You must enter the name of a data set."))
#             return()
#           }
#           if (!is.valid.name(dsnameValue)){
#             errorCondition(recall=getdsname,
#                            message=paste('"', dsnameValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
#             return()
#           }
#           if (is.element(dsnameValue, listDataSets())) {
#             if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
#               getdsname()
#               return()
#             }
#           }
#           doItAndPrint(paste(dsnameValue, " <- .Datasets", sep=""))
#           logger("remove(.Datasets)")
#           remove(".Datasets", envir=.GlobalEnv)
#           activeDataSet(dsnameValue)
#         }
#         OKCancelHelp()
#         tkgrid(labelRcmdr(top, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="w")
#         tkgrid(buttonsFrame, columnspan="2", sticky="w")
#         tkgrid.configure(entryDsname, sticky="w")
#         dialogSuffix(rows=2, columns=2, focus=entryDsname)
#       }
#       getdsname()
#     }
#     else {
#       fmt <- grep("^FORMAT", names(.Datasets))
#       if (length(fmt) >= 1) gassign(".Datasets", .Datasets[-fmt])
#       if (length(.Datasets) == 1){
#         dsname <- capwords(names(.Datasets))
#         if (is.element(dsname, listDataSets())) {
#           if ("no" == tclvalue(checkReplace(dsname, gettextRcmdr("Data set")))){
#             importSAS()
#             return()
#           }
#         }
#         #                 assign(dsname, .Datasets[[1]], envir=.GlobalEnv)
#         #                 logger(paste(dsname, " <- .Datasets[[1]]", sep=""))
#         doItAndPrint(paste(dsname, " <- .Datasets[[1]]", sep=""))
#         doItAndPrint(paste("colnames(", dsname, ") <- ", "tolower(colnames(", 
#                            dsname, "))", sep=""))
#         logger("remove(.Datasets)")
#         remove(".Datasets", envir=.GlobalEnv)
#         activeDataSet(dsname)
#       }
#       else {
#         dsnames <- capwords(names(.Datasets))
#         datasets <- listDataSets()
#         initializeDialog(title=gettextRcmdr("Select Dataset"))
#         datasetsBox <- variableListBox(top, dsnames, 
#                                        title=gettextRcmdr("Datasets in file (pick one)"),
#                                        initialSelection=0)
#         onOK <- function() {
#           dsname <- getSelection(datasetsBox)
#           for (ds in 1:length(dsnames)){
#             if (is.element(dsnames[ds], datasets)) {
#               if ("no" == tclvalue(checkReplace(dsnames[ds], gettextRcmdr("Data set")))){
#                 next()
#               }
#             }
#             #                         assign(dsnames[ds], .Datasets[[ds]], envir=.GlobalEnv)
#             #                         logger(paste(dsnames[ds], " <- .Datasets[[", ds, "]]", sep=""))
#             doItAndPrint(paste(dsnames[ds], " <- .Datasets[[", ds, "]]", sep=""))
#             doItAndPrint(paste("colnames(", dsnames[ds], ") <- ", "tolower(colnames(", 
#                                dsnames[ds], "))", sep=""))
#           }
#           logger("remove(.Datasets)")
#           remove(".Datasets", envir=.GlobalEnv)
#           activeDataSet(dsname)
#           closeDialog()
#           tkfocus(CommanderWindow())
#         }
#         OKCancelHelp(helpSubject="read.xport")
#         tkgrid(getFrame(datasetsBox), sticky="w")
#         tkgrid(buttonsFrame, sticky="w")
#         dialogSuffix(rows=2, columns=1)
#       }
#     }
#   }
#   tkfocus(CommanderWindow())
# }
# 
# # import data set from Minitab file
# importMinitab <- function() {
#   Library("foreign")
#   initializeDialog(title=gettextRcmdr("Import Minitab Data Set"))
#   dsname <- tclVar(gettextRcmdr("Dataset"))
#   entryDsname <- ttkentry(top, width="20", textvariable=dsname)
#   onOK <- function(){
#     closeDialog()
#     dsnameValue <- trim.blanks(tclvalue(dsname))
#     if (dsnameValue == ""){
#       errorCondition(recall=importMinitab,
#                      message=gettextRcmdr("You must enter the name of a data set."))
#       return()
#     }
#     if (!is.valid.name(dsnameValue)){
#       errorCondition(recall=importMinitab,
#                      message=paste('"', dsnameValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
#       return()
#     }
#     if (is.element(dsnameValue, listDataSets())) {
#       if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
#         importMinitab()
#         return()
#       }
#     }
#     file <- tclvalue(tkgetOpenFile(
#       filetypes=gettextRcmdr('{"All Files" {"*"}} {"Minitab portable files" {".mtp" ".MTP"}}')))
#     if (file == "") {
#       tkfocus(CommanderWindow())
#       return()
#     }
#     command <- paste('read.mtp("', file,'")', sep="")
#     datalist <- justDoIt(command)
#     lengths <- sapply(datalist, length)
#     datalist <- datalist[lengths != 0]
#     lengths <- lengths[lengths != 0]
#     if (!all(lengths == length(datalist[[1]]))){
#       Message(message=
#                 paste(gettextRcmdr("Minitab data set contains elements of unequal length.\nData set cannot be converted.")),
#               type="error")
#       tkdestroy(top)
#       tkfocus(CommanderWindow())
#       return()
#     }
#     #   	assign(dsnameValue, as.data.frame(datalist), envir=.GlobalEnv)
#     # 		logger(paste(dsnameValue, " <- as.data.frame(", command, ")", sep=""))
#     doItAndPrint(paste(dsnameValue, " <- as.data.frame(", command, ")", sep=""))
#     activeDataSet(dsnameValue)
#     tkfocus(CommanderWindow())
#   }
#   OKCancelHelp(helpSubject="read.mtp")
#   tkgrid(labelRcmdr(top, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="e")
#   tkgrid(buttonsFrame, columnspan="2", sticky="w")
#   tkgrid.configure(entryDsname, sticky="w")
#   dialogSuffix(rows=2, columns=2, focus=entryDsname)
# }
# 
# # imports data set from STATA file
# # the following function was contributed by Michael Ash (modified by J. Fox 2 Feb 05)
# importSTATA <- function() {
#   Library("foreign")
#   initializeDialog(title=gettextRcmdr("Import STATA Data Set"))
#   dsname <- tclVar(gettextRcmdr("Dataset"))
#   entryDsname <- ttkentry(top, width="20", textvariable=dsname)
#   asFactor <- tclVar("1")
#   asFactorCheckBox <- tkcheckbutton(top, variable=asFactor)
#   asDate <- tclVar("1")
#   asDateCheckBox <- tkcheckbutton(top, variable=asDate)
#   asMissingType <- tclVar("1")
#   asMissingTypeCheckBox <- tkcheckbutton(top, variable=asMissingType)
#   asConvertUnderscore <- tclVar("1")
#   asConvertUnderscoreCheckBox <- tkcheckbutton(top, variable=asConvertUnderscore)
#   asWarnMissingLabels <- tclVar("1")
#   asWarnMissingLabelsCheckBox <- tkcheckbutton(top, variable=asWarnMissingLabels)
#   onOK <- function(){
#     closeDialog()
#     dsnameValue <- trim.blanks(tclvalue(dsname))
#     if (dsnameValue == ""){
#       errorCondition(recall=importSTATA,
#                      message=gettextRcmdr("You must enter the name of a data set."))
#       return()
#     }
#     if (!is.valid.name(dsnameValue)){
#       errorCondition(recall=importSTATA,
#                      message=paste('"', dsnameValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
#       return()
#     }
#     if (is.element(dsnameValue, listDataSets())) {
#       if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
#         importSTATA()
#         return()
#       }
#     }
#     file <- tclvalue(tkgetOpenFile(
#       filetypes=gettextRcmdr('{"All Files" {"*"}} {"STATA datasets" {".dta" ".DTA"}}')))
#     if (file == "") {
#       tkfocus(CommanderWindow())
#       return()
#     }
#     convert.date <- tclvalue(asDate) == "1"
#     factor <- tclvalue(asFactor) == "1"
#     missingtype <- tclvalue(asMissingType) == "1"
#     convertunderscore <- tclvalue(asConvertUnderscore) == "1"
#     warnmissinglabels <- tclvalue(asWarnMissingLabels) == "1"
#     command <- paste('read.dta("', file,'", convert.dates=', convert.date,
#                      ", convert.factors=", factor, ", missing.type=", missingtype,
#                      ", convert.underscore=", convertunderscore, ", warn.missing.labels=TRUE)", sep="")
#     logger(paste(dsnameValue, " <- ", command, sep=""))
#     result <- justDoIt(command)
#     if (class(result)[1] !=  "try-error"){
#       gassign(dsnameValue, result)
#       activeDataSet(dsnameValue)
#     }
#     tkfocus(CommanderWindow())
#   }
#   OKCancelHelp(helpSubject="read.dta")
#   tkgrid(labelRcmdr(top, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="w")
#   tkgrid(labelRcmdr(top, text=gettextRcmdr("Convert value labels\nto factor levels"), justify="left"),
#          asFactorCheckBox, sticky="w")
#   tkgrid(labelRcmdr(top, text=gettextRcmdr("Convert dates to R format"), justify="left"),
#          asDateCheckBox, sticky="w")
#   tkgrid(labelRcmdr(top, text=gettextRcmdr("Multiple missing types (>=Stata 8)"), justify="left"),
#          asMissingTypeCheckBox, sticky="w")
#   tkgrid(labelRcmdr(top, text=gettextRcmdr("Convert underscore to period"), justify="left"),
#          asConvertUnderscoreCheckBox, sticky="w")
#   tkgrid(labelRcmdr(top, text=gettextRcmdr("Warn on missing labels"), justify="left"),
#          asWarnMissingLabelsCheckBox, sticky="w")
#   tkgrid(buttonsFrame, columnspan="2", sticky="w")
#   tkgrid.configure(entryDsname, sticky="w")
#   tkgrid.configure(asFactorCheckBox, sticky="w")
#   tkgrid.configure(asDateCheckBox, sticky="w")
#   tkgrid.configure(asMissingTypeCheckBox, sticky="w")
#   tkgrid.configure(asWarnMissingLabelsCheckBox, sticky="w")
#   dialogSuffix(rows=4, columns=2, focus=entryDsname)
# }
# 
# # import data set from Excel file
# importExcel <- function(){
#   Library("XLConnect")
#   initializeDialog(title = gettextRcmdr("Import Excel Data Set"))
#   dsname <- tclVar(gettextRcmdr("Dataset"))
#   entryDsname <- ttkentry(top, width = "35", textvariable = dsname)
#   onOK <- function(){
#     closeDialog()
#     dsnameValue <- trim.blanks(tclvalue(dsname))
#     if(dsnameValue == ""){
#       errorCondition(recall = importExcel,
#                      message = gettextRcmdr("You must enter the name of a data set."))
#       return()
#     }
#     if(!is.valid.name(dsnameValue)){
#       errorCondition(recall = importExcel,
#                      message = paste('"', dsnameValue, '" ',
#                                      gettextRcmdr("is not a valid name."), sep = ""))
#       return()
#     }
#     if(is.element(dsnameValue, listDataSets())){
#       if("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
#         importExcel()
#         return()
#       }
#     }
#     File <- tclvalue(tkgetOpenFile(filetypes = gettextRcmdr(
#       '{"All Files" {"*"}} {"MS Excel 2007 file" {".xlsx" ".XLSX"}} {"MS Excel file" {".xls" ".XLS"}}'
#     ), parent=CommanderWindow()))
#     if(File == ""){
#       tkfocus(CommanderWindow())
#       return()
#     }
#     command <- paste('loadWorkbook("', File, '")', sep="")
#     #         logger(paste(".Workbook <- ", command, sep=""))
#     #         justDoIt(paste('assign(".Workbook", ', command, ", envir=.GlobalEnv)", sep=""))
#     doItAndPrint(paste(".Workbook <- ", command, sep=""))
#     worksheets <- getSheets(.Workbook)
#     if(length(worksheets)>1)
#       worksheet <- tk_select.list(worksheets,
#                                   title = gettextRcmdr("Select one table"))
#     else
#       worksheet <- worksheets
#     if(worksheet == ""){
#       errorCondition(message=gettextRcmdr("No table selected"))
#       return()
#     }
#     command <- paste('readWorksheet(.Workbook, "', worksheet, '")', sep="")
#     logger(paste(dsnameValue, " <- ", command, sep=""))
#     result <- justDoIt(command)
#     if (class(result)[1] !=  "try-error"){
#       gassign(dsnameValue, result)
#       activeDataSet(dsnameValue)
#     }
#     logger("remove(.Workbook)")
#     justDoIt("remove(.Workbook, envir=.GlobalEnv)")
#   }
#   OKCancelHelp(helpSubject="readWorksheet")
#   tkgrid(labelRcmdr(top, text=gettextRcmdr("Enter name of data set:  ")),
#          entryDsname, sticky="e")
#   tkgrid(buttonsFrame, columnspan="2", sticky="w")
#   tkgrid.configure(entryDsname, sticky="w")
#   dialogSuffix(rows=2, columns=2, focus=entryDsname)
# }
# 
# # export data set
# exportDataSet <- function() {
#   dsname <- activeDataSet()
#   initializeDialog(title=gettextRcmdr("Export Active Data Set"))
#   checkBoxes(frame="optionsFrame", boxes=c("colnames", "rownames", "quotes"),
#              initialValues=rep(1,3), labels=gettextRcmdr(c("Write variable names:", "Write row names:", "Quotes around character values:")))
#   missingVariable <- tclVar("NA")
#   missingEntry <- ttkentry(optionsFrame, width="8", textvariable=missingVariable)
#   radioButtons(name="delimiter", buttons=c("spaces", "tabs", "commas"), labels=gettextRcmdr(c("Spaces", "Tabs", "Commas")),
#                title=gettextRcmdr("Field Separator"))
#   otherButton <- ttkradiobutton(delimiterFrame, variable=delimiterVariable, value="other")
#   otherVariable <- tclVar("")
#   otherEntry <- ttkentry(delimiterFrame, width="4", textvariable=otherVariable)
#   onOK <- function(){
#     closeDialog()
#     col <- tclvalue(colnamesVariable) == 1
#     row <- tclvalue(rownamesVariable) == 1
#     quote <- tclvalue(quotesVariable) == 1
#     delim <- tclvalue(delimiterVariable)
#     missing <- tclvalue(missingVariable)
#     sep <- if (delim == "tabs") "\\t"
#     else if (delim == "spaces") " "
#     else if (delim == "commas") ","
#     else trim.blanks(tclvalue(otherVariable))
#     saveFile <- tclvalue(tkgetSaveFile(filetypes=gettextRcmdr('{"All Files" {"*"}} {"Text Files" {".txt" ".TXT" ".dat" ".DAT" ".csv" ".CSV"}}'),
#                                        defaultextension="txt",
#                                        initialfile=paste(dsname, if (delim == "commas") ".csv" else ".txt", sep=""),
#                                        parent=CommanderWindow()))
#     if (saveFile == "") {
#       tkfocus(CommanderWindow())
#       return()
#     }
#     command <- paste("write.table(", dsname, ', "', saveFile, '", sep="', sep,
#                      '", col.names=', col, ", row.names=", row, ", quote=", quote,
#                      ', na="', missing, '")', sep="")
#     justDoIt(command)
#     logger(command)
#     Message(paste(gettextRcmdr("Active dataset exported to file"), saveFile), type="note")
#     tkfocus(CommanderWindow())
#   }
#   OKCancelHelp(helpSubject="write.table")
#   tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Missing values:")), missingEntry, sticky="w")
#   tkgrid(optionsFrame, sticky="w")
#   tkgrid(labelRcmdr(delimiterFrame, text=gettextRcmdr("Other")), otherButton,
#          labelRcmdr(delimiterFrame, text=gettextRcmdr("  Specify:")), otherEntry, sticky="w")
#   tkgrid(delimiterFrame, stick="w")
#   tkgrid(buttonsFrame, sticky="w")
#   dialogSuffix(rows=3, columns=1)
# }
# 
# # save active data set
# saveDataSet <- function() {
#   file <- tclvalue(tkgetSaveFile(filetypes=
#                                    gettextRcmdr('{"All Files" {"*"}} {"R Data Files" {".RData" ".rda" ".Rda" ".RDA"}}'),
#                                  defaultextension=".RData", initialfile=paste(activeDataSet(), ".RData", sep="")))
#   if (file == "") return()
#   command <- paste('save("', activeDataSet(), '", file="', file, '")', sep="")
#   justDoIt(command)
#   logger(command)
# }