# # Copied from R-forge Dec 1, 2015 by Jessica Reese
# # all credit goes to John Fox and those he lists as contributers
# ## Edit Data
# 
# # recoding variables within data set
# RecodeDialog <- function () {
#   require("car")
#   processRecode <- function(recode) {
#     parts <- strsplit(recode, "=")[[1]]
#     if (length(grep(",", parts[1])) > 0) 
#       paste("c(", parts[1], ") = ", parts[2], sep = "")
#     else paste(parts, collapse = "=")
#   }
#   dataSet <- activeDataSet()
#   defaults <- list (initial.asFactor = 1, initial.variables = NULL, initial.name = gettextRcmdr ("variable"),
#                     initial.recode.directives="")
#   dialog.values <- getDialog ("RecodeDialog", defaults)
#   initializeDialog(title = gettextRcmdr("Recode Variables"))
#   variablesBox <- variableListBox(top, Variables(), selectmode = "multiple", 
#                                   title = gettextRcmdr("Variables to recode (pick one or more)"),
#                                   initialSelection = varPosn (dialog.values$initial.variables, "all"))
#   variablesFrame <- tkframe(top)
#   newVariableName <- tclVar(dialog.values$initial.name)
#   newVariable <- ttkentry(variablesFrame, width = "20", textvariable = newVariableName)
#   recodesFrame <- tkframe(top)
#   recodes <- tktext(recodesFrame, bg = "white", font = getRcmdr("logFont"), 
#                     height = "5", width = "40", wrap = "none")
#   recodesXscroll <- ttkscrollbar(recodesFrame, orient = "horizontal", 
#                                  command = function(...) tkxview(recodes, ...))
#   recodesYscroll <- ttkscrollbar(recodesFrame, command = function(...) tkyview(recodes, 
#                                                                                ...))
#   tkconfigure(recodes, xscrollcommand = function(...) tkset(recodesXscroll, 
#                                                             ...))
#   tkconfigure(recodes, yscrollcommand = function(...) tkset(recodesYscroll, 
#                                                             ...))
#   tkinsert(recodes, "1.0", dialog.values$initial.recode.directives)
#   asFactorFrame <- tkframe(top)
#   asFactorVariable <- tclVar(dialog.values$initial.asFactor)
#   asFactorCheckBox <- tkcheckbutton(asFactorFrame, variable = asFactorVariable)
#   onOK <- function() {
#     asFactor <- tclvalue(asFactorVariable) == "1"
#     save.recodes <- trim.blanks(tclvalue(tkget(recodes, "1.0", "end")))
#     recode.directives <- gsub("\n", "; ", save.recodes)
#     check.empty <- gsub(";", "", gsub(" ", "", recode.directives))
#     if ("" == check.empty) {
#       errorCondition(recall = RecodeDialog, message = gettextRcmdr("No recode directives specified."))
#       return()
#     }
#     if (0 != length(grep("'", recode.directives))) {
#       errorCondition(recall = RecodeDialog, message = gettextRcmdr("Use only double-quotes (\" \") in recode directives"))
#       return()
#     }
#     recode.directives <- strsplit(recode.directives, ";")[[1]]
#     recode.directives <- paste(sapply(recode.directives, 
#                                       processRecode), collapse = ";")
#     recode.directives <- sub(" *; *$", "", recode.directives)
#     variables <- getSelection(variablesBox)
#     closeDialog()
#     if (length(variables) == 0) {
#       errorCondition(recall = RecodeDialog, message = gettextRcmdr("You must select a variable."))
#       return()
#     }
#     multiple <- if (length(variables) > 1) 
#       TRUE
#     else FALSE
#     name <- trim.blanks(tclvalue(newVariableName))
#     #        save.recodes <- gsub("; ", "\\\n", trim.blanks(recode.directives))  
#     putDialog ("RecodeDialog", list (initial.asFactor = asFactor, initial.variables = variables,
#                                      initial.name = name, initial.recode.directives=save.recodes))
#     for (variable in variables) {
#       newVar <- if (multiple) 
#         paste(name, variable, sep = "")
#       else name
#       if (!is.valid.name(newVar)) {
#         errorCondition(recall = RecodeDialog, message = paste("\"", 
#                                                               newVar, "\" ", gettextRcmdr("is not a valid name."), 
#                                                               sep = ""))
#         return()
#       }
#       if (is.element(newVar, Variables())) {
#         if ("no" == tclvalue(checkReplace(newVar))) {
#           RecodeDialog()
#           return()
#         }
#       }
#       cmd <- paste("Recode(", dataSet, "$", variable, ", '", 
#                    recode.directives, "', as.factor.result=", asFactor, 
#                    ")", sep = "")
#       logger(paste(dataSet, "$", newVar, " <- ", cmd, sep = ""))
#       result <- justDoIt(paste(dataSet, "$", newVar, " <- ", 
#                                cmd, sep = ""))
#       if (class(result)[1] != "try-error") 
#         activeDataSet(dataSet, flushModel = FALSE, flushDialogMemory = FALSE)
#       tkfocus(CommanderWindow())
#     }
#   }
#   OKCancelHelp(helpSubject = "RecodeDialog", reset = "RecodeDialog")
#   tkgrid(getFrame(variablesBox), sticky = "nw")
#   tkgrid(labelRcmdr(variablesFrame, text = ""))
#   tkgrid(labelRcmdr(variablesFrame, text = gettextRcmdr("New variable name or prefix for multiple recodes: ")), 
#          newVariable, sticky = "w")
#   tkgrid(labelRcmdr(asFactorFrame, text = gettextRcmdr("Make (each) new variable a factor")), 
#          asFactorCheckBox, sticky = "w")
#   tkgrid(labelRcmdr(asFactorFrame, text = ""))
#   tkgrid(labelRcmdr(recodesFrame, text = gettextRcmdr("Enter recode directives"), 
#                     fg = "blue"), sticky = "w")
#   tkgrid(recodes, recodesYscroll, sticky = "nw")
#   tkgrid(recodesXscroll)
#   tkgrid(variablesFrame, sticky = "w")
#   tkgrid(asFactorFrame, sticky = "w")
#   tkgrid(recodesFrame, sticky = "w")
#   tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#   tkgrid.configure(recodesXscroll, sticky = "ew")
#   tkgrid.configure(recodesYscroll, sticky = "ns")
#   dialogSuffix(rows = 4, columns = 2, bindReturn = FALSE)
# }
# 
# # computes a new variable within data set
# Compute <- function(){
#   onDoubleClick <-function(){
#     var <- trim.blanks(getSelection(variablesBox))
#     word <- paste("\\[", gettextRcmdr("factor"), "\\]", sep="")
#     if (length(grep(word, var)) == 1)
#       var <- trim.blanks(sub(word, "",  var))
#     tkfocus(compute)
#     expr <- tclvalue(computeVar)
#     tclvalue(computeVar) <- if (expr == "") var
#     else paste(expr, var, sep=if (rev(strsplit(expr, "")[[1]])[1] =="(" ) "" else " ")
#     tkicursor(compute, "end")
#     tkxview.moveto(compute, "1")
#   }
#   defaults <- list(initial.name = gettextRcmdr("variable"), initial.expression = "")
#   dialog.values <- getDialog("Compute", defaults)
#   dataSet <- activeDataSet()
#   initializeDialog(title=gettextRcmdr("Compute New Variable"))
#   .variables <- Variables()
#   variables <- paste(.variables, ifelse(is.element(.variables, Factors()), gettextRcmdr("[factor]"), ""))
#   variablesBox <- variableListBox(top, variables, title=gettextRcmdr("Current variables (double-click to expression)"))
#   tkbind(variablesBox$listbox, "<Double-ButtonPress-1>", onDoubleClick)
#   variablesFrame <- tkframe(top)
#   newVariableName <- tclVar(dialog.values$initial.name)
#   newVariable <- ttkentry(variablesFrame, width="20", textvariable=newVariableName)
#   computeFrame <- tkframe(top)
#   computeVar <- tclVar(dialog.values$initial.expression)
#   compute <- ttkentry(computeFrame, font=getRcmdr("logFont"), width="30", textvariable=computeVar)
#   computeXscroll <- ttkscrollbar(computeFrame,
#                                  orient="horizontal", command=function(...) tkxview(compute, ...))
#   tkconfigure(compute, xscrollcommand=function(...) tkset(computeXscroll, ...))
#   onOK <- function(){
#     closeDialog()
#     newVar <- trim.blanks(tclvalue(newVariableName))
#     if (!is.valid.name(newVar)){
#       errorCondition(recall=Compute,
#                      message=paste('"', newVar, '" ', gettextRcmdr("is not a valid name."), sep=""))
#       return()
#     }
#     express <- tclvalue(computeVar)
#     check.empty <- gsub(";", "", gsub(" ", "", express))
#     if ("" == check.empty) {
#       errorCondition(recall=Compute,
#                      message=gettextRcmdr("No expression specified."))
#       return()
#     }
#     putDialog("Compute", list(initial.name=newVar, initial.expression=express))
#     if (is.element(newVar, Variables())) {
#       if ("no" == tclvalue(checkReplace(newVar, gettextRcmdr("Variable")))){
#         Compute()
#         return()
#       }
#     }
#     command <-  paste(dataSet,"$",newVar, " <- with(", ActiveDataSet(),
#                       ", ", express, ")", sep="")
#     logger(command)
#     result <- justDoIt(command)
#     if (class(result)[1] !=  "try-error") activeDataSet(dataSet, flushModel=FALSE, flushDialogMemory=FALSE)
#     tkfocus(CommanderWindow())
#   }
#   OKCancelHelp(helpSubject="Compute", reset = "Compute")
#   tkgrid(getFrame(variablesBox), sticky="nw", columnspan=2)
#   tkgrid(labelRcmdr(variablesFrame, text=gettextRcmdr("New variable name")), sticky="w")
#   tkgrid(newVariable, labelRcmdr(variablesFrame, text="     "), sticky="w")
#   tkgrid(labelRcmdr(computeFrame, text=gettextRcmdr("Expression to compute")), sticky="w")
#   tkgrid(compute, sticky="w")
#   tkgrid(computeXscroll, sticky="ew")
#   tkgrid(variablesFrame, computeFrame, sticky="nw")
#   tkgrid(buttonsFrame, sticky="w", columnspan=2)
#   dialogSuffix(rows=3, columns=2, focus=compute)
# }
# 
# # standardize variables
# standardize <- function(X){
#   initializeDialog(title=gettextRcmdr("Standardize Variables"))
#   xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variables (pick one or more)"),
#                           selectmode="multiple")
#   onOK <- function(){
#     x <- getSelection(xBox)
#     closeDialog()
#     if (length(x) == 0) {
#       errorCondition(recall=standardize, message=gettextRcmdr("You must select one or more variables."))
#       return()
#     }
#     xx <- paste('"', x, '"', sep="")
#     .activeDataSet <- ActiveDataSet()
#     command <- paste("scale(", .activeDataSet, "[,c(", paste(xx, collapse=","),
#                      ")])", sep="")
#     result <- justDoIt(command)
#     gassign(".Z", result)
#     logger(paste(".Z <- ", command, sep=""))
#     for (i in 1:length(x)){
#       Z <- paste("Z.", x[i], sep="")
#       if (is.element(Z, Variables())) {
#         if ("no" == tclvalue(checkReplace(Z))){
#           if (GrabFocus()) tkgrab.release(top)
#           tkdestroy(top)
#           next
#         }
#       }
#       justDoIt(paste(.activeDataSet, "$", Z, " <- .Z[,", i, "]", sep=""))
#       logger(paste(.activeDataSet, "$", Z, " <- .Z[,", i, "]", sep=""))
#     }
#     remove(.Z, envir=.GlobalEnv)
#     logger("remove(.Z)")
#     if (class(result)[1] !=  "try-error") activeDataSet(.activeDataSet, flushModel=FALSE, flushDialogMemory=FALSE)
#     tkfocus(CommanderWindow())
#   }
#   OKCancelHelp(helpSubject="scale")
#   tkgrid(getFrame(xBox), sticky="w")
#   tkgrid(buttonsFrame, sticky="w")
#   dialogSuffix(rows=2, columns=1)
# }
# 
# # convert numeric variable to factor
# numericToFactor <- function(){
#   initializeDialog(title=gettextRcmdr("Convert Numeric Variables to Factors"))
#   variableBox <- variableListBox(top, Numeric(), selectmode="multiple",
#                                  title=gettextRcmdr("Variables (pick one or more)"))
#   radioButtons(name="levels", buttons=c("names", "numbers"),
#                labels=gettextRcmdr(c("Supply level names", "Use numbers")), title=gettextRcmdr("Factor Levels"))
#   factorName <- tclVar(gettextRcmdr("<same as variables>"))
#   factorNameField <- ttkentry(top, width="20", textvariable=factorName)
#   onOK <- function(){
#     variables <- getSelection(variableBox)
#     closeDialog()
#     if (length(variables) == 0) {
#       errorCondition(recall=numericToFactor, message=gettextRcmdr("You must select a variable."))
#       return()
#     }
#     facname <- trim.blanks(tclvalue(factorName))
#     .activeDataSet <- ActiveDataSet()
#     cmd <- paste("apply(", .activeDataSet, "[c(", paste(
#       paste('"', variables, '"', sep=""),
#       collapse=","), ")], 2, function(x) sort(unique(x)))", sep="")
#     levs <- eval(parse(text=cmd), envir=.GlobalEnv)
#     sameLevels <- (length(variables) == 1) ||
#       ((is.matrix(levs)) && (all(0 == apply(levs, 1, var))))
#     for (name in variables){
#       fname <- if (facname == gettextRcmdr("<same as variables>")) name
#       else if (length(variables) == 1) facname
#       else paste(facname, name, sep="")
#       if (!is.valid.name(fname)){
#         errorCondition(recall=numericToFactor,
#                        message=paste('"', fname, '" ', gettextRcmdr("is not a valid name."), sep=""))
#         return()
#       }
#       if (is.element(fname, Variables())) {
#         if ("no" == tclvalue(checkReplace(fname))){
#           numericToFactor()
#           return()
#         }
#       }
#       levelsType <- tclvalue(levelsVariable)
#       env <- environment()
#       if (((name == variables[1]) || (!sameLevels)) && (levelsType == "names")){
#         values <- sort(unique(eval(parse(text=paste(.activeDataSet, "$", name, sep="")),
#                                    envir=.GlobalEnv)))
#         nvalues <- length(values)
#         if (nvalues > 30) {
#           errorCondition(recall=numericToFactor,
#                          message=sprintf(gettextRcmdr("Number of levels (%d) too large."), nvalues))
#           return()
#         }
#         initializeDialog(subdialog,
#                          title=paste(gettextRcmdr("Level Names for"),
#                                      if(sameLevels && length(variables) > 1) "Factors" else fname))
#         names <- rep("", nvalues)
#         onOKsub <- function() {
#           closeDialog(subdialog)
#           for (i in 1:nvalues){
#             names[i] <- eval(parse(text=paste("tclvalue(levelName", i, ")", sep="")))
#           }
#           if (length(unique(names)) != nvalues){
#             errorCondition(recall=numericToFactor,
#                            message=gettextRcmdr("Levels names are not unique."))
#             return()
#           }
#           if (any(names == "")){
#             errorCondition(recall=numericToFactor,
#                            message=gettextRcmdr("A level name is empty."))
#             return()
#           }
#           assign("labels", paste(paste("'", names, "'", sep=""), collapse=","),
#                  envir=env)
#         }
#         subOKCancelHelp()
#         tkgrid(labelRcmdr(subdialog, text=gettextRcmdr("Numeric value")), labelRcmdr(subdialog, text=gettextRcmdr("Level name")), sticky="w")
#         for (i in 1:nvalues){
#           valVar <- paste("levelName", i, sep="")
#           assign(valVar, tclVar(""))
#           assign(paste("entry", i, sep=""), ttkentry(subdialog, width="20",
#                                                      textvariable=get(valVar)))
#           #                        textvariable=eval(parse(text=valVar))))
#           tkgrid(labelRcmdr(subdialog, text=values[i]), get(paste("entry", i, sep="")), sticky="w")
#           #                    tkgrid(labelRcmdr(subdialog, text=values[i]), eval(parse(text=paste("entry", i, sep=""))), sticky="w")
#         }
#         tkgrid(subButtonsFrame, sticky="w", columnspan=2)
#         dialogSuffix(subdialog, rows=nvalues+2, columns=2, focus=entry1, onOK=onOKsub)
#       }
#       if (levelsType == "names"){
#         if (!exists("labels", mode="character")) return()
#         command <- paste("factor(", .activeDataSet, "$", name,
#                          ", labels=c(", labels, "))", sep="")
#         result <- justDoIt(paste(.activeDataSet, "$", fname, " <- ", command, sep=""))
#         logger(paste(.activeDataSet,"$", fname," <- ", command, sep=""))
#         if (class(result)[1] !=  "try-error") activeDataSet(.activeDataSet)
#         tkfocus(CommanderWindow())
#       }
#       else{
#         command <- paste("as.factor(", .activeDataSet, "$", name, ")", sep="")
#         result <- justDoIt(paste(.activeDataSet, "$", fname, " <- ", command, sep=""))
#         logger(paste(.activeDataSet, "$", fname," <- ", command, sep=""))
#         if (class(result)[1] !=  "try-error") activeDataSet(.activeDataSet, flushModel=FALSE, flushDialogMemory=FALSE)
#         tkfocus(CommanderWindow())
#       }
#     }
#   }
#   OKCancelHelp(helpSubject="factor")
#   tkgrid(getFrame(variableBox), levelsFrame, sticky="nw")
#   tkgrid(labelRcmdr(top,
#                     text=gettextRcmdr("New variable name or prefix for multiple variables:")),
#          factorNameField, sticky="w")
#   tkgrid(buttonsFrame, sticky="w", columnspan=2)
#   tkgrid.configure(numbersButton, sticky="w")
#   tkgrid.configure(namesButton, sticky="w")
#   dialogSuffix(rows=4, columns=2, preventGrabFocus=TRUE)
# }
# 
# # bin (section off/categorize) numeric variables
# binVariable <- function () {
#   # Author: Dan Putler (revision by J. Fox, 2 Feb 05)
#   defaults <- list (initial.levels = "specify", initial.bins = "3", initial.varName = NULL, 
#                     initial.newVar = gettextRcmdr("variable"), initial.method = "intervals")
#   dialog.values <- getDialog ("binVariable", defaults)
#   env <- environment()
#   initializeDialog(title = gettextRcmdr("Bin a Numeric Variable"))
#   variableFrame <- tkframe(top)
#   variableBox <- variableListBox(variableFrame, Numeric(), 
#                                  title = gettextRcmdr("Variable to bin (pick one)"), 
#                                  initialSelection = varPosn (dialog.values$initial.varName, "numeric"))
#   newVariableFrame <- tkframe(variableFrame)
#   newVariableName <- tclVar(dialog.values$initial.newVar)
#   newVariable <- ttkentry(newVariableFrame, width = "18", textvariable = newVariableName)
#   binsFrame <- tkframe(top)
#   binsVariable <- tclVar(dialog.values$initial.bins)
#   slider <- tkscale(binsFrame, from = 2, to = 20, showvalue = TRUE, 
#                     variable = binsVariable, resolution = 1, orient = "horizontal")
#   optionsFrame <- tkframe(top)
#   radioButtons(optionsFrame, name = "levels", buttons = c("specify", 
#                                                           "numbers", "ranges"), labels = gettextRcmdr(c("Specify names", 
#                                                                                                         "Numbers", "Ranges")), title = gettextRcmdr("Level Names"),
#                initialValue = dialog.values$initial.levels)
#   radioButtons(optionsFrame, name = "method", buttons = c("intervals", 
#                                                           "proportions", "natural"), labels = gettextRcmdr(c("Equal-width bins", 
#                                                                                                              "Equal-count bins", "Natural breaks\n(from K-means clustering)")), 
#                title = gettextRcmdr("Binning Method"), 
#                initialValue = dialog.values$initial.method)
#   onOK <- function() {
#     levels <- tclvalue(levelsVariable)
#     bins <- as.numeric(tclvalue(binsVariable))
#     varName <- getSelection(variableBox)
#     closeDialog()
#     if (length(varName) == 0) {
#       errorCondition(recall = binVariable, message = gettextRcmdr("You must select a variable."))
#       return()
#     }
#     newVar <- tclvalue(newVariableName)
#     if (is.element(newVar, Variables())) {
#       if ("no" == tclvalue(checkReplace(newVar))) {
#         binVariable()
#         return()
#       }
#     }
#     if (!is.valid.name(newVar)) {
#       errorCondition(message = paste("\"", newVar, "\" ", 
#                                      gettextRcmdr("is not a valid name."), sep = ""), 
#                      recall = binVariable)
#       return()
#     }
#     method <- tclvalue(methodVariable)
#     putDialog ("binVariable", list (initial.levels = levels, initial.bins = bins, initial.varName = varName, 
#                                     initial.newVar = newVar, initial.method = method))
#     if (levels == "specify") {
#       initializeDialog(subdialog, title = gettextRcmdr("Bin Names"))
#       onOKsub <- function() {
#         closeDialog(subdialog)
#         level <- character(bins)
#         for (i in 1:bins) {
#           level[i] <- eval(parse(text = paste("tclvalue(levelName", 
#                                               i, ")", sep = "")))
#         }
#         if (length(unique(level)) != length(level)) {
#           errorCondition(window = subdialog, message = gettextRcmdr("Level names must be unique."), 
#                          recall = onOK)
#           return()
#         }
#         assign("levelNames", level, envir = env)
#       }
#       subOKCancelHelp()
#       tkgrid(labelRcmdr(subdialog, text = gettextRcmdr("Bin"), 
#                         fg = "blue"), labelRcmdr(subdialog, text = gettextRcmdr("Name"), 
#                                                  fg = "blue"), sticky = "w")
#       for (i in 1:bins) {
#         valVar <- paste("levelName", i, sep = "")
#         assign(valVar, tclVar(i))
#         assign(paste("entry", i, sep = ""), ttkentry(subdialog, 
#                                                      width = "20", textvariable = get(valVar)))
#         tkgrid(labelRcmdr(subdialog, text = as.character(i)), 
#                get(paste("entry", i, sep = "")), sticky = "w")
#       }
#       tkgrid(subButtonsFrame, sticky = "w", columnspan = 2)
#       dialogSuffix(subdialog, focus = entry1, rows = bins + 
#                      1, columns = 2, bindReturn = FALSE)
#     }
#     labels <- if (levels == "numbers") 
#       "FALSE"
#     else if (levels == "ranges") 
#       "NULL"
#     else {
#       if (!exists("levelNames")) {
#         onCancel()
#         binVariable()
#         return()
#       }
#       paste("c('", paste(levelNames, collapse = "','"), 
#             "')", sep = "")
#     }
#     .activeDataSet <- ActiveDataSet()
#     command <- paste(.activeDataSet, "$", newVar, " <- ", 
#                      "bin.var(", .activeDataSet, "$", varName, ", bins=", 
#                      bins, ", method=", "'", method, "', labels=", labels, 
#                      ")", sep = "")
#     logger(command)
#     result <- justDoIt(command)
#     if (class(result)[1] != "try-error") 
#       activeDataSet(.activeDataSet, flushModel = FALSE, 
#                     flushDialogMemory = FALSE)
#     tkfocus(CommanderWindow())
#   }
#   OKCancelHelp(helpSubject = "bin.var", reset = "binVariable")
#   tkgrid(labelRcmdr(newVariableFrame, text = gettextRcmdr("New variable name"), 
#                     fg = "blue"), sticky = "w")
#   tkgrid(newVariable, sticky = "w")
#   tkgrid(getFrame(variableBox), labelRcmdr(variableFrame, text = "    "), 
#          newVariableFrame, sticky = "nw")
#   tkgrid(variableFrame, sticky = "w")
#   tkgrid(labelRcmdr(binsFrame, text = gettextRcmdr("Number of bins:")), 
#          slider, sticky = "s")
#   tkgrid(binsFrame, sticky = "w")
#   tkgrid(levelsFrame, labelRcmdr(optionsFrame, text = "    "), 
#          methodFrame, sticky = "nw")
#   tkgrid(optionsFrame, sticky = "w")
#   tkgrid(buttonsFrame, sticky = "w")
#   dialogSuffix(rows = 4, columns = 1)
# }
# 
# # reorder factor levels
# reorderFactor <- function(){
#   initializeDialog(title=gettextRcmdr("Reorder Factor Levels"))
#   variableBox <- variableListBox(top, Factors(), title=gettextRcmdr("Factor (pick one)"))
#   orderedFrame <- tkframe(top)
#   orderedVariable <- tclVar("0")
#   orderedCheckBox <- tkcheckbutton(orderedFrame, variable=orderedVariable)
#   factorName <- tclVar(gettextRcmdr("<same as original>"))
#   factorNameField <- ttkentry(top, width="20", textvariable=factorName)
#   onOK <- function(){
#     variable <- getSelection(variableBox)
#     closeDialog()
#     if (length(variable) == 0) {
#       errorCondition(recall=reorderFactor, message=gettextRcmdr("You must select a variable."))
#       return()
#     }
#     name <- trim.blanks(tclvalue(factorName))
#     if (name == gettextRcmdr("<same as original>")) name <- variable
#     if (!is.valid.name(name)){
#       errorCondition(recall=reorderFactor,
#                      message=paste('"', name, '" ', gettextRcmdr("is not a valid name."), sep=""))
#       return()
#     }
#     if (is.element(name, Variables())) {
#       if ("no" == tclvalue(checkReplace(name))){
#         reorderFactor()
#         return()
#       }
#     }
#     .activeDataSet <- ActiveDataSet()
#     old.levels <- eval(parse(text=paste("levels(", .activeDataSet, "$", variable, ")",
#                                         sep="")), envir=.GlobalEnv)
#     nvalues <- length(old.levels)
#     ordered <- tclvalue(orderedVariable)
#     if (nvalues > 30) {
#       errorCondition(recall=reorderFactor,
#                      message=sprintf(gettextRcmdr("Number of levels (%d) too large."), nvalues))
#       return()
#     }
#     initializeDialog(subdialog, title=gettextRcmdr("Reorder Levels"))
#     order <- 1:nvalues
#     onOKsub <- function() {
#       closeDialog(subdialog)
#       opt <- options(warn=-1)
#       for (i in 1:nvalues){
#         order[i] <- as.numeric(eval(parse(text=paste("tclvalue(levelOrder", i, ")", sep=""))))
#       }
#       options(opt)
#       if (any(sort(order) != 1:nvalues) || any(is.na(order))){
#         errorCondition(recall=reorderFactor,
#                        message=paste(gettextRcmdr("Order of levels must include all integers from 1 to "), nvalues, sep=""))
#         return()
#       }
#       levels <- old.levels[order(order)]
#       ordered <- if (ordered == "1") ", ordered=TRUE" else ""
#       command <- paste("factor(", .activeDataSet, "$", variable,
#                        ", levels=c(", paste(paste("'", levels, "'", sep=""), collapse=","), ")",
#                        ordered, ")", sep="")
#       result <- justDoIt(paste(.activeDataSet, "$", name, " <- ", command, sep=""))
#       logger(paste(.activeDataSet,"$", name," <- ", command, sep=""))
#       if (class(result)[1] !=  "try-error") activeDataSet(.activeDataSet, flushModel=FALSE, flushDialogMemory=FALSE)
#     }
#     subOKCancelHelp()
#     tkgrid(labelRcmdr(subdialog, text=gettextRcmdr("Old Levels"), fg="blue"),
#            labelRcmdr(subdialog, text=gettextRcmdr("New order"), fg="blue"), sticky="w")
#     for (i in 1:nvalues){
#       valVar <- paste("levelOrder", i, sep="")
#       assign(valVar, tclVar(i))
#       assign(paste("entry", i, sep=""), ttkentry(subdialog, width="2",
#                                                  textvariable=get(valVar)))
#       tkgrid(labelRcmdr(subdialog, text=old.levels[i]), get(paste("entry", i, sep="")), sticky="w")
#     }
#     tkgrid(subButtonsFrame, sticky="w", columnspan=2)
#     dialogSuffix(subdialog, focus=entry1, rows=nvalues+1, columns=2)
#   }
#   OKCancelHelp(helpSubject="factor")
#   tkgrid(getFrame(variableBox), sticky="nw")
#   tkgrid(labelRcmdr(top, text=gettextRcmdr("Name for factor")), sticky="w")
#   tkgrid(factorNameField, sticky="w")
#   tkgrid(labelRcmdr(orderedFrame, text=gettextRcmdr("Make ordered factor")), orderedCheckBox, sticky="w")
#   tkgrid(orderedFrame, sticky="w")
#   tkgrid(buttonsFrame, sticky="w")
#   dialogSuffix(rows=5, columns=1, preventGrabFocus=TRUE)
# }
# 
# # set contrasts for a factor
# setContrasts <- function(){
#   initializeDialog(title=gettextRcmdr("Set Contrasts for Factor"))
#   variableBox <- variableListBox(top, Factors(), title=gettextRcmdr("Factor (pick one)"))
#   radioButtons(name="contrasts", buttons=c("treatment", "sum", "helmert", "poly", "specify"),
#                values=c("contr.Treatment", "contr.Sum", "contr.helmert", "contr.poly", "specify"),
#                labels=gettextRcmdr(c("Treatment (dummy) contrasts", "Sum (deviation) contrasts", "Helmert contrasts",
#                                      "Polynomial contrasts", "Other (specify)")), title=gettextRcmdr("Contrasts"))
#   onOK <- function(){
#     variable <- getSelection(variableBox)
#     closeDialog()
#     if (length(variable) == 0) {
#       errorCondition(recall=setContrasts, message=gettextRcmdr("You must select a variable."))
#       return()
#     }
#     contrasts <- tclvalue(contrastsVariable)
#     if (contrasts != "specify"){
#       command <- paste("contrasts(", ActiveDataSet(), "$", variable, ') <- "', contrasts, '"', sep="")
#       result <- justDoIt(command)
#       logger(command)
#       if (class(result)[1] !=  "try-error") activeDataSet(ActiveDataSet())
#       tkfocus(CommanderWindow())
#     }
#     else{
#       initializeDialog(subdialog, title=gettextRcmdr("Specify Contrasts"))
#       tkgrid(labelRcmdr(subdialog, text=gettextRcmdr("Enter Contrast Coefficients"), fg="blue"), sticky="w")
#       env <- environment()
#       tableFrame <- tkframe(subdialog)
#       row.names <- eval(parse(text=paste("levels(", ActiveDataSet(), "$", variable, ")")))
#       row.names <- substring(paste(abbreviate(row.names, 12), "            "), 1, 12)
#       nrows <- length(row.names)
#       ncols <- nrows - 1
#       make.col.names <- paste("labelRcmdr(tableFrame, text='", gettextRcmdr("Contrast Name:"), "')", sep="")
#       for (j in 1:ncols) {
#         varname <- paste(".col.", j, sep="")
#         assign(varname, tclVar(paste(".", j, sep="")), envir=env)
#         make.col.names <- paste(make.col.names, ", ",
#                                 "ttkentry(tableFrame, width='12', textvariable=", varname, ")", sep="")
#       }
#       eval(parse(text=paste("tkgrid(", make.col.names, ", sticky='w')", sep="")), envir=env)
#       for (i in 1:nrows){
#         make.row <- paste("labelRcmdr(tableFrame, text='", row.names[i], "')")
#         for (j in 1:ncols){
#           varname <- paste(".tab.", i, ".", j, sep="")
#           assign(varname, tclVar("0"), envir=env)
#           make.row <- paste(make.row, ", ", "ttkentry(tableFrame, width='5', textvariable=",
#                             varname, ")", sep="")
#         }
#         eval(parse(text=paste("tkgrid(", make.row, ", sticky='w')", sep="")), envir=env)
#       }
#       tkgrid(tableFrame, sticky="w")
#       onOKsub <- function(){
#         closeDialog(subdialog)
#         cell <- 0
#         values <- rep(NA, nrows*ncols)
#         for (j in 1:ncols){
#           for (i in 1:nrows){
#             cell <- cell + 1
#             varname <- paste(".tab.", i, ".", j, sep="")
#             values[cell] <- as.numeric(eval(parse(text=paste("tclvalue(", varname,")", sep=""))))
#           }
#         }
#         values <- na.omit(values)
#         if (length(values) != nrows*ncols){
#           errorCondition(subdialog, recall=setContrasts,
#                          message=sprintf(gettextRcmdr(
#                            "Number of valid entries in contrast matrix(%d)\nnot equal to number of levels (%d) * number of contrasts (%d)."), length(values), nrows, ncols))
#           return()
#         }
#         if (qr(matrix(values, nrows, ncols))$rank < ncols) {
#           errorCondition(subdialog, recall=setContrasts, message=gettextRcmdr("Contrast matrix is not of full column rank"))
#           return()
#         }
#         contrast.names <- rep("", ncols)
#         for (j in 1:ncols){
#           varname <- paste(".col.", j, sep="")
#           contrast.names[j] <- eval(parse(text=paste("tclvalue(", varname,")", sep="")))
#         }
#         if (length(unique(contrast.names)) < ncols) {
#           errorCondition(subdialog, recall=setContrasts, message=gettextRcmdr("Contrast names must be unique"))
#           return()
#         }
#         command <- paste("matrix(c(", paste(values, collapse=","), "), ", nrows, ", ", ncols,
#                          ")", sep="")
#         #     		assign(".Contrasts", justDoIt(command), envir=.GlobalEnv)
#         # 				logger(paste(".Contrasts <- ", command, sep=""))
#         doItAndPrint(paste(".Contrasts <- ", command, sep=""))
#         command <- paste("colnames(.Contrasts) <- c(",
#                          paste("'", contrast.names, "'", sep="", collapse=", "), ")", sep="")
#         justDoIt(command)
#         logger(command)
#         command <- paste("contrasts(", ActiveDataSet(), "$", variable, ") <- .Contrasts", sep="")
#         result <- justDoIt(command)
#         logger(command)
#         justDoIt("remove(.Contrasts, envir=.GlobalEnv)")
#         logger("remove(.Contrasts)")
#         if (class(result)[1] !=  "try-error") activeDataSet(ActiveDataSet(), flushModel=FALSE, flushDialogMemory=FALSE)
#         tkfocus(CommanderWindow())
#       }
#       subOKCancelHelp(helpSubject="contrasts")
#       tkgrid(tableFrame, sticky="w")
#       tkgrid(labelRcmdr(subdialog, text=""))
#       tkgrid(subButtonsFrame, sticky="w")
#       dialogSuffix(subdialog, rows=5, columns=1, focus=subdialog)
#     }
#   }
#   OKCancelHelp(helpSubject="contrasts")
#   tkgrid(getFrame(variableBox), sticky="nw")
#   tkgrid(contrastsFrame, sticky="w")
#   tkgrid(buttonsFrame, sticky="w")
#   dialogSuffix(rows=4, columns=1)
# }
# 
# # rename existing variable
# renameVariables <- function(){
#   initializeDialog(title=gettextRcmdr("Rename Variables"))
#   variableBox <- variableListBox(top, Variables(), title=gettextRcmdr("Variables (pick one or more)"),
#                                  selectmode="multiple", initialSelection=NULL)
#   onOK <- function(){
#     variables <- getSelection(variableBox)
#     closeDialog()
#     nvariables <- length(variables)
#     if (nvariables < 1) {
#       errorCondition(recall=renameVariables, message=gettextRcmdr("No variables selected."))
#       return()
#     }
#     .activeDataSet <- ActiveDataSet()
#     unordered.names <- names(get(.activeDataSet))
#     which.variables <- match(variables, unordered.names)
#     initializeDialog(subdialog, title=gettextRcmdr("Variable Names"))
#     newnames <- rep("", nvariables)
#     onOKsub <- function() {
#       closeDialog(subdialog)
#       for (i in 1:nvariables){
#         newnames[i] <- eval(parse(text=paste("tclvalue(newName", i, ")", sep="")))
#       }
#       if (any(newnames == "")){
#         errorCondition(recall=renameVariables, message=gettextRcmdr("A variable name is empty."))
#         return()
#       }
#       test.names <- newnames == make.names(newnames)
#       if (!all(test.names)){
#         errorCondition(recall=renameVariables,
#                        message=paste(gettextRcmdr("The following variable names are not valid:\n"),
#                                      paste(newnames[!test.names], collapse=", ")))
#         return()
#       }
#       all.names <- names(get(.activeDataSet))
#       #            all.names <- eval(parse(text=paste("names(", .activeDataSet, ")")))
#       all.names[which.variables] <- newnames
#       if (length(unique(all.names)) != length(all.names)){
#         errorCondition(recall=renameVariables, message=gettextRcmdr("Variable names are not unique"))
#         return()
#       }
#       command <- paste("names(", .activeDataSet, ")[c(", paste(which.variables, collapse=","),
#                        ")] <- c(", paste('"', newnames, '"', collapse=",", sep=""), ")", sep="")
#       result <- justDoIt(command)
#       logger(command)
#       if (class(result)[1] !=  "try-error") activeDataSet(.activeDataSet, flushModel=FALSE, flushDialogMemory=FALSE)
#       tkfocus(CommanderWindow())
#     }
#     subOKCancelHelp()
#     tkgrid(labelRcmdr(subdialog, text=gettextRcmdr("Old Name"), fg="blue"),
#            labelRcmdr(subdialog, text=gettextRcmdr("New name"), fg="blue"), sticky="w")
#     for (i in 1:nvariables){
#       valVar <- paste("newName", i, sep="")
#       assign(valVar, tclVar(""))
#       assign(paste("entry", i, sep=""), ttkentry(subdialog, width="20",
#                                                  textvariable=get(valVar)))
#       tkgrid(labelRcmdr(subdialog, text=variables[i]), get(paste("entry", i, sep="")), sticky="w")
#     }
#     tkgrid(subButtonsFrame, sticky="w", columnspan=2)
#     dialogSuffix(subdialog, rows=nvariables+2, columns=2, focus=entry1, onOK=onOKsub)
#   }
#   OKCancelHelp(helpSubject="names")
#   tkgrid(getFrame(variableBox), sticky="nw")
#   tkgrid(buttonsFrame, sticky="w")
#   dialogSuffix(rows=2, columns=1)
# }
# 
# # deletes a variable within a data set
# deleteVariable <- function(){
#   dataSet <- activeDataSet()
#   initializeDialog(title=gettextRcmdr("Delete Variables"))
#   variablesBox <- variableListBox(top, Variables(),
#                                   title=gettextRcmdr("Variable(s) to delete (pick one or more)"), selectmode="multiple",
#                                   initialSelection=NULL)
#   onOK <- function(){
#     variables <- getSelection(variablesBox)
#     closeDialog()
#     if (length(variables) == 0) {
#       errorCondition(recall=deleteVariable, message=gettextRcmdr("You must select one or more variables."))
#       return()
#     }
#     if (length(variables) == 1){
#       response <- tclvalue(RcmdrTkmessageBox(message=sprintf(gettextRcmdr("Delete %s?\nPlease confirm."), variables), icon="warning", type="okcancel", default="cancel"))
#       if (response == "cancel") {
#         onCancel()
#         return()
#       }
#     }
#     else{
#       response <- tclvalue(RcmdrTkmessageBox(message=
#                                                sprintf(gettextRcmdr("Delete %d variables?\nPlease confirm."), length(variables)),
#                                              icon="warning", type="okcancel", default="cancel"))
#       if (response == "cancel") {
#         onCancel()
#         return()
#       }
#     }
#     for (variable in variables){
#       eval(parse(text=paste(dataSet, "$", variable, "<- NULL", sep="")), envir=.GlobalEnv)
#       logger(paste(dataSet, "$", variable, " <- NULL", sep=""))
#     }
#     activeDataSet(dataSet, flushModel=FALSE, flushDialogMemory=FALSE)
#     tkfocus(CommanderWindow())
#   }
#   OKCancelHelp(helpSubject="NULL")
#   tkgrid(getFrame(variablesBox), sticky="nw")
#   tkgrid(buttonsFrame, sticky="w")
#   dialogSuffix(rows=2, columns=1)
# }
# 
# # Remove rows in active dataset
# RemoveRows <- function(){
#   dataSet <- activeDataSet()
#   initializeDialog(title=gettextRcmdr("Remove Rows from Active Data Set"))
#   removeVariable <- tclVar(gettextRcmdr(""))
#   removeFrame <- tkframe(top)
#   removeEntry <- ttkentry(removeFrame, width="60", textvariable=removeVariable)
#   removeScroll <- ttkscrollbar(removeFrame, orient="horizontal",
#                                command=function(...) tkxview(removeEntry, ...))
#   tkconfigure(removeEntry, xscrollcommand=function(...) tkset(removeScroll, ...))
#   newDataSetName <- tclVar(gettextRcmdr("<same as active data set>"))
#   dataSetNameFrame <- tkframe(top)
#   dataSetNameEntry <- ttkentry(dataSetNameFrame, width="25", textvariable=newDataSetName)
#   onOK <- function(){
#     newName <- trim.blanks(tclvalue(newDataSetName))
#     if (newName == gettextRcmdr("<same as active data set>")) newName <- ActiveDataSet()
#     if (!is.valid.name(newName)){
#       errorCondition(recall=RemoveRows,
#                      message=paste('"', newName, '" ', gettextRcmdr("is not a valid name."), sep=""))
#       return()
#     }
#     if (is.element(newName, listDataSets())) {
#       if ("no" == tclvalue(checkReplace(newName, type=gettextRcmdr("Data set")))){
#         closeDialog()
#         RemoveRows()
#         return()
#       }
#     }
#     remove <- tclvalue(removeVariable)
#     if (remove==""){
#       errorCondition(recall=RemoveRows,
#                      message="No rows to remove")
#       closeDialog()
#       return()
#     }
#     removeRows <- paste("c(", gsub(" ", ",", remove), ")", sep="")
#     remove <- try(eval(parse(text=removeRows)), silent=TRUE)
#     if (class(remove) == "try-error"){
#       errorCondition(recall=RemoveRows,
#                      message=remove)
#       closeDialog()
#       return()
#     }
#     closeDialog()
#     removeRows <- if (is.numeric(remove)) paste("-", removeRows, sep="") 
#     else paste("!(rownames(", ActiveDataSet(), ") %in% ", removeRows, ")", sep="")
#     command <- paste(newName, " <- ", ActiveDataSet(), "[", removeRows, ",]", sep="")
#     logger(command)
#     result <- justDoIt(command)
#     if (class(result)[1] !=  "try-error") activeDataSet(newName)
#     tkfocus(CommanderWindow())
#   }
#   OKCancelHelp(helpSubject="[.data.frame")
#   tkgrid(labelRcmdr(removeFrame, text=gettextRcmdr("Indices or quoted names of row(s) to remove"),
#                     foreground=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w")
#   tkgrid(removeEntry, sticky="w")
#   tkgrid(removeScroll, sticky="ew")
#   tkgrid(removeFrame, sticky="w")
#   tkgrid(labelRcmdr(dataSetNameFrame, text=gettextRcmdr("Name for new data set")), sticky="w")
#   tkgrid(dataSetNameEntry, sticky="w")
#   tkgrid(dataSetNameFrame, sticky="w")
#   tkgrid(buttonsFrame, sticky="w")
#   dialogSuffix()
# }