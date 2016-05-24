cchalandise.DEP0 <- function () {
Library("maps")
defaults <- list(initial.response = NULL, initial.alternative = "two.sided", 
initial.confidenceLevel = ".25", initial.variances = "FALSE")
dialog.values <- getDialog("cchalandise.DEP0", defaults)
initializeDialog(title = gettextRcmdr(paste("Chalandise par d", "\U00E9", "partement", 
            sep = "")))
variablesFrame <- tkframe(top)
responseBox <- variableListBox(variablesFrame, Numeric(), 
title = gettextRcmdr(paste("Code d", "\U00E9", "partement (ex : 69) ou code postal (ex : 69002)", 
            sep = "")),
initialSelection = varPosn(dialog.values$initial.response, "numeric"))
onOK <- function() {
response <- getSelection(responseBox)
if (length(response) == 0) {
errorCondition(recall = cchalandise.DEP0, 
message = gettextRcmdr("You must select a response variable."))
return()
}
alternative <- as.character(tclvalue(alternativeVariable))
level <- tclvalue(confidenceLevel)
variances <- as.character(tclvalue(variancesVariable))
putDialog ("cchalandise.DEP0", list (initial.response = response, initial.alternative = alternative, 
initial.confidenceLevel = level, initial.variances = variances))        
closeDialog()


command <- paste("chalandise.DEP0(", ActiveDataSet(), 
              "$", response, ",", level, ",arrondi=",variances,")", sep = "")
doItAndPrint(command)


tkfocus(CommanderWindow())
}
OKCancelHelp(helpSubject = "map", reset = "cchalandise.DEP0")
optionsFrame <- tkframe(top)
radioButtons(optionsFrame, name = "alternative", buttons = c("twosided", 
"less", "greater"), values = c("two.sided", "less", "greater"), 
labels = gettextRcmdr(c("Two-sided", "Difference < 0", 
"Difference > 0")), title = gettextRcmdr("Alternative Hypothesis"),
initialValue = dialog.values$initial.alternative)
confidenceFrame <- tkframe(optionsFrame)
confidenceLevel <- tclVar(dialog.values$initial.confidenceLevel)
confidenceField <- ttkentry(confidenceFrame, width = "6", 
textvariable = confidenceLevel)
radioButtons(optionsFrame, name = "variances", buttons = c("yes", 
"no"), values = c("TRUE", "FALSE"),  
labels = gettextRcmdr(c("Code postal", paste("D", "\U00E9", "partement", 
            sep = ""))),title = gettextRcmdr("Code?"),
initialValue = dialog.values$initial.variances)
tkgrid(getFrame(responseBox), labelRcmdr(variablesFrame, text = "    "), 
sticky = "nw")
tkgrid(variablesFrame, sticky = "nw")
tkgrid(labelRcmdr(confidenceFrame, text = gettextRcmdr("Taille des disques"), 
fg = "blue"), sticky = "w")
tkgrid(confidenceField, sticky = "w")
tkgrid(   variancesFrame,
confidenceFrame, labelRcmdr(optionsFrame, text = "    "), 
sticky = "nw")
tkgrid(optionsFrame, sticky = "nw")
tkgrid(buttonsFrame, sticky = "w")
dialogSuffix(rows = 4, columns = 1)
}
