# Some Rcmdr useful extensions

# Note (last modified: 2013-01-24 by J. Fox): the following function (with contributions from Richard Heiberger and Milan Bouchet-Valat) 
# can be included in any Rcmdr plug-in package to cause the package to load
# the Rcmdr if it is not already loaded
.onAttach <- function(libname, pkgname){
    if (!interactive()) return()
    putRcmdr("slider.env", new.env())    
    Rcmdr <- options()$Rcmdr
    plugins <- Rcmdr$plugins
    if (!pkgname %in% plugins) {
        Rcmdr$plugins <- c(plugins, pkgname)
        options(Rcmdr=Rcmdr)
        if("package:Rcmdr" %in% search()) {
            if(!getRcmdr("autoRestart")) {
                closeCommander(ask=FALSE, ask.save=TRUE)
                Commander()
            }
        }
        else {
            Commander()
        }
    }
}

# Function to be called by Rcmdr to test for randomness using runs.test from tseries packages
randomnessFTest <- function()
  {
  # To ensure that menu name is included in pot file
  gettext("Randomness test for two level factor...", domain="R-RcmdrPlugin.UCA")
  # Build dialog
  initializeDialog(title=gettext("Randomness test for two level factor", domain="R-RcmdrPlugin.UCA"))
  variablesBox <- variableListBox(top, TwoLevelFactors(), selectmode="single", initialSelection=NULL, title=gettextRcmdr("Variable (pick one)"))
	onOK <- function(){
		x <- getSelection(variablesBox)
		if (length(x) == 0) {
                    errorCondition(recall=randomnessNTest, message=gettextRcmdr("No variables were selected."))
                    return()
                }
		closeDialog()
                # Apply test
                doItAndPrint(paste("with(", ActiveDataSet(), ", twolevelfactor.runs.test(", x, "))", sep = ""))
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="runs.test", reset = "randomnessFTest", apply = "randomnessFTest")
	tkgrid(getFrame(variablesBox), sticky="nw")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=6, columns=1)
}

# Function to be called by Rcmdr to test for randomness using runs.test from randtest packages
randomnessNTest <- function()
  {
  # To ensure that menu name is included in pot file
  gettext("Randomness test for numeric variable...", domain="R-RcmdrPlugin.UCA")
  # Build dialog
  initializeDialog(title=gettext("Randomness test for numeric variable", domain="R-RcmdrPlugin.UCA"))
  variablesBox <- variableListBox(top, Numeric(), selectmode="single", initialSelection=NULL, title=gettextRcmdr("Variable (pick one)"))
	onOK <- function(){
		x <- getSelection(variablesBox)
		if (length(x) == 0) {
                    errorCondition(recall=randomnessNTest, message=gettextRcmdr("No variables were selected."))
                    return()
                }
		closeDialog()
                # Apply test
                doItAndPrint(paste("with(", ActiveDataSet(), ", numeric.runs.test(", x, "))", sep = ""))
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="runs.test", reset = "randomnessNTest", apply = "randomnessNTest")
	tkgrid(getFrame(variablesBox), sticky="nw")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=6, columns=1)
}

numeric.runs.test <- function(...) randtest.runs.test(...)
twolevelfactor.runs.test <- function(...) tseries.runs.test(...)
