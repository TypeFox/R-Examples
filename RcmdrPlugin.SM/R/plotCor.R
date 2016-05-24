
plotCor <- function () {
	Library("ade4")
	
defaults <- list(initial.variables = NULL,initial.mean="0") 

	dialog.values <- getDialog("plotCor", defaults)
	

	initializeDialog(title = gettextRcmdr("Analyse en Composantes Principales"))
	variablesBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Select variables (three or more)"), 
			selectmode = "multiple", initialSelection = varPosn (dialog.values$initial.variables, "numeric"))

checkBoxes(frame="checkBoxFrame", boxes=c("mean"), 
               initialValues=c(dialog.values$initial.mean), 
               labels=gettextRcmdr(c("ACP sur matrice de covariances")))

	onOK <- function() {
		variables <- getSelection(variablesBox)

meanVar <- tclvalue(meanVariable)
putDialog ("plotCor", list (initial.variables = variables, initial.mean=meanVar))

		closeDialog()
		if (length(variables) < 3) {
			errorCondition(recall =plotCor, message = gettextRcmdr("Fewer than 3 variable selected."))
			return()
		}
				
		.activeDataSet <- ActiveDataSet() 


		listvar<-paste(variables,collapse="\",\"")


if (meanVar == 0){
			command <- paste("barplot(dudi.pca(na.omit(",.activeDataSet,"[c(\"",listvar,"\")]),scannf=FALSE,nf=2)$eig)", sep = "")
			logger(command)
			justDoIt(command)

			command <- "windows()"
			logger(command)
			justDoIt(command)


			command <- paste("s.corcircle(dudi.pca(na.omit(",.activeDataSet,"[c(\"",listvar,"\")]),scannf=FALSE,nf=2)$co)", sep = "")
			logger(command)
			justDoIt(command)
}
else{

command <- paste("barplot(dudi.pca(na.omit(",.activeDataSet,"[c(\"",listvar,"\")]),scannf=FALSE,nf=2,scale=FALSE)$eig)", sep = "")
			logger(command)
			justDoIt(command)

			command <- "windows()"
			logger(command)
			justDoIt(command)


			command <- paste("s.arrow(dudi.pca(na.omit(",.activeDataSet,"[c(\"",listvar,"\")]),scannf=FALSE,nf=2,scale=FALSE)$co)", sep = "")
			logger(command)
			justDoIt(command)
}
		
		activateMenus()
		tkfocus(CommanderWindow())
	}
	groupsBox(scatterPlot)
	OKCancelHelp(helpSubject = "dotchart", reset = "plotCor")
	tkgrid(getFrame(variablesBox), sticky = "nw")
	
tkgrid(checkBoxFrame, sticky="w")

	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	dialogSuffix(rows = 6, columns = 2)
}

