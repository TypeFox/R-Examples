######## Basic Simple Random sample window
######## Last edited 7 November 2012
######## Author Susie Jentoft

SRS.notebook<-function() {
	initializeDialog(title=gettextRcmdr("One-Stage Simple Random Sample"))
	.activeDataSet <- ActiveDataSet()
	sizedefault <- ifelse(nrow(get(.activeDataSet))>=200, 100, round(nrow(get(.activeDataSet))/2))
    	sizeVar <- tclVar(sizedefault)
    	sizeEntry <- tkentry(top, width="6", textvariable=sizeVar)

	onOK <- function(){
		closeDialog()
        		nsize <- as.numeric(tclvalue(sizeVar))
		nvar <- paste("Sample_", format(Sys.time(), "%d%m%Y"), "_", format(Sys.time(), "%H%M"), sep="")
        		if (is.na(nsize)){
            			errorCondition(recall=SRS.notebook, message="Sample size must be a number.")
            			return()
            		}
		command1 <- paste(nvar, " <- sample(x = nrow(", .activeDataSet, "), size = ", nsize, ", replace = FALSE, prob = NULL)", sep="")
		command2 <- paste(nvar, " <- getdata(", .activeDataSet, ", ", nvar, ")", sep="")
		doItAndPrint(command1)
		doItAndPrint(command2)
		doItAndPrint(paste("View(", nvar, ")", sep=""))
		return()	
	}
	OKCancelHelp(helpSubject="sample", reset="SRS.notebook")
	tkgrid(tklabel(top, text="Please choose a sample size "), sticky="w", columnspan=3)
	tkgrid(tklabel(top, text=""))
	tkgrid(tklabel(top, text=gettextRcmdr("Sample size: ")), sizeEntry, sticky="w")
	tkgrid(tklabel(top, text="")) 
   	tkgrid(buttonsFrame, sticky="w", columnspan=8)
	dialogSuffix(rows=7, columns=2, focus=sizeEntry)
}

