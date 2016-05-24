# Package: FFD
# Internal functions for FFD_GUI()
# 
# Author: Ian Kopacka
###############################################################################

importDataFun <- function(inputDataVar, herdSizeColVar, riskGroupColVar){
	fileName <- tclvalue(tkgetOpenFile(filetypes = 
		"{{CSV} {.csv}}")) 	
    if (nchar(fileName) > 0){
        tclvalue(inputDataVar) <- fileName	
	    tclvalue(herdSizeColVar) <- ""
        tclvalue(riskGroupColVar) <- ""
    }
}

chooseColumnName <- function(inputDataVar, columnVar, dispText){
	fileName <- tclvalue(inputDataVar)	
	if (nchar(fileName) > 0){	
		dataFileHead <- try(suppressWarnings(read.csv2(file = fileName, 
			header = TRUE, nrows = 1)), silent = TRUE)
		if (class(dataFileHead) == "data.frame"){
			dataFileHead <- read.csv2(file = fileName, header = TRUE, nrows = 1)
			colNamesVec <- names(dataFileHead)
			
			## Select appropriate column:
			mainColumns <- tktoplevel(borderwidth = 10)
			tkwm.title(mainColumns, paste("FFD - ", dispText, sep = ""))
			tkwm.resizable(mainColumns, FALSE, FALSE)   # Window not resizable
			## Text frame:
			frameText <- tkframe(mainColumns, relief = "groove", borderwidth = 0)			
			tkgrid(tklabel(frameText, 
				text = paste("Please specify the column containing the ", 
				dispText, ":", sep = "")))
			tkgrid(tklabel(frameText, text = ""))
			columnsBox <- tkwidget(frameText, "ComboBox", editable = TRUE,
				values = colNamesVec)
			tkgrid(columnsBox)
			tkgrid(tklabel(frameText, text = ""))
			
			## Buttons frame:
			frameButtons <- tkframe(mainColumns, relief = "groove", borderwidth = 0)			
			ok.but <- tkbutton(frameButtons, text = " OK ", 
				command = function() onOkChooseCol(colNamesVec, columnsBox, 
				columnVar, mainColumns), padx = 15, pady = 2)
			cancelHS.but <- tkbutton(frameButtons, text = " Cancel ", 
				command = function() tkdestroy(mainColumns), padx = 15, pady = 2)
			tkgrid(ok.but, tklabel(frameButtons, text = " "), cancelHS.but)	
			
			## Build frames:
			tkgrid(frameText)
			tkgrid(frameButtons)
			
		} else {
			tkmessageBox(message = "Error! Could not open farm data file: No such file or directory.", 
				icon = "error", type = "ok")			
		}
	} else {
		tkmessageBox(message = "Error! No data file selected.", 
			icon = "error", type = "ok")
	}
}

onOkChooseCol <- function(colNamesVec, columnsBox, columnVar, mainColumns){
	columnName <- colNamesVec[as.numeric(tclvalue(tcl(columnsBox,"getvalue")))+1]
	tclvalue(columnVar) <- columnName	
	tkdestroy(mainColumns)
}

###############################################################################

## Function to assign the risk group parameters. 
## Is called by clicking on "Set risk group parameters" on the "Parameters"
## tab of the main window.
setRGParameters <- function(inputDataVar, riskGroupColVar, riskGroupTabVar,
	nRiskGroups, mainEnvir){
	
	## Check if window is already open. If window is already open do nothing:
    ID <- get("setRGWindow_ID", envir = mainEnvir)
#	ID <- get("setRGWindow_ID", envir = .GlobalEnv)
	if (exists(ID, envir = mainEnvir, inherits = FALSE)) return(NULL)
	
	## Read in population data:		
	if (file.exists(tclvalue(inputDataVar))){
	    populationData <- read.csv2(file = tclvalue(inputDataVar))
	} else {
		tkmessageBox(message = 
			"Error! Could not open farm data file: No such file or directory.", 
			icon = "error", type = "ok")
	    return(NULL)		
	}    
	
	## Extract riskGroupVec:
	riskGroupColName <- tclvalue(riskGroupColVar)
	if (riskGroupColName == ""){
		riskGroupVec <- character()
		tkmessageBox(message = "Error! No risk groups specified.", 
			icon = "error", type = "ok")
	    return(NULL)
	} else {
		if (riskGroupColName %in% names(populationData)){
	        riskGroupVec <- as.character(populationData[,riskGroupColName])
	    } else {
		    tkmessageBox(message = 
			    "Error! Could not import risk groups: Incorrect column name.", 
			    icon = "error", type = "ok")
	        return(NULL)
	    }
	}
	
	setRGWindow <- tktoplevel(padx = 10, pady = 10)	
	ID <- .Tk.ID(setRGWindow)
	## Save id of window to global variable (for the check if the window is 
	## already open or not):
	assign("setRGWindow_ID", ID, envir = mainEnvir)
#    assign("setRGWindow_ID", ID, envir = .GlobalEnv)
	tkwm.title(setRGWindow, "FFD - Risk Group Data")
	tkwm.resizable(setRGWindow, FALSE, FALSE)
	rgTop <- tkframe(setRGWindow, padx = 5, pady = 5)
	rgBottom <- tkframe(setRGWindow, padx = 5, pady = 5)

    ## Assign Values to Table:
	rgNames <- sort(unique(riskGroupVec))
	nRiskGroupsVal <- length(rgNames)
	nRiskGroupsValOld <- as.numeric(tclvalue(nRiskGroups))
	tclvalue(nRiskGroups) <- as.character(nRiskGroupsVal) 
	namesRows <- rgNames
	namesCols <- c("Risk groups:", "Risk_values", "Fixed_sample_sizes", 
		"Sample_probabilty")
	maxWidth <- max(c(nchar(namesRows),nchar(namesCols)))
	
	## Delete old values:	
	if (nRiskGroupsValOld > nRiskGroupsVal){
		for (ii in (nRiskGroupsVal+1):nRiskGroupsValOld){
		    for (jj in seq(along = namesRows)){
				riskGroupTabVar[[ii,jj-1]] <- character()
			}			
		}		
	}
	## Set new values:
    for(ii in seq(along = namesCols)) riskGroupTabVar[[0,ii-1]] <- namesCols[ii]
    for(jj in seq(along = namesRows)) riskGroupTabVar[[jj,0]] <- namesRows[jj]
    tableRg <- tkwidget(rgTop, "table", resizeborders = "none",
			selectmode = "extended", multiline = "0", rowseparator = "\n",
			colseparator = "\t", variable = riskGroupTabVar, rows = length(namesRows) + 1, 
			cols = length(namesCols),
			titlerows = 1, titlecols = 1, colwidth = maxWidth, 
			background = "white")	
	
	## Buttons:
	save.but <- tkbutton(rgBottom, text = " Save and close ", 
		padx = 15, pady = 2, command = function() tkdestroy(setRGWindow))
    reset.but <- tkbutton(rgBottom, text = " Reset ", 
		padx = 15, pady = 2, command = function() resetArray(riskGroupTabVar, 
		length(namesRows) + 1, length(namesCols)))

    ## Create Frames and buttons:
	tkgrid(tableRg)
	tkgrid(save.but, tklabel(rgBottom, text = " "), reset.but)
	tkgrid(rgTop)
	tkgrid(rgBottom)
}

resetArray <- function(arrayVar, nRows, nCols){
    for(ii in 1:(nRows-1)){
		for (jj in 1:(nCols-1)){
		    arrayVar[[ii,jj]] <- character()	
		}		
	} 	
} 
		

###############################################################################

## Summaries speichern:
saveText <- function(textVec){
	fileName <- tclvalue(tkgetSaveFile(filetypes = 
	    "{{Text Files} {.txt}} {{All files} *}", defaultextension = ".txt"))
    if (nchar(fileName) > 0){
        outFile <- file(fileName, "w+")
	    for (ii in seq(along = textVec)) cat(paste(textVec[ii], "\n"), file = outFile)
	    close(outFile)
    }
}

###############################################################################

## Graphiken speichern:
savePlot <- function(mySamplingSummary){
	fileName <- tclvalue(tkgetSaveFile(filetypes = 
 		"{{JPG} {.jpg}} {{PNG} {.png}} {{PDF} {.pdf}} {{All files} *}",
              defaultextension = ".png"))
    if (nchar(fileName) > 0){
        ## Extension:
	    fileType <- substr(fileName, nchar(fileName)-2, nchar(fileName))
	    ## If-Abfrage bzgl. File-Type:
		switch(fileType,
			"jpg" = jpeg(filename = fileName, width = 800, height = 800, 
				pointsize = 14),
			"png" = png(filename = fileName, width = 800, height = 800, 
				pointsize = 14),
			"pdf" = pdf(file = fileName, pointsize = 10), 
			png(filename = fileName, width = 800, height = 800))
		plot(mySamplingSummary)
		dev.off()
	}
}

###############################################################################

## Create mySurvey object:
buildMySurvey <- function(infoList){
	## Required Parameters:
	#######################	
	## population data:		
	if (file.exists(tclvalue(infoList$inputDataVar))){
	    populationData <- read.csv2(file = tclvalue(infoList$inputDataVar))
	} else {
		tkmessageBox(message = 
			"Error! Could not open farm data file: No such file or directory.", 
			icon = "error", type = "ok")
	    return(NULL)		
	}    
	
	## nAnimalVec:
	if (tclvalue(infoList$herdSizeColVar) %in% names(populationData)){
	    nAnimalVec <- populationData[,tclvalue(infoList$herdSizeColVar)]
	} else {
		tkmessageBox(message = 
			"Error! Could not import herd sizes: Incorrect column name.", 
			icon = "error", type = "ok")
	    return(NULL)
	}
	
	## riskGroupVec:
	riskGroupColName <- tclvalue(infoList$riskGroupColVar)
	if (riskGroupColName == ""){
		riskGroupVec <- character()
		nRiskGroups <- 0
	} else {
		if (riskGroupColName %in% names(populationData)){
	        riskGroupVec <- as.character(populationData[,riskGroupColName])
			nRiskGroups <- length(unique(riskGroupVec))
	    } else {
		    tkmessageBox(message = 
			    "Error! Could not import risk groups: Incorrect column name.", 
			    icon = "error", type = "ok")
	        return(NULL)
	    }
	}

	## Risk group data:
	if (nRiskGroups > 0){
		riskGroup <- rep(NA, nRiskGroups)
		riskValues <- rep(NA, nRiskGroups)
		for (ii in 1:nRiskGroups){
		    ## riskValueData: column riskGroup	
			if (length(infoList$riskGroupTabVar[[ii,0]]) > 0){
				riskGroup[ii] <- tclvalue(infoList$riskGroupTabVar[[ii,0]])
			} else {
				riskGroup[ii] <- NA
			}
		    ## riskValueData: column riskValues	
			if (length(infoList$riskGroupTabVar[[ii,1]]) > 0){
				riskValues[ii] <- tclvalue(infoList$riskGroupTabVar[[ii,1]])
			} else {
				riskValues[ii] <- NA
			}				
		}		
		riskGroup[riskGroup == ""] <- NA
		riskValues[riskValues == ""] <- NA
		if (any(is.na(c(riskGroup, riskValues)))){
			riskValueData <- data.frame()
			tkmessageBox(message = 
			    "Warning! Could not import risk group parameters.", 
			    icon = "warning", type = "ok")
		} else {
			riskValueData <- data.frame(riskGroup = as.character(riskGroup),
		        riskValues = as.numeric(riskValues))
	    }
	} else {
		riskValueData <- data.frame()
	}	
	
	## pi: 
	designPrevalence <- as.numeric(tclvalue(infoList$pi))/100
	## alpha:
	alphaVal <- as.numeric(tclvalue(infoList$alpha))
	## pi_IH:
	intraHerdPrevalence <- as.numeric(tclvalue(infoList$piIH))/100
	## Se:
	diagSensitivity <- as.numeric(tclvalue(infoList$Se))/100
	
	## Optional Parameters:
	#######################
	costHerdVal <- tclvalue(infoList$costHerd)
	if (costHerdVal == ""){
		costHerdVal <- numeric()
	} else {
		costHerdVal <- as.numeric(costHerdVal)
	}
	costAnimalVal <- tclvalue(infoList$costAnimal)
	if (costAnimalVal == ""){
		costAnimalVal <- numeric()
	} else {
		costAnimalVal <- as.numeric(costAnimalVal)
	}	
	
	## Create surveyData-Object:
	mySurvey <- try(surveyData(nAnimalVec = nAnimalVec,
		riskGroupVec = riskGroupVec, riskValueData = riskValueData,
		populationData = populationData, designPrevalence = designPrevalence,
		alpha = alphaVal, intraHerdPrevalence = intraHerdPrevalence,
		diagSensitivity = diagSensitivity, costHerd = costHerdVal, 
		costAnimal = costAnimalVal), silent = TRUE)

    if (class(mySurvey) == "try-error") {
	    tkmessageBox(message = as.character(mySurvey), icon = "error", type = "ok")
	    return(NULL)
	} else {
		return(mySurvey)
	}	
}

###############################################################################

importRiskBasedSamplingParameters <- function(mySurvey, infoList){
	nRiskGroups <- dim(mySurvey@riskValueData)[1]
	nSampleFixVec <- rep(NA, nRiskGroups)
	probVec <- rep(NA, nRiskGroups)
	for (ii in 1:nRiskGroups){
		## riskValueData: column riskGroup	
		if (length(infoList$riskGroupTabVar[[ii,2]]) > 0){
			nSampleFixVec[ii] <- tclvalue(infoList$riskGroupTabVar[[ii,2]])
		} else {
			nSampleFixVec[ii] <- NA
		}
		## riskValueData: column riskValues	
		if (length(infoList$riskGroupTabVar[[ii,3]]) > 0){
			probVec[ii] <- tclvalue(infoList$riskGroupTabVar[[ii,3]])
		} else {
			probVec[ii] <- NA
		}				
	}	
	nSampleFixVec[nSampleFixVec == ""] <- NA
	probVec[probVec == ""] <- NA
	nSampleFixVec <- as.numeric(nSampleFixVec)
	probVec <- as.numeric(probVec)
	
	## Is the model over-specified, i.e., nSampleFix AND prob specified for a 
	## risk group:
	if (any(!is.na(nSampleFixVec) & !is.na(probVec))){
	    tkmessageBox(message = paste("Warning! Fixed sample size and sample probability\n",
			"both specified for a risk group.\nIgnoring sample probability.", sep = ""), 
			icon = "warning", type = "ok")	
	}
	
	## Sample size for all but one risk group is fixed:
	if (sum(is.na(nSampleFixVec)) == 1){
		return(list(nSampleFixVec = nSampleFixVec, probVec = 1))
	}
	
	## Check for missing data for a risk group:
	if (any(is.na(nSampleFixVec) & is.na(probVec))){
		tkmessageBox(message = paste("Warning! Risk group parameters incomplete.\n",
			"Ignoring risk groups.", sep = ""), 
			icon = "warning", type = "ok")
		return(list(nSampleFixVec = NULL, probVec = NULL))
	}
	
	probVec <- probVec[is.na(nSampleFixVec)]		
	## Return value:
	return(list(nSampleFixVec = nSampleFixVec, probVec = probVec))
}


###############################################################################

## Function: Tab Sample Size, Calculate Sample size
OnCalcSS <- function(infoList, sampSizeVar, main){
	## Create mySurvey object:
	mySurvey <- buildMySurvey(infoList)     
	
	if (!is.null(mySurvey)){	
	## Which sampling strategy:		
	boxValue <- infoList$sampStratVec[as.numeric(tclvalue(tcl(infoList$sampStratBox,
		"getvalue")))+1]	
	if (length(boxValue) > 0){
		## With or without risk groups:
		if (dim(mySurvey@riskValueData)[1] > 0){
			useRiskGroups <- tkmessageBox(message = "Should risk groups be considered?",
					icon = "question", type = "yesnocancel", default = "yes")
			if (tclvalue(useRiskGroups) == "cancel") return(NULL)
			if (tclvalue(useRiskGroups) == "yes"){
				riskDataList <- importRiskBasedSamplingParameters(mySurvey, infoList)
				nSampleFixVec <- riskDataList$nSampleFixVec
				probVec <- riskDataList$probVec
#				nRiskGroups <- dim(mySurvey@riskValueData)[1]
#				nSampleFixVec <- rep(NA, nRiskGroups)
#				probVec <- rep(NA, nRiskGroups)
#				for (ii in 1:nRiskGroups){
#					## riskValueData: column riskGroup	
#					if (length(infoList$riskGroupTabVar[[ii,2]]) > 0){
#						nSampleFixVec[ii] <- tclvalue(infoList$riskGroupTabVar[[ii,2]])
#					} else {
#						nSampleFixVec[ii] <- NA
#					}
#					## riskValueData: column riskValues	
#					if (length(infoList$riskGroupTabVar[[ii,3]]) > 0){
#						probVec[ii] <- tclvalue(infoList$riskGroupTabVar[[ii,3]])
#					} else {
#						probVec[ii] <- NA
#					}				
#				}	
#				nSampleFixVec[nSampleFixVec == ""] <- NA
#				probVec[probVec == ""] <- NA
#				class(nSampleFixVec) <- "numeric"	
#				probVec <- probVec[is.na(nSampleFixVec)]
#				probVec <- as.numeric(probVec)		
#				if (any(is.na(probVec))){
#				    nSampleFixVec <- NULL
#				    probVec <- NULL
#					tkmessageBox(message = 
#			            "Warning! Could not import risk group parameters.", 
#			            icon = "warning", type = "ok")
#				}				
			} else {
				nSampleFixVec <- NULL
				probVec <- NULL
			}			
		} else {
			nSampleFixVec <- NULL
			probVec <- NULL
		}	
		
		if (boxValue == "limited sampling"){
			kVal <- suppressWarnings(as.numeric(tclvalue(sampSizeVar)))
			if (!is.na(kVal)){
			    ## Compute sample size:
			    tkconfigure(main,cursor = "watch")
				mySampling <- try(ltdSampling(survey.Data = mySurvey,
				    sampleSizeLtd = kVal, nSampleFixVec = nSampleFixVec,
					probVec = probVec), silent = TRUE)		
			    tkconfigure(main,cursor = "arrow")
			    tttitle <- "FFD - GUI: Limited Sampling"
			} else {
			    tkmessageBox(message = 
					"Error! Invalid value for sample limit.", 
			        icon = "error", type = "ok")
				return(NULL)
			}
		} else {
			herdSensVal <- as.numeric(tclvalue(sampSizeVar))/100
			if (!is.na(herdSensVal)){
			    ## Compute sample size:
			    tkconfigure(main,cursor = "watch")
				mySampling <- try(indSampling(survey.Data = mySurvey,
				    herdSensitivity = herdSensVal, nSampleFixVec = nSampleFixVec,
					probVec = probVec), silent = TRUE)
			    tkconfigure(main,cursor = "arrow")
			    tttitle <- "FFD - GUI: Individual Sampling"
			} else {
			    tkmessageBox(message = 
					"Error! Invalid value for herd sensitivity.", 
			        icon = "error", type = "ok")
				return(NULL)
			}
		}		
		
	    if (class(mySampling) == "try-error"){
			tkmessageBox(message = as.character(mySampling), 
			        icon = "error", type = "ok")
			return(NULL)
		}	
	
	## Build widget:
	tt  <- tktoplevel()
	tkwm.resizable(tt, FALSE, FALSE)   # Window not resizable
	tkwm.title(tt, tttitle)
	## Main Frame:
	##########################################################
	frameOverall <- tkframe(tt)
	
	## Text Frame:
	##########################################################
	frameText <- tkframe(frameOverall,padx = 10, pady = 10,
			relief = "groove", borderwidth = 1)
	## Scrollbar:
	xscr <- tkscrollbar(frameText, repeatinterval=5,orient="horizontal",
			command=function(...)tkxview(txt,...))
	yscr <- tkscrollbar(frameText, repeatinterval=5,
			command=function(...)tkyview(txt,...))
	txt <- tktext(frameText,bg="white",
			xscrollcommand=function(...)tkset(xscr,...),yscrollcommand = 
			function(...)tkset(yscr,...), wrap="none")
	tkgrid(txt,yscr)
	tkgrid(xscr)
	tkgrid.configure(yscr,sticky="ns")
	tkgrid.configure(xscr,sticky="ew")
	## Write text:
	outVec <- as.vector(summary(mySampling, "percent"))
	for (ii in seq(along = outVec)) tkinsert(txt, "end", paste(outVec[ii], 
						"\n"))
	## Make it non-editable:		
	tkconfigure(txt, state = "disabled")
	## For resizing:
	tkgrid.columnconfigure(tt,0,weight=1)
	tkgrid.rowconfigure(tt,0,weight=1)
	
	## Button Frame:
	##########################################################
	frameButton <- tkframe(frameOverall,padx = 10, pady = 10,
		relief = "groove", borderwidth = 0)
	close.but <- tkbutton(frameButton, text = " Close ", 
		command = function() tkdestroy(tt), padx = 15, pady = 2)
	save.but <- tkbutton(frameButton, text = " Save ", 
		command = function() saveText(outVec), padx = 15, pady = 2)	
	tkgrid(save.but, tklabel(frameButton,text = " "), close.but)	
	
	## Draw frames:
	##########################################################
	tkgrid(frameText)
	tkgrid(frameButton)
	tkgrid(frameOverall)
	tkfocus(txt)	
	} else {
		tkmessageBox(message = "Error! No sampling strategy selected.",icon="error",type="ok")
	}
	}
}
	
###############################################################################

OnSampleSS <- function(infoList, sampSizeVar, main){
	## Which sampling strategy:		
	boxValue <- infoList$sampStratVec[as.numeric(tclvalue(tcl(infoList$sampStratBox,
									"getvalue")))+1]
	
	## Is sampling strategy selected?
	if (length(boxValue) > 0){
		if (boxValue == "limited sampling"){
		    sampleText <- "using Limited Sampling"
		} else {
			sampleText <- "using Individual Sampling"
		}
		
		## Build widget:
		tt  <- tktoplevel(borderwidth = 15)
		tkwm.resizable(tt, FALSE, FALSE)   # Window not resizable
		fontBold <- tkfont.create(size = 8, weight = "bold")	
		tkwm.title(tt, "FFD - GUI: Sampling")
		
		## Main Frame:
		##########################################################
		frameOverall <- tkframe(tt)
		tkgrid(tklabel(frameOverall, text = " Draw sample from Population ", 
			font = fontBold))
		tkgrid(tklabel(frameOverall, text = sampleText, 
			font = fontBold))
		tkgrid(tklabel(frameOverall, text = ""))		
		
		## Text Frame:
		##########################################################
		frameText <- tkframe(frameOverall,padx = 10, pady = 10,
			relief = "groove", borderwidth = 1, width = 200)
		## Sample size:
		sampSizeVec <- c("fixed", "dynamic")
		sampSizeBox <- tkwidget(frameText, "ComboBox", editable = FALSE,
			values = sampSizeVec, width = "8")
		tkgrid(tklabel(frameText,text = "Sample size: "), 
			sampSizeBox, sticky = "e")
		
		## Seed:
		seedValue <- sample.int(n = 100000, size = 1)
		seedVar <- tclVar(as.character(seedValue))
		entry.seed <- tkentry(frameText, width = "11", textvariable = seedVar)
		tkgrid(tklabel(frameText,text = "Seed: "), 
			entry.seed, sticky = "e")  
		
		## Output Frame:
		##########################################################
		frameOutput <- tkframe(frameOverall,padx = 10, pady = 10,
			relief = "groove", borderwidth = 1, width = 200)
		## Variables:
		## Computed sample size:
	    compSampSizeVar <- tclVar("")
		compSampSizeLabel <- tklabel(frameOutput, 
			text = tclvalue(compSampSizeVar), width = "9")
		tkconfigure(compSampSizeLabel, textvariable = compSampSizeVar)	
		## Actual sample size:
		actSampSizeVar <- tclVar("")
		sampSizeLabel <- tklabel(frameOutput, 
			text = tclvalue(actSampSizeVar), width = "9")
		tkconfigure(sampSizeLabel, textvariable = actSampSizeVar)	
		## A-posteriori alpha error:
		alphaErrorVar <- tclVar("")
		alphaErrorLabel <- tklabel(frameOutput, text = tclvalue(alphaErrorVar), width = "9")
		tkconfigure(alphaErrorLabel, textvariable = alphaErrorVar)
		
		tkgrid(tklabel(frameOutput,text = "Computed sample size : "), 
			compSampSizeLabel, sticky = "e")
		tkgrid(tklabel(frameOutput,text = "Actual sample size : "), sampSizeLabel,
			sticky = "e")
		tkgrid(tklabel(frameOutput,text = "Alpha error: "), alphaErrorLabel,
			sticky = "e")
		
		## Button Frame:
		##########################################################
		frameButton <- tkframe(frameOverall,padx = 10, pady = 10,
			relief = "groove", borderwidth = 0)
		close.but <- tkbutton(frameButton, text = " Close ", 
			command = function() tkdestroy(tt), padx = 15, pady = 2)
		sample.but <- tkbutton(frameButton, text = " Sample ", 
			command = function() sampleFarms(tt, main, infoList, sampSizeVar, boxValue, 
			sampSizeVec, sampSizeBox, seedVar, compSampSizeVar, actSampSizeVar, 
			alphaErrorVar), 
            padx = 15, pady = 2)	
		tkgrid(sample.but, tklabel(frameButton,text = " "), close.but)	
		
		## Draw frames:
		##########################################################
		tkgrid(frameText, pady = 10)
		tkgrid(frameOutput, pady = 10)
		tkgrid(frameButton)
		tkgrid(frameOverall)
		tkfocus(tt)
		
	} else {
		## No sampling strategy selected:
		tkmessageBox(message = "Error! No sampling strategy selected.",
			icon = "error", type = "ok")		
	}	
}

###############################################################################

sampleFarms <- function(tt, main, infoList, sampSizeVar, boxValue, 
	sampSizeVec, sampSizeBox, seedVar, compSampSizeVar, 
	actSampSizeVar, alphaErrorVar){
	## Seed:
	seedValue <- suppressWarnings(as.integer(tclvalue(seedVar)))
	if (is.na(seedValue)){
		tkmessageBox(message = "Error! Seed must be an integer.", 
			icon = "error", type = "ok")
		return(NULL)
	}
	## Sample size strategy:
	sampSizeScheme <- sampSizeVec[as.numeric(tclvalue(tcl(sampSizeBox,"getvalue")))+1]
	## Is a sample size strategy selected?
	if (length(sampSizeScheme) > 0){	
		## Sampling:
		############
		
		# 1) Create mySurvey-Object:
		mySurvey <- buildMySurvey(infoList)
		## Did it work?
		if (!is.null(mySurvey)){				
			## 2) Choose: With or without risk groups:
			if (dim(mySurvey@riskValueData)[1] > 0){
				useRiskGroups <- tkmessageBox(message = "Should risk groups be considered?",
					icon = "question", type = "yesnocancel", default = "yes")
				if (tclvalue(useRiskGroups) == "cancel") return(NULL)
				if (tclvalue(useRiskGroups) == "yes"){
					riskDataList <- importRiskBasedSamplingParameters(mySurvey, infoList)
					nSampleFixVec <- riskDataList$nSampleFixVec
					probVec <- riskDataList$probVec
#				
#					nRiskGroups <- dim(mySurvey@riskValueData)[1]
#					nSampleFixVec <- rep(NA, nRiskGroups)
#					probVec <- rep(NA, nRiskGroups)
#					for (ii in 1:nRiskGroups){
#						## riskValueData: column riskGroup	
#						if (length(infoList$riskGroupTabVar[[ii,2]]) > 0){
#							nSampleFixVec[ii] <- tclvalue(infoList$riskGroupTabVar[[ii,2]])
#						} else {
#							nSampleFixVec[ii] <- NA
#						}
#						## riskValueData: column riskValues	
#						if (length(infoList$riskGroupTabVar[[ii,3]]) > 0){
#							probVec[ii] <- tclvalue(infoList$riskGroupTabVar[[ii,3]])
#						} else {
#							probVec[ii] <- NA
#						}				
#					}	
#					nSampleFixVec[nSampleFixVec == ""] <- NA
#					probVec[probVec == ""] <- NA
#					class(nSampleFixVec) <- "numeric"	
#					probVec <- probVec[is.na(nSampleFixVec)]
#					probVec <- as.numeric(probVec)		
#					if (any(is.na(probVec))){
#						nSampleFixVec <- NULL
#						probVec <- NULL
#						tkmessageBox(message = 
#							"Warning! Could not import risk group parameters.", 
#							icon = "warning", type = "ok")
#					}				
				} else {
					nSampleFixVec <- NULL
					probVec <- NULL
				}			
			} else {
				nSampleFixVec <- NULL
				probVec <- NULL
			}
			
			## Save location:
			fileName <- tclvalue(tkgetSaveFile(filetypes = 
									"{{CSV} {.csv}} {{All files} *}", defaultextension = ".csv"))
			if (nchar(fileName) > 0){
				
				# 3) Create mySampling-Object:
				if (boxValue == "limited sampling"){
					kVal <- suppressWarnings(as.numeric(tclvalue(sampSizeVar)))					
					if (!is.na(kVal)){
						## Compute sample size:
						tkconfigure(main,cursor = "watch")						
						mySampling <- try(suppressWarnings(ltdSampling(survey.Data = mySurvey,
							sampleSizeLtd = kVal, nSampleFixVec = nSampleFixVec,
							probVec = probVec)), silent = TRUE)	
						tkconfigure(main,cursor = "arrow")
						sampleText <- "using Limited Sampling"
					} else {
						tkmessageBox(message = 
							"Error! Invalid value for sample limit.", 
							icon = "error", type = "ok")
						return(NULL)
					}
				} else {
					herdSensVal <- as.numeric(tclvalue(sampSizeVar))/100
					if (!is.na(herdSensVal)){
						## Compute sample size:
						tkconfigure(main,cursor = "watch")
						try(suppressWarnings(mySampling <- indSampling(survey.Data = mySurvey,
							herdSensitivity = herdSensVal, nSampleFixVec = nSampleFixVec,
							probVec = probVec)), silent = TRUE)	
						tkconfigure(main,cursor = "arrow")
						sampleText <- "using Individual Sampling"
					} else {
						tkmessageBox(message = 
							"Error! Invalid value for herd sensitivity.", 
							icon = "error", type = "ok")
						return(NULL)
					}
				}
				
				# 4) Did that work?
				if (class(mySampling) == "try-error"){
					tkmessageBox(message = as.character(mySampling), 
						icon = "error", type = "ok")
					return(NULL)
				}
				
				# 5) Sample:
				set.seed(seedValue)		
				sampleList <- sample(x = mySampling, size = sampSizeScheme)
				saveCheck <- suppressWarnings(try(write.csv2(sampleList$sample, 
					file = fileName, row.names = FALSE), silent = TRUE))
                if (class(saveCheck) == "NULL"){
					## Saving worked:
					tclvalue(compSampSizeVar) <- as.character(mySampling@nHerds)
					tclvalue(actSampSizeVar) <- as.character(length(sampleList$indexSample))
					tclvalue(alphaErrorVar) <- paste(as.character(round(sampleList$aPostAlpha*100,3)),
						"%", sep = " ")
			    } else {
					## Saving did not work:
					tkmessageBox(message = paste("Error! Could not write to file:\n",
						fileName, sep = ""), 
			            icon = "error", type = "ok")					
				}
			} else {
				## It did not work:
				tkdestroy(tt)
			}
		}
	} else {
		## No sample size strategy selected:
		tkmessageBox(message = "Error! No sample size strategy selected.", 
			icon = "error", type = "ok")
	} 
	return(NULL)
}


###############################################################################

# Hilfsfunktion:
plotWidget <- function(mySamplingSummary){
	plotFunction <- function()
	{
		params <- par(bg = "white")
		plot(mySamplingSummary)
		par(params)
	}
	plotWin <- tktoplevel()
	tkwm.resizable(plotWin, FALSE, FALSE)   # Window not resizable
	tkwm.title(plotWin, "Diagnostic Plot")
	## Plot frame:
	framePlot <- tkframe(plotWin, padx = 10, pady = 10,
		relief = "groove", borderwidth = 0)	
	img <- tkrplot(framePlot, fun = plotFunction, 
		hscale = 1.7, vscale = 1.7)
	tkgrid(img)
	## Buttons frame:
	frameButton <- tkframe(plotWin, padx = 10, pady = 10,
		relief = "groove", borderwidth = 0)	
	
	close.but <- tkbutton(frameButton, text=" Close ", 
			command = function() tkdestroy(plotWin), padx = 15, pady = 2)
	save.but <- tkbutton(frameButton, text=" Save ", 
		command = function() savePlot(mySamplingSummary), padx = 15, pady = 2)
	#tkgrid(close.but)
	tkgrid(save.but, tklabel(frameButton,text=" "), close.but)
	tkgrid(framePlot)
	tkgrid(frameButton)	
}
	
## Function: Tab Sample Size, Calculate Sample size
OnCalcDiag <- function(infoList, diagVar, main){
	## Create mySurvey object:
	mySurvey <- buildMySurvey(infoList)  
	
	if (!is.null(mySurvey)){
	## Which sampling strategy:		
	boxValue <- infoList$sampStratVec[as.numeric(tclvalue(tcl(infoList$sampStratBox,"getvalue")))+1]	
	if (length(boxValue) > 0){
		## With or without risk groups:
		if (dim(mySurvey@riskValueData)[1] > 0){
			useRiskGroups <- tkmessageBox(message = "Should risk groups be considered?",
					icon = "question", type = "yesnocancel", default = "yes")
			if (tclvalue(useRiskGroups) == "cancel") return(NULL)
			if (tclvalue(useRiskGroups) == "yes"){
				riskDataList <- importRiskBasedSamplingParameters(mySurvey, infoList)
				nSampleFixVec <- riskDataList$nSampleFixVec
				probVec <- riskDataList$probVec
#				nRiskGroups <- dim(mySurvey@riskValueData)[1]
#				nSampleFixVec <- rep(NA, nRiskGroups)
#				probVec <- rep(NA, nRiskGroups)
#				for (ii in 1:nRiskGroups){
#					## riskValueData: column riskGroup	
#					if (length(infoList$riskGroupTabVar[[ii,2]]) > 0){
#						nSampleFixVec[ii] <- tclvalue(infoList$riskGroupTabVar[[ii,2]])
#					} else {
#						nSampleFixVec[ii] <- NA
#					}
#					## riskValueData: column riskValues	
#					if (length(infoList$riskGroupTabVar[[ii,3]]) > 0){
#						probVec[ii] <- tclvalue(infoList$riskGroupTabVar[[ii,3]])
#					} else {
#						probVec[ii] <- NA
#					}				
#				}	
#				nSampleFixVec[nSampleFixVec == ""] <- NA
#				probVec[probVec == ""] <- NA
#				class(nSampleFixVec) <- "numeric"	
#				probVec <- probVec[is.na(nSampleFixVec)]
#				probVec <- as.numeric(probVec)		
#				if (any(is.na(probVec))){
#				    nSampleFixVec <- NULL
#				    probVec <- NULL
#					tkmessageBox(message = 
#			            "Warning! Could not import risk group parameters.", 
#			            icon = "warning", type = "ok")
#				}				
			} else {
				nSampleFixVec <- NULL
				probVec <- NULL
			}			
		} else {
			nSampleFixVec <- NULL
			probVec <- NULL
		}
		
		
		if (boxValue == "limited sampling"){
			sampleSizeLtdMax <- suppressWarnings(as.numeric(tclvalue(diagVar)))			
			if (!is.na(sampleSizeLtdMax)){
                ## Compute sample size:
			    ## Wait cursor: Cursor=Sanduhr waehrend gerechnet wird
			    tkconfigure(main, cursor = "watch")
			    mySamplingSummary <- 
					try(suppressWarnings(ltdSamplingSummary(survey.Data = mySurvey,
					sampleSizeLtdMax = sampleSizeLtdMax,
					nSampleFixVec = nSampleFixVec,
					probVec = probVec)), silent = TRUE)
			    tkconfigure(main,cursor = "arrow")
			    tttitle <- "FFD - GUI: Limited Sampling"
			} else {
			    tkmessageBox(message = 
					"Error! Invalid value for max. sample limit.", 
			        icon = "error", type = "ok")
				return(NULL)
			}			
		} else {
			stepSize <- suppressWarnings(as.numeric(tclvalue(diagVar))/100)
			if (!is.na(stepSize)){
                ## Compute sample size:
			    ## Wait cursor: Cursor=Sanduhr waehrend gerechnet wird
			    tkconfigure(main,cursor = "watch")
			    mySamplingSummary <- 
					try(suppressWarnings(indSamplingSummary(survey.Data = mySurvey,
					stepSize = stepSize,
					nSampleFixVec = nSampleFixVec,
					probVec = probVec)), silent = TRUE)	
			    tkconfigure(main,cursor = "arrow")					
			    tttitle <- "FFD - GUI: Individual Sampling"
			} else {
			    tkmessageBox(message = 
					"Error! Invalid value for herd sens. step size.", 
			        icon = "error", type = "ok")
				return(NULL)
			}			
		}
	if (class(mySamplingSummary) == "try-error"){
		tkmessageBox(message = 
			as.character(mySamplingSummary), 
			icon = "error", type = "ok")
		return(NULL)
	}
		
	outVec <- as.vector(summary(mySamplingSummary, "percent"))
	## Build widget:
	tt  <- tktoplevel()
	tkwm.resizable(tt, FALSE, FALSE)   # Window not resizable
	tkwm.title(tt, tttitle)
	## Main Frame:
	##########################################################
	frameOverall <- tkframe(tt)
	
	## Text Frame:
	##########################################################
	frameText <- tkframe(frameOverall,padx = 10, pady = 10,
			relief = "groove", borderwidth = 1)
	## Scrollbar:
	xscr <- tkscrollbar(frameText, repeatinterval = 5, orient = "horizontal",
			command = function(...)tkxview(txt,...))
	yscr <- tkscrollbar(frameText, repeatinterval = 5,
			command = function(...)tkyview(txt,...))
	txt <- tktext(frameText, bg = "white",
			xscrollcommand = function(...)tkset(xscr,...), yscrollcommand = 
			function(...)tkset(yscr,...), wrap = "none")
	tkgrid(txt,yscr)
	tkgrid(xscr)
	tkgrid.configure(yscr,sticky="ns")
	tkgrid.configure(xscr,sticky="ew")
	## Write text:
	for (ii in seq(along = outVec)) tkinsert(txt, "end", paste(outVec[ii], 
						"\n"))
	## Make it non-editable:		
	tkconfigure(txt, state = "disabled")
	## For resizing:
	tkgrid.columnconfigure(tt, 0, weight = 1)
	tkgrid.rowconfigure(tt, 0, weight = 1)
	
	
	## Button Frame:
	##########################################################
	frameButton <- tkframe(frameOverall,padx = 10, pady = 10,
			relief = "groove", borderwidth = 0)
	close.but <- tkbutton(frameButton, text=" Close ", 
			command = function() tkdestroy(tt), padx = 15, pady = 2)
	save.but <- tkbutton(frameButton, text=" Save ", 
			command = function() saveText(outVec), padx = 15, pady = 2)
	plot.but <- tkbutton(frameButton, text=" Plot ", 
			command = function() plotWidget(mySamplingSummary), padx = 15, 
			pady = 2)
	tkgrid(save.but, tklabel(frameButton,text=" "), close.but,
			tklabel(frameButton,text=" "), plot.but)	
	
	## Draw frames:
	##########################################################
	tkgrid(frameText)
	tkgrid(frameButton)
	tkgrid(frameOverall)
	tkfocus(txt)	
	} else {
		tkmessageBox(message = "Error! No sampling strategy selected.",
			icon = "error", type = "ok")
	}
	}
}

###############################################################################
FFDmanual <- function(){
    pdfManual <- system.file("doc/FFD-intro.pdf", package = "FFD")
    if (.Platform$OS.type == "windows"){
        shell.exec(pdfManual)
    } else {
    system(paste(shQuote(getOption("pdfviewer")), shQuote(pdfManual)), wait = FALSE)
    }
}

###############################################################################
saveFFD <- function(inputDataVar, herdSizeColVar, riskGroupColVar, pi, alpha, 
	piIH, Se, costHerd, costAnimal, sampSizeVar, diagVar, sampStratBox, 
	nRiskGroups, riskGroupTabVar){

    ## Risk group data:
	###################
	nRiskGroupsVal <- tclvalue(nRiskGroups)
	if (as.numeric(nRiskGroupsVal) > 0){	
		riskGroupTabVarVal <- matrix(NA, as.numeric(nRiskGroupsVal)+1, 4)
		for (ii in 1:(as.numeric(nRiskGroupsVal)+1)){
		    for (jj in 1:4){
			    if (length(riskGroupTabVar[[ii-1,jj-1]]) > 0){
					riskGroupTabVarVal[ii,jj] <- tclvalue(riskGroupTabVar[[ii-1,jj-1]])
				} 				
		    }
	    }		
	} else {
		riskGroupTabVarVal <- matrix(numeric(), 0, 0)		
	}	

	## Collect data:
	################
    FFD_InfoList <- list(inputDataVar = tclvalue(inputDataVar),
		herdSizeColVar = tclvalue(herdSizeColVar),
		riskGroupColVar = tclvalue(riskGroupColVar),
		pi = tclvalue(pi),
		alpha = tclvalue(alpha),
		piIH = tclvalue(piIH),
		Se = tclvalue(Se),
		costHerd = tclvalue(costHerd),
		costAnimal = tclvalue(costAnimal),
		sampSizeVar = tclvalue(sampSizeVar),
		diagVar = tclvalue(diagVar),
		sampStratBox = tclvalue(tcl(sampStratBox, "getvalue")),
		nRiskGroups = nRiskGroupsVal,
		riskGroupTabVar = riskGroupTabVarVal,
		type = "FFD_Info_List_Save_File"
		)
		
	## Save data:
	#############
	fileName <- tclvalue(tkgetSaveFile(filetypes = 
	    "{{FFD-R Data} {.ffd}} {{All files} *}", defaultextension = ".ffd"))
    if (nchar(fileName) > 0){
        save(FFD_InfoList, file = fileName)
    }		
}

loadFFD <- function(inputDataVar, herdSizeColVar, riskGroupColVar, pi, alpha, 
	piIH, Se, costHerd, costAnimal, sampSizeVar, diagVar, sampStratBox, tn,
	nRiskGroups, riskGroupTabVar){
    fileName <- tclvalue(tkgetOpenFile(filetypes = 
		"{{FFD-R Data} {.ffd}} {{All files} *}")) 	
    if (nchar(fileName) > 0){
        
		## Load file:
		##############
		FFD_InfoList <- NULL
		load(fileName)
		
		## Validity check:
		##################
#		if (!("FFD_InfoList" %in% ls())){
#			tkmessageBox(message = "Error! Invalid File.", icon = "error",
#				type = "ok")
#		    return(NULL)
#		}
		if (is.null(FFD_InfoList)){
			tkmessageBox(message = "Error! Invalid File.", icon = "error",
				type = "ok")
		    return(NULL)
		}
		namesVec <- c("inputDataVar", "herdSizeColVar", "riskGroupColVar", 
			"pi", "alpha", "piIH", "Se", "costHerd", "costAnimal", 
			"sampSizeVar", "diagVar", "sampStratBox", "nRiskGroups",
			"riskGroupTabVar", "type")
	    if (!all(namesVec %in% names(FFD_InfoList))){
			tkmessageBox(message = "Error! Invalid File.", icon = "error",
				type = "ok")
		    return(NULL)			
		} 
		if (FFD_InfoList$type != "FFD_Info_List_Save_File"){
			tkmessageBox(message = "Error! Invalid File.", icon = "error",
				type = "ok")
		    return(NULL)			
		} 
		
		## Set values:
		##############
		tclvalue(inputDataVar) <- FFD_InfoList$inputDataVar
		tclvalue(herdSizeColVar) <- FFD_InfoList$herdSizeColVar
		tclvalue(riskGroupColVar) <- FFD_InfoList$riskGroupColVar
		tclvalue(pi) <- FFD_InfoList$pi
		tclvalue(alpha) <- FFD_InfoList$alpha
		tclvalue(piIH) <- FFD_InfoList$piIH
		tclvalue(Se) <- FFD_InfoList$Se
		tclvalue(costHerd) <- FFD_InfoList$costHerd
		tclvalue(costAnimal) <- FFD_InfoList$costAnimal
		if (FFD_InfoList$sampStratBox != "-1"){
			tcl(sampStratBox, "setvalue", paste("@", FFD_InfoList$sampStratBox,
				sep = ""))			
		} else {
			tcl(sampStratBox, "clearvalue")
		}	
		tcl(tn, "select", 1)
		tkfocus(sampStratBox)
		tclvalue(sampSizeVar) <- FFD_InfoList$sampSizeVar
		tclvalue(diagVar) <- FFD_InfoList$diagVar
		
		## Risk groups:
		nRiskGroupsOldVal <- as.numeric(tclvalue(nRiskGroups))
		nRiskGroupsNewVal <- as.numeric(FFD_InfoList$nRiskGroups)
		tclvalue(nRiskGroups) <- FFD_InfoList$nRiskGroups
		riskGroupTabVarVal <- FFD_InfoList$riskGroupTabVar
		## Delete old values
		if (nRiskGroupsOldVal > nRiskGroupsNewVal){
			for (ii in (nRiskGroupsNewVal+1):nRiskGroupsOldVal){
		        for (jj in 1:4){
				    riskGroupTabVar[[ii,jj-1]] <- character()
			    }			
		    }			
		}
		## Set new values:
		if (nRiskGroupsNewVal > 0){	
		    for (ii in 1:(nRiskGroupsNewVal+1)){
		        for (jj in 1:4){
			        if (!is.na(riskGroupTabVarVal[ii,jj])){
					    riskGroupTabVar[[ii-1,jj-1]] <- riskGroupTabVarVal[ii,jj] 
				    } else {
						riskGroupTabVar[[ii-1,jj-1]] <- character()
					}				
		        }
	        }		
	    }
 		
    }
}

resetFFD <- function(inputDataVar, herdSizeColVar, riskGroupColVar, pi, 
	alpha, piIH, Se, costHerd, costAnimal, sampSizeVar, diagVar, sampStratBox, 
	labelText, labelText1b, labelText2b, labelText2, entry.ssV, entry.diagV, 
	tn, riskGroupTabVar, nRiskGroups, FFD_DefaultList){
    tclvalue(inputDataVar) <- FFD_DefaultList$inputDataVar
	tclvalue(herdSizeColVar) <- FFD_DefaultList$herdSizeColVar
	tclvalue(riskGroupColVar) <- FFD_DefaultList$riskGroupColVar
	tclvalue(pi) <- FFD_DefaultList$pi
	tclvalue(alpha) <- FFD_DefaultList$alpha
	tclvalue(piIH) <- FFD_DefaultList$piIH
	tclvalue(Se) <- FFD_DefaultList$Se
	tclvalue(costHerd) <- FFD_DefaultList$costHerd
	tclvalue(costAnimal) <- FFD_DefaultList$costAnimal
	tclvalue(sampSizeVar) <- FFD_DefaultList$sampSizeVar
	tclvalue(diagVar) <- FFD_DefaultList$diagVar
	tcl(sampStratBox, "clearvalue")
	
	tkconfigure(entry.ssV, state = "disabled")
	tkconfigure(entry.ssV, width = FFD_DefaultList$resetWidth)
	tkconfigure(entry.diagV, state = "disabled")
	tkconfigure(entry.diagV, width = FFD_DefaultList$resetWidth)
	tclvalue(labelText) <- FFD_DefaultList$textNoStrat
	tclvalue(labelText2) <- FFD_DefaultList$textNoStrat
	tclvalue(labelText1b) <- ""
	tclvalue(labelText2b) <- ""
    tcl(tn, "select", 0)
    tcl(tn, "tab", 2, state = "disabled")
	
	nRiskGroupsVal <- as.numeric(tclvalue(nRiskGroups))
	if (nRiskGroupsVal > 0){
		for (ii in 1:(nRiskGroupsVal+1)){
		    for (jj in 1:4){
			    riskGroupTabVar[[ii-1,jj-1]] <- character()
		    }
	    }		
	}
	tclvalue(nRiskGroups) <- "0"
}

## Example: Brucella melitensis without risk groups:
loadExample <- function(inputDataVar, herdSizeColVar, riskGroupColVar, pi, 
	alpha, piIH, Se, costHerd, costAnimal, sampSizeVar, diagVar, 
	sampStratBox, tn, nRiskGroups, riskGroupTabVar){
    tclvalue(inputDataVar) <- file.path(system.file(package = "FFD"), 
		"sheepData.csv", fsep = .Platform$file.sep)
	tclvalue(herdSizeColVar) <- "nSheep"
    tclvalue(riskGroupColVar) <- ""
	tclvalue(pi) <- "0.2"
	tclvalue(alpha) <- "0.05"
	tclvalue(piIH) <- "12"
	tclvalue(Se) <- "90"
	tclvalue(costHerd) <- "30"
	tclvalue(costAnimal) <- "10"
	tcl(sampStratBox, "setvalue", "@1")
	tcl(tn, "select", 1)
	tkfocus(sampStratBox)
	tclvalue(sampSizeVar) <- "90"
	tclvalue(diagVar) <- "5"

	## Risk groups:
	nRiskGroupsOldVal <- as.numeric(tclvalue(nRiskGroups))
	if (nRiskGroupsOldVal > 0){ 
		for (ii in 1:(nRiskGroupsOldVal+1)){
			for (jj in 1:4){
				riskGroupTabVar[[ii-1,jj-1]] <- character()
	    	}				    
		}			
    }
	tclvalue(nRiskGroups) <- "0"
}

## Example: Brucella melitensis with risk groups:
loadExample2 <- function(inputDataVar, herdSizeColVar, riskGroupColVar, pi, 
	alpha, piIH, Se, costHerd, costAnimal, sampSizeVar, diagVar, 
	sampStratBox, tn, nRiskGroups, riskGroupTabVar){
    tclvalue(inputDataVar) <- file.path(system.file(package = "FFD"), 
		"sheepData.csv", fsep = .Platform$file.sep)
	tclvalue(herdSizeColVar) <- "nSheep"
    tclvalue(riskGroupColVar) <- "import"
	tclvalue(pi) <- "0.2"
	tclvalue(alpha) <- "0.05"
	tclvalue(piIH) <- "12"
	tclvalue(Se) <- "90"
	tclvalue(costHerd) <- "30"
	tclvalue(costAnimal) <- "10"
	tcl(sampStratBox, "setvalue", "@1")
	tcl(tn, "select", 1)
	tkfocus(sampStratBox)
	tclvalue(sampSizeVar) <- "90"
	tclvalue(diagVar) <- "5"

	## Risk groups
	nRiskGroupsOldVal <- as.numeric(tclvalue(nRiskGroups))
	nRiskGroupsNewVal <- 2
	tclvalue(nRiskGroups) <- as.character(nRiskGroupsNewVal)
	riskGroupTabVarVal <- matrix(c("Risk groups", "import", "no_import", 
		"Risk_values", "2", "1", 
		"Fixed_sample_sizes", NA, NA,
		"Sample_probability", "4", "1"), 3, 4)
	## Delete old values
	if (nRiskGroupsOldVal > nRiskGroupsNewVal){
		for (ii in (nRiskGroupsNewVal+1):nRiskGroupsOldVal){
	        for (jj in 1:4){
			    riskGroupTabVar[[ii,jj-1]] <- character()
		    }			
	    }			
	}
	## Set new values:
	for (ii in 1:(nRiskGroupsNewVal+1)){
	    for (jj in 1:4){
	        if (!is.na(riskGroupTabVarVal[ii,jj])){
			    riskGroupTabVar[[ii-1,jj-1]] <- riskGroupTabVarVal[ii,jj] 
		    } else {
				riskGroupTabVar[[ii-1,jj-1]] <- character()
			}				
		}
	}	
}

