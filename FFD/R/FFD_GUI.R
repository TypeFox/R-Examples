## Main GUI-Widget:
###############################################################################
###############################################################################
FFD_GUI <- function(){
	tclRequire("BWidget")  ## For dropdown menu
	tclRequire("Tktable") ## For Tables
	main <- tktoplevel(borderwidth = 10)
	#cat(str(main))
	mainEnvir <- get("parent", envir = main$env)$env
	tkwm.resizable(main, FALSE, FALSE)   # Window not resizable
	tkwm.title(main, "FFD - GUI")
	fontBold <- tkfont.create(size = 8, weight = "bold")	
	
	## Default values:
	FFD_DefaultList <- list(inputDataVar = "",
		herdSizeColVar = "",
		riskGroupColVar = "",
		pi = "",
		alpha = "",
		piIH = "",
		Se = "",
		costHerd = "",
		costAnimal = "",
		sampSizeVar = "",
		diagVar = "",
		textNoStrat = "        No sampling strategy specified        ",
		resetWidth = "0")

    ## Global variables:
#	riskGroupTabVar <<- tclArray()
#	isPopWindowOpen <<- FALSE
	
	## Menu:
	##########################################################
	##########################################################	
	FFDmenu <- tkmenu(main)
    tkconfigure(main, menu = FFDmenu)
    fileMenu <- tkmenu(FFDmenu, tearoff = FALSE)
    helpMenu <- tkmenu(FFDmenu, tearoff = FALSE)
	examplesMenu <- tkmenu(fileMenu, tearoff=FALSE)
	## File Menu:
	tkadd(fileMenu, "command", label = "Load", 
		command = function() loadFFD(inputDataVar, herdSizeColVar, 
		riskGroupColVar, pi, alpha, piIH, Se, costHerd, costAnimal, 
		sampSizeVar, diagVar, sampStratBox, tn, nRiskGroups, riskGroupTabVar))
	tkadd(fileMenu, "command", label = "Save", 
		command = function() saveFFD(inputDataVar, herdSizeColVar, 
		riskGroupColVar, pi, alpha, piIH, Se, costHerd, costAnimal, 
		sampSizeVar, diagVar, sampStratBox, nRiskGroups, riskGroupTabVar))
	tkadd(fileMenu, "separator")
	tkadd(fileMenu, "cascade", label = "Load examples", menu=examplesMenu)
    tkadd(examplesMenu, "command", label = "Brucella melitensis without risk groups", 
		command = function() loadExample(inputDataVar, herdSizeColVar, 
	    riskGroupColVar, pi, alpha, piIH, Se, costHerd, costAnimal, 
		sampSizeVar, diagVar, sampStratBox, tn, nRiskGroups, riskGroupTabVar))
    tkadd(examplesMenu, "command", label = "Brucella melitensis with risk groups", 
		command = function() loadExample2(inputDataVar, herdSizeColVar, 
	    riskGroupColVar, pi, alpha, piIH, Se, costHerd, costAnimal, 
		sampSizeVar, diagVar, sampStratBox, tn, nRiskGroups, riskGroupTabVar))    
    tkadd(fileMenu,"separator")
	tkadd(fileMenu, "command", label = "Reset", 
		command = function() resetFFD(inputDataVar, herdSizeColVar, 
		riskGroupColVar, pi, alpha, 
		piIH, Se, costHerd, costAnimal, sampSizeVar, diagVar, sampStratBox,
		labelText, labelText1b, labelText2b, labelText2, entry.ssV, 
		entry.diagV, tn, riskGroupTabVar, nRiskGroups, FFD_DefaultList))
	tkadd(FFDmenu, "cascade", label = "File", menu = fileMenu)
	##Help Menu:	
	tkadd(helpMenu, "command", label = "R Documentation", command = function() 
		print(help("FFD_GUI", package = "FFD", help_type = "text")))
    tkadd(helpMenu, "command", label = "FFD Manual (PDF)", 
		command = function() FFDmanual())	
	packageInfo <- packageDescription("FFD")
	infoMessage <- paste("FFD v", packageInfo$Version, " (", 
		packageInfo$Date, ")\n", packageInfo$URL, 
		"\nDeveloped and maintained by:\n",
		packageInfo$Maintainer, sep = "")
	tkadd(helpMenu,"separator")
	tkadd(helpMenu, "command", label = "About FFD_GUI", command = function(){ 
		tkmessageBox(message = infoMessage, title = "About FFD_GUI")})	
	tkadd(FFDmenu, "cascade", label = "Help", menu = helpMenu)
	
	
	## Main Frame:
	##########################################################
	frameOverall <- tkframe(main)
	
	## Header Frame:
	##########################################################
	##########################################################
	frameHeader <- tkframe(frameOverall, padx = 10, pady = 10,
		relief = "groove", borderwidth = 0)


    #### Logografik:
    image1 <- tclVar()
	imageFile <- file.path(system.file(package = "FFD"), "GUI_Logo4.ppm", 
		fsep = .Platform$file.sep)
    tcl("image", "create", "photo", image1, file = imageFile)
    #imgAsLabel <- tklabel(frameHeader, image = image1, bg = "white")
#	imgAsLabel <- tklabel(frameHeader, image = image1, bg = "grey82")
	imgAsLabel <- tklabel(frameHeader, image = image1, bg = "#D4D0C8")
	tkgrid(imgAsLabel)	
	
	## Tabs Frame:
	#########################################################
	##########################################################
	tn <- ttknotebook(frameOverall)
	
	## Tab Data Input:
	##########################################################
	tabDataInput <- ttkframe(tn, borderwidth = 17)	
	frameDataInput <- tkframe(tabDataInput, padx = 20,
		relief = "groove", borderwidth = 0)
    frameDataInputInner <- tkframe(frameDataInput, padx = 10, pady = 10,
		relief = "groove", borderwidth = 1)
	tkgrid(tklabel(frameDataInputInner, text = "Farm data: ", font=fontBold), sticky = "w")
	
	## Farm Data (csv):
    inputDataVar <- tclVar(FFD_DefaultList$inputDataVar)
	herdSizeColVar <- tclVar(FFD_DefaultList$herdSizeColVar)
	riskGroupColVar <- tclVar(FFD_DefaultList$riskGroupColVar)
	LoadData.but <- tkbutton(frameDataInputInner, text = "...", 
		command = function() importDataFun(inputDataVar, herdSizeColVar,
		riskGroupColVar)) 
    entry.inputData <- tkentry(frameDataInputInner, width = "21", 
		textvariable = inputDataVar)
	tkgrid(tklabel(frameDataInputInner, text = "Data file: "), entry.inputData,
		LoadData.but, sticky = "e")
    	
	## Herd size vector:	
	ChooseHerdSizeCol.but <- tkbutton(frameDataInputInner, text = "...", 
		command = function() chooseColumnName(inputDataVar, herdSizeColVar, "herd sizes")) 	
	entry.herdSizeCol <- tkentry(frameDataInputInner, width = "21", 
		textvariable = herdSizeColVar)
	tkgrid(tklabel(frameDataInputInner, text = "Herd sizes column: "), 
		entry.herdSizeCol, ChooseHerdSizeCol.but, sticky = "e")

    tkgrid(tklabel(frameDataInputInner, text = ""))

	frameDataInputInnerRG <- tkframe(frameDataInput, padx = 10, pady = 10,
		relief = "groove", borderwidth = 1)
#	tkgrid(tklabel(frameDataInputInnerRG, text = "Risk group data (optional): ", 
#		font=fontBold), sticky = "e")
	tkgrid(tklabel(frameDataInputInnerRG, text = "Risk group data (optional): ", 
		font=fontBold), sticky = "w", columnspan = 2)
#	tkgrid(tklabel(frameDataInputInnerRG, text = "Risk group data", 
#		font=fontBold), tklabel(frameDataInputInnerRG, text = "(optional):", 
#		font=fontBold), sticky = "w")
	## Risk Group vector:	
	ChooseRGCol.but <- tkbutton(frameDataInputInnerRG, text = "...", 
		command = function() chooseColumnName(inputDataVar, riskGroupColVar, "risk groups")) 	
	entry.rgCol <- tkentry(frameDataInputInnerRG, width = "21", 
		textvariable = riskGroupColVar)
	tkgrid(tklabel(frameDataInputInnerRG, text = "Risk group column: "), 
		entry.rgCol, ChooseRGCol.but, sticky = "e")
    tkgrid(tklabel(frameDataInputInnerRG, text = ""))
	
	
	## Tab Parameters:
	##########################################################
	tabParameters <- ttkframe(tn, borderwidth = 0)		
	frameParameters <- tkframe(tabParameters, padx = 20,
		relief = "groove", borderwidth = 0)
	
    ## Required parameters:
    frameInput <- tkframe(frameParameters,padx = 10, pady = 10,
		relief = "groove", borderwidth = 1)
    tkgrid(tklabel(frameInput, text = "Required parameters: ", font=fontBold), sticky = "e")
    ## Design prevalence:
    pi <- tclVar(FFD_DefaultList$pi)
    entry.pi <- tkentry(frameInput, width = "10", textvariable = pi)
	tkgrid(tklabel(frameInput,text = "Design prevalence: "), entry.pi, 
		tklabel(frameInput,text ="%"), sticky = "e")
	## Intraherd prevalence:
	piIH <- tclVar(FFD_DefaultList$piIH)
    entry.piIH <- tkentry(frameInput, width = "10", textvariable = piIH)
	tkgrid(tklabel(frameInput,text = "Intraherd prevalence: "), 
		entry.piIH, tklabel(frameInput,text ="%"), sticky = "e")
    ## Test sensitivity:
	Se <- tclVar(FFD_DefaultList$Se)
    entry.Se <- tkentry(frameInput, width="10", textvariable = Se)
	tkgrid(tklabel(frameInput,text = "Test sensitivity: "), 
		entry.Se, tklabel(frameInput,text ="%"), sticky = "e")  
    ## Significance level:
	alpha <- tclVar(FFD_DefaultList$alpha)
    entry.alpha <- tkentry(frameInput, width="10", textvariable = alpha)
	tkgrid(tklabel(frameInput,text = "Type I error level (alpha): "), entry.alpha, sticky = "e")
	
    ## Sampling strategy:
    sampStratVec <- c("limited sampling", "individual sampling")
	sampStratBox <- tkwidget(frameInput, "ComboBox", editable = FALSE,
		values = sampStratVec, width = "19")
    tkgrid(tklabel(frameInput,text = "Sampling strategy: "), 
		sampStratBox, sticky = "e")
    tkbind(sampStratBox,"<FocusIn>", 
		function() testFocusFun())
    tkbind(sampStratBox,"<FocusOut>", 
		function() changeTextSampSize(sampStratBox, sampStratVec))


    ## Optional parameters:
    frameOptionalInput <- tkframe(frameParameters,padx = 10, pady = 10,
		relief = "groove", borderwidth = 1)
    tkgrid(tklabel(frameOptionalInput, text = "Optional parameters: ", 
		font = fontBold), sticky = "e")
    ## Cost per herd:
    costHerd <- tclVar(FFD_DefaultList$costHerd)
    entry.costHerd <- tkentry(frameOptionalInput, width="10", 
		textvariable = costHerd)
	tkgrid(tklabel(frameOptionalInput, text = "Cost per herd: "), entry.costHerd, 
		sticky = "e")
	## Cost per animal:
    costAnimal <- tclVar(FFD_DefaultList$costAnimal)
    entry.costAnimal <- tkentry(frameOptionalInput, width="10", 
		textvariable = costAnimal)
	tkgrid(tklabel(frameOptionalInput, text = "Cost per animal: "), entry.costAnimal, 
		sticky = "e")	
	## Specify risk group parameters:
	riskGroupTabVar <- tclArray()
	nRiskGroups <- tclVar("0")
	## Set global variable for the id of the "set risk groups"-window:
    assign("setRGWindow_ID", "windowNotYetInitialized", envir = mainEnvir)
#    assign("setRGWindow_ID", "windowNotYetInitialized", envir = .GlobalEnv)
    setRGParameters.but <- tkbutton(frameOptionalInput, text = " Set risk group parameters ", 
		command = function() setRGParameters(inputDataVar, riskGroupColVar, 
		riskGroupTabVar, nRiskGroups, mainEnvir), 
		padx = 15, pady = 2)    
    tkgrid(setRGParameters.but, columnspan = 2)


	## Bundle Information:
	infoList <- list(pi = pi, alpha = alpha, piIH = piIH, Se = Se, 
		costHerd = costHerd, costAnimal = costAnimal, 
		sampStratVec = sampStratVec, sampStratBox = sampStratBox,
		inputDataVar = inputDataVar, herdSizeColVar = herdSizeColVar,
		riskGroupColVar = riskGroupColVar, riskGroupTabVar = riskGroupTabVar,
		nRiskGroups = nRiskGroups)

	
	## Tab Sample Size:
	##########################################################
	tabCalculations <- ttkframe(tn, borderwidth = 10)
	frameSampleSize <- tkframe(tabCalculations, padx = 40,
		relief = "groove", borderwidth = 0)    
    
    ## Computation of Sample Size:
    frameSampleSizeInner <- tkframe(frameSampleSize,padx = 10, pady = 10,
		relief = "groove", borderwidth = 1)    
    frameSampleSizeInnerTop <- tkframe(frameSampleSizeInner, 
		relief = "groove", borderwidth = 0)    
	labelText <- tclVar(FFD_DefaultList$textNoStrat)
    label1 <- tklabel(frameSampleSizeInnerTop, text = tclvalue(labelText))
    tkconfigure(label1, textvariable = labelText)
	labelText1b <- tclVar("")
    label1b <- tklabel(frameSampleSizeInnerTop, text = tclvalue(labelText1b), width="1")
    tkconfigure(label1b, textvariable = labelText1b)
	sampSizeVar <- tclVar(FFD_DefaultList$sampSizeVar)
    entry.ssV <- tkentry(frameSampleSizeInnerTop, width="10", 
		textvariable = sampSizeVar)
    tkconfigure(entry.ssV, state = "disabled")
	tkconfigure(entry.ssV, width = FFD_DefaultList$resetWidth)
	
	tkgrid(tklabel(frameSampleSizeInnerTop, 
		text = "Compute sample size: ", font = fontBold), sticky = "w")
    tkgrid(label1, entry.ssV, label1b, sticky = "e")
	tkgrid(tklabel(frameSampleSizeInnerTop, text = ""))
	tkgrid(frameSampleSizeInnerTop)
	
	frameSampleSizeInnerBottom <- tkframe(frameSampleSizeInner, 
		relief = "groove", borderwidth = 0) 
	CalcSS.but <- tkbutton(frameSampleSizeInnerBottom, text = " Calculate ", 
		command = function() OnCalcSS(infoList,sampSizeVar,main), padx = 15, pady = 2)  
    SampleSS.but <- tkbutton(frameSampleSizeInnerBottom, text = " Sample ", 
		command = function() OnSampleSS(infoList,sampSizeVar,main), padx = 15, pady = 2) 
    #tkgrid(CalcSS.but, columnspan = 2)
	tkgrid(CalcSS.but, tklabel(frameSampleSizeInnerBottom, text = " "), SampleSS.but)
    tkgrid(frameSampleSizeInnerBottom)
	
	## Diagnostics:
	frameSampleSizeDiag <- tkframe(frameSampleSize,padx = 10, pady = 10,
		relief = "groove", borderwidth = 1)    
	labelText2 <- tclVar(FFD_DefaultList$textNoStrat)
    label2 <- tklabel(frameSampleSizeDiag, text = tclvalue(labelText2))
    tkconfigure(label2, textvariable = labelText2)
	labelText2b <- tclVar("")
    label2b <- tklabel(frameSampleSizeDiag, text = tclvalue(labelText2b), width="1")
    tkconfigure(label2b, textvariable = labelText2b)
	diagVar <- tclVar(FFD_DefaultList$diagVar)
    entry.diagV <- tkentry(frameSampleSizeDiag, width="10", 
		textvariable = diagVar)
    tkconfigure(entry.diagV, state = "disabled")
	tkconfigure(entry.diagV, width = FFD_DefaultList$resetWidth)
	tkgrid(tklabel(frameSampleSizeDiag, 
		text = "Sample size diagnostics: ", font = fontBold), sticky = "w")
    tkgrid(label2, entry.diagV, label2b, sticky = "e")
	tkgrid(tklabel(frameSampleSizeDiag, text = ""))
	
	CalcDiag.but <- tkbutton(frameSampleSizeDiag, text = " Calculate ", 
		command = function() OnCalcDiag(infoList,diagVar,main), 
		padx = 15, pady = 2)    
    tkgrid(CalcDiag.but, columnspan = 2)
	

	## Build Tabs:
	##########################################################
	##########################################################
	tkadd(tn, tabDataInput, text = " Data Input ") ### tabid=0
    tkadd(tn, tabParameters, text = " Parameters ") ### tabid=1  
    tkadd(tn, tabCalculations, text = " Calculations ")  ### tabid=2    
		
	## Buttons Frame:
	##########################################################
	##########################################################
	frameButton <- tkframe(frameOverall, relief = "groove", padx = 10, pady = 10,
		borderwidth = 0)	
    Close.but <- tkbutton(frameButton, text = " Close ", 
		command = function() tkdestroy(main), padx = 15, pady = 2)
    tkgrid(Close.but)	

    ## Internal Function:
	##########################################################
	##########################################################
	
	testFocusFun <- function() tkfocus(entry.costHerd)
	
    changeTextSampSize <- function(textBox, textVec){	
    	boxValue <- textVec[as.numeric(tclvalue(tcl(textBox,"getvalue")))+1]
        if (length(boxValue) > 0){
		    if (boxValue == "limited sampling"){
    		    if (tclvalue(labelText) == "Herd sensitivity: "){				
				    tclvalue(sampSizeVar) <- FFD_DefaultList$sampSizeVar
	                tclvalue(diagVar) <- FFD_DefaultList$diagVar
			    }
    		    if (tclvalue(labelText) != "Sample limit: "){				
				    tclvalue(labelText) <- "Sample limit: "
                    tclvalue(labelText1b) <- ""
                    tclvalue(labelText2b) <- ""
                    tclvalue(labelText2) <- "Max. sample limit: "   
					tkconfigure(entry.ssV, state = "normal")
	                tkconfigure(entry.ssV, width = "10")
					tkconfigure(entry.diagV, state = "normal")
	                tkconfigure(entry.diagV, width = "10")
					tcl(tn, "tab", 2, state="normal")
			    }
	        } else {
				if (tclvalue(labelText) == "Sample limit: "){
    		        tclvalue(sampSizeVar) <- FFD_DefaultList$sampSizeVar
	                tclvalue(diagVar) <- FFD_DefaultList$diagVar
			    }
				if (tclvalue(labelText) != "Herd sensitivity: "){
    		        tclvalue(labelText) <- "Herd sensitivity: "
                    tclvalue(labelText1b) <- "%"
                    tclvalue(labelText2b) <- "%"
				    tclvalue(labelText2) <- "Herd sens. step size: "  
					tkconfigure(entry.ssV, state = "normal")
	                tkconfigure(entry.ssV, width = "10")
					tkconfigure(entry.diagV, state = "normal")
	                tkconfigure(entry.diagV, width = "10")
					tcl(tn, "tab", 2, state="normal")
			    }
	        }				
			#tkfocus(entry.costHerd)
	    }
    }	
	
	## Build Main Frame:
	##########################################################
	##########################################################	
	tkgrid(frameHeader, pady=10)
	tkgrid(frameDataInputInner, pady=10)
	tkgrid(frameDataInputInnerRG, pady=10)
	tkgrid(frameDataInput)	
	tkgrid(frameInput, pady=10)
	tkgrid(frameOptionalInput, pady = 10)
	tkgrid(frameParameters,  sticky = "e")	
	tkgrid(frameSampleSizeInner, pady=10)
	tkgrid(frameSampleSizeDiag, pady=10)
	tkgrid(frameSampleSize)	
	tkgrid(tn)
	tkgrid(frameButton)
	tkgrid(frameOverall)
	tcl(tn, "tab", 2, state="disabled")
	tkfocus(main)	
}
