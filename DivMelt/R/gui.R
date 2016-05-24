require(tcltk)

DivMelt_gui <- function()
{
    tt <- tktoplevel()
    tkwm.title(tt,"HRM Diversity Assay Analysis Tool")
    tkwm.geometry(tt,"400x350+0+0")

    launchTime <- Sys.time()
    fontForHelp <- tkfont.create(family="times",size=12,weight="bold",slant="italic")

    # Set all variable setting defaults
    epcValue <- tclVar("1")
    earlyValue <- tclVar("3")
    lpcValue <- tclVar("1")
    lateValue <- tclVar("3")
    hhtValue <- tclVar("1")
    shoulderValue <- tclVar("10")
    t1Value <- tclVar("50")
    t1HoldValue <- tclVar("1")
    t1SlopeWinValue <- tclVar("1")
    t2Value <- tclVar("30")
    t2HoldValue <- tclVar("1")
    t2SlopeWinValue <- tclVar("1")
    dirName <- tclVar(getwd())
    fileName <- tclVar("")
    sampleName <- tclVar("")
    plotTypeValue <- tclVar("PDF")
    fluorSmoothCBValue <- tclVar("0")
    dfluorSmoothCBValue <- tclVar("1")
    outdirName <- tclVar(getwd())
    outfileName <- tclVar("")
    statsName <- tclVar("")
    scoreName <- tclVar("")
    lassoValue <- tclVar("1")
    modelName <- tclVar("")
    xrangeValue <- tclVar("0")
    xminValue <- tclVar("70.")
    xmaxValue <- tclVar("95.")
    pdfVName <- tclVar("acroread")
    pngVName <- tclVar("display")

    loadSettings <- function(settingsFile)
    {
	if (!nchar(settingsFile) || settingsFile == "")
	{
	    return
	}
	data <- matrix (scan (file=settingsFile, skip=0, what='character',
		    sep="="), ncol=2, byrow=TRUE)
	for (i in (1:nrow(data)))
	{
	    if (data[i,1]=="sample")
		tclvalue(sampleName) <- data[i,2]
	    else if (data[i,1]=="plottype")
		tclvalue(plotTypeValue) <- data[i,2]
	    else if (data[i,1]=="pdfViewer")
		tclvalue(pdfVName) <- data[i,2]
	    else if (data[i,1]=="pngViewer")
		tclvalue(pngVName) <- data[i,2]
	    else if (data[i,1]=="smoothFluor")
		tclvalue(fluorSmoothCBValue) <- data[i,2]
	    else if (data[i,1]=="smoothdFluor")
		tclvalue(dfluorSmoothCBValue) <- data[i,2]
	    else if (data[i,1]=="excludeEarly")
		tclvalue(epcValue) <- data[i,2]
	    else if (data[i,1]=="earlyCutoff")
		tclvalue(earlyValue) <- data[i,2]
	    else if (data[i,1]=="excludeLate")
		tclvalue(lpcValue) <- data[i,2]
	    else if (data[i,1]=="lateCutoff")
		tclvalue(lateValue) <- data[i,2]
	    else if (data[i,1]=="includeShoulders")
		tclvalue(hhtValue) <- data[i,2]
	    else if (data[i,1]=="shoulderHeight")
		tclvalue(shoulderValue) <- data[i,2]
	    else if (data[i,1]=="useLasso")
		tclvalue(lassoValue) <- data[i,2]
	    else if (data[i,1]=="lassoModel")
		tclvalue(modelName) <- data[i,2]
	    else if (data[i,1]=="theta1")
		tclvalue(t1Value) <- data[i,2]
	    else if (data[i,1]=="theta1Hold")
		tclvalue(t1HoldValue) <- data[i,2]
	    else if (data[i,1]=="theta1Window")
		tclvalue(t1SlopeWinValue) <- data[i,2]
	    else if (data[i,1]=="theta2")
		tclvalue(t2Value) <- data[i,2]
	    else if (data[i,1]=="theta2Hold")
		tclvalue(t2HoldValue) <- data[i,2]
	    else if (data[i,1]=="theta2Window")
		tclvalue(t2SlopeWinValue) <- data[i,2]
	    else if (data[i,1]=="useXrange")
		tclvalue(xrangeValue) <- data[i,2]
	    else if (data[i,1]=="minX")
		tclvalue(xminValue) <- data[i,2]
	    else if (data[i,1]=="maxX")
		tclvalue(xmaxValue) <- data[i,2]
	    else 
	        tkmessageBox(message=paste("Unrecognized keyword",data[i,1],"in file",settingsFile,"ignored.",sep=" "),icon="warning")
		
	}
    }

    if (file.exists("settings.txt"))
    {
	loadSettings("settings.txt")
    }

    runCommand <- function()
    {
	showOptions(1)
    }
    showCommand <- function()
    {
	showOptions(0)
    }
    showOptions <- function(doRun)
    {
	dlg <- tktoplevel()
	tkwm.deiconify(dlg)
	tkgrab.set(dlg)
	tkfocus(dlg)
	title = "HRM Diversity Assay Analysis Tool"
	tkwm.title(dlg, title)
	#tkwm.geometry(dlg,"400x200+0+0")
	tkwm.geometry(dlg,"+0+0")

        runFrame <- tkframe(dlg, borderwidth=2, relief="groove")
	tkgrid.columnconfigure(runFrame,1,weight=1)
	tkgrid.columnconfigure(runFrame,2,weight=1)
	tkgrid.configure(runFrame,sticky="nswe")

	outdirN <- as.character(tclvalue(outdirName))
	dirN <- as.character(tclvalue(dirName))
	fileN <- as.character(tclvalue(fileName))
	sampleN <- as.character(tclvalue(sampleName))
	plotTypeVal <- tolower(as.character(tclvalue(plotTypeValue)))
	fsmCB <- as.character(tclvalue(fluorSmoothCBValue))
	dfsmCB <- as.character(tclvalue(dfluorSmoothCBValue))
	epcV <- as.character(tclvalue(epcValue))
	lpcV <- as.character(tclvalue(lpcValue))
	hhtV <- as.character(tclvalue(hhtValue))
	lassoV <- as.character(tclvalue(lassoValue))
	earlyV <- as.character(tclvalue(earlyValue))
	lateV <- as.character(tclvalue(lateValue))
	shoulderV <- as.character(tclvalue(shoulderValue))
	modelF <- as.character(tclvalue(modelName))
	outF <- as.character(tclvalue(outfileName))
	scoreF <- as.character(tclvalue(scoreName))
	statsF <- as.character(tclvalue(statsName))
	t1V <- as.character(tclvalue(t1Value))
	t1HoldV <- as.character(tclvalue(t1HoldValue))
	t1SlopeW <- as.character(tclvalue(t1SlopeWinValue))
	t2V <- as.character(tclvalue(t2Value))
	t2HoldV <- as.character(tclvalue(t2HoldValue))
	t2SlopeW <- as.character(tclvalue(t2SlopeWinValue))
	xrange <- as.character(tclvalue(xrangeValue))
	xmin <- as.character(tclvalue(xminValue))
	xmax <- as.character(tclvalue(xmaxValue))

	if (doRun)
	    resp= tklabel(runFrame,text="Perform HRM Diversity Assay Analysis with the following options?",font=fontForHelp)
	else
	    resp= tklabel(runFrame,text="You have selected the following options:",font=fontForHelp)

	tkgrid(resp)

	opts<-paste("Input directory: ",dirN)
	info <- tklabel(runFrame,text=opts)
	tkgrid(info)
	if (fileN != "")
	{
	    opts<-paste("Filename: ",fileN)
	    info <- tklabel(runFrame,text=opts)
	    tkgrid(info)
	}

	if (sampleN != "")
	{
	    opts<-paste("Sample name/mask: ",sampleN)
	    info <- tklabel(runFrame,text=opts)
	    tkgrid(info)
	}

	opts<-paste("Plot type: ",plotTypeVal)

	if (fsmCB != 0)
	    opts<-paste("Show Fluorescence smoothing: on")
	else
	    opts<-paste("Show Fluorescence smoothing: off")
	info <- tklabel(runFrame,text=opts)
	tkgrid(info)

	if (dfsmCB != 0)
	    opts<-paste("Show dFluorescence smoothing: on")
	else
	    opts<-paste("Show dFluorescence smoothing: off")
	info <- tklabel(runFrame,text=opts)
	tkgrid(info)

	if (epcV != 0)
	{
	    opts<-paste("Early peak cutoff: ",earlyV," degrees C")
	    info <- tklabel(runFrame,text=opts)
	    tkgrid(info)
	}

	if (lpcV != 0)
	{
	    opts<-paste("Late peak cutoff: ",lateV," degrees C")
	    info <- tklabel(runFrame,text=opts)
	    tkgrid(info)
	}
	
	if (hhtV != 0)
	{
	    opts<-paste("Shoulder height threshold: ",shoulderV," percent")
	    info <- tklabel(runFrame,text=opts)
	    tkgrid(info)
	}
	
	if (lassoV != 0)
	{
	    if (modelF != "")
		opts<-paste("Lasso model: ",modelF)
	    else
		opts<-paste("Default Lasso model")
	    info <- tklabel(runFrame,text=opts)
	    tkgrid(info)
	}
	opts<-paste("Theta 1: ",t1V," degrees")
	info <- tklabel(runFrame,text=opts)
	tkgrid(info)
	opts<-paste("T1 hold interval: ",t1HoldV," degrees C")
	info <- tklabel(runFrame,text=opts)
	tkgrid(info)
	opts<-paste("T1 slope window: ",t1SlopeW," degrees C")
	info <- tklabel(runFrame,text=opts)
	tkgrid(info)
	opts<-paste("Theta 2: ",t2V," degrees")
	info <- tklabel(runFrame,text=opts)
	tkgrid(info)
	opts<-paste("T2 hold interval: ",t2HoldV," degrees C")
	info <- tklabel(runFrame,text=opts)
	tkgrid(info)
	opts<-paste("T2 slope window: ",t2SlopeW," degrees C")
	info <- tklabel(runFrame,text=opts)
	tkgrid(info)
	if (xrange != 0)
	{
	   opts<-paste("Plot X min: ",xmin)
	   info <- tklabel(runFrame,text=opts)
	   tkgrid(info)
	   opts<-paste("Plot X max: ",xmax)
	   info <- tklabel(runFrame,text=opts)
	   tkgrid(info)
	}
	
	opts<-paste("Output directory: ",outdirN)
	info <- tklabel(runFrame,text=opts)
	tkgrid(info)

	if (outF != "")
	{
	    opts<-paste("Output file: ",outF)
	    info <- tklabel(runFrame,text=opts)
	    tkgrid(info)
	}

	if (statsF != "")
	{
	    opts<-paste("Stats file: ",statsF)
	    info <- tklabel(runFrame,text=opts)
	    tkgrid(info)
	}
	
	if (scoreF != "")
	{
	    opts<-paste("HRM score file: ",scoreF)
	    info <- tklabel(runFrame,text=opts)
	    tkgrid(info)
	}
	

	tkgrid(runFrame,sticky="nswe")

	done <- tclVar(0)   # tclVar() creates a Tcl variable

        butFrame <- tkframe(dlg, borderwidth=2, relief="groove")
	tkgrid.configure(butFrame,sticky="nswe")
	tkgrid.columnconfigure(butFrame,0,weight=1)

	OK.but <- tkbutton(butFrame, text = "OK",
	    command = function() tclvalue(done) <- 1)


	if (doRun)
	{
	    Cancel.but <- tkbutton(butFrame, text = "Cancel",
		command = function() tclvalue(done) <- 2)

	    tkgrid(OK.but,Cancel.but)
	    tkgrid.configure(Cancel.but,sticky="nswe")
	    tkgrid.columnconfigure(butFrame,1,weight=1)
	}
	else
	    tkgrid(OK.but)

	tkgrid.configure(OK.but,sticky="nswe")

	tkwait.variable(done)

	if (tclvalue(done) == 1 && doRun) 
	{
	    launchTime <- Sys.time()
	    pb <- tkProgressBar("Analysis progress bar", "progress in %",
				0, 100, 0)

	    updateBar<-function(fractDone)
	    {
		i <- round(100*fractDone)
		    
		info <- sprintf("Analyzing samples (%d%% done)...", i)

		setTkProgressBar(pb, i, title="Progress of analysis...",label=sprintf("%s", info))
		tcl("update", "idletasks")
	    }
	    setTkProgressBar(pb, 0, title="Progress of analysis...",label="Preparing for analysis...")

	    tryCatch({
		processDivMelt(dir=dirN,
		    singleFile=fileN,
		    format=plotTypeVal,
		    excludeEarly=(epcV!=0),earlyCutoff=as.double(earlyV),
		    excludeLate=(lpcV!=0),lateCutoff=as.double(lateV),
		    includeShoulders=(hhtV!=0),shoulderCutoff=as.double(shoulderV)/100.0,
		    showFlorSmooth=(fsmCB != 0),
		    showDFlorSmooth=(dfsmCB != 0),
		    useLasso=(lassoV!=0),
		    modelFile=modelF,
		    dataFile=outF,
		    statsFile=statsF,
		    scoreFile=scoreF,
		    sampleMask=sampleN,
		    theta1=as.double(t1V),
		    theta2=as.double(t2V),
		    t1Cutoff=as.double(t1HoldV),
		    t2Cutoff=as.double(t2HoldV),
		    t1SlopeWindow=as.double(t1SlopeW),
		    t2SlopeWindow=as.double(t2SlopeW),
		    xrangeOverride=(xrange != 0),
		    xrangemin=as.double(xmin),
		    xrangemax=as.double(xmax),
		    progressFunc=updateBar)

		tkmessageBox(message=paste("Processing complete.  Click 'View new files' to see generated files.",sep=" "))

	    }, interrupt = function(ex) {
	        print(paste("Interrupt:",ex,"\n"))
	        tkmessageBox(message="Processing interrupted!",icon="warning")
	    }, error = function(ex) {
	        print(paste("Error:",ex,"\n"))
	        tkmessageBox(message="An error occurred.  Check R console output for details.",icon="warning")
	    }, finally = {
		close(pb)
	    }) # tryCatch()
	}
        tkgrab.release(dlg)
        tkdestroy(dlg)
        tkfocus(tt)
    }

    inputDlg <- function()
    {
	header="Input Options"
	dlg <- tktoplevel()
	tkwm.deiconify(dlg)
	tkgrab.set(dlg)
	tkfocus(dlg)
	title = "HRM Diversity Assay Analysis Tool"
	tkwm.title(dlg, title)
	#tkwm.geometry(dlg,"400x200+0+0")
	tkwm.geometry(dlg,"+0+0")

        frameSources <- tkwidget(dlg, "labelframe",  text=header)
	tkgrid.columnconfigure(frameSources,1,weight=1)
	tkgrid.columnconfigure(frameSources,2,weight=1)
	tkgrid.configure(frameSources,sticky="nswe")

	dirInfo <- tklabel(frameSources,text="i",font=fontForHelp)
	dir.but <- tkbutton(frameSources,text="Input directory",command=selectDir)
	tkbind(dirInfo, "<Button-1>",dirHelp)

	dirValueText <- tkentry(frameSources,textvariable=dirName)
	tkconfigure(dirValueText,width="40",textvariable=dirName)

	fileInfo <- tklabel(frameSources,text="i",font=fontForHelp)
	file.but <- tkbutton(frameSources,text="Input file (opt)",command=selectFile)
	
	tkbind(fileInfo, "<Button-1>",fileHelp)

	fileValueText <- tkentry(frameSources,width="40",textvariable=fileName)

	sampleInfo <- tklabel(frameSources,text="i",font=fontForHelp)
	sampleLab <- tklabel(frameSources,text="Sample name (opt)")
	sampleValueText <- tkentry(frameSources,width="40",textvariable=sampleName)
	tkbind(sampleInfo, "<Button-1>",sampleHelp)

	tkgrid(dirInfo,dir.but,dirValueText)
	tkgrid.configure(dirInfo,sticky="w")
	tkgrid.configure(dir.but,sticky="we")
	tkgrid.configure(dirValueText,sticky="we")

	tkgrid(fileInfo,file.but,fileValueText)
	tkgrid.configure(fileInfo,sticky="w")
	tkgrid.configure(file.but,sticky="we")
	tkgrid.configure(fileValueText,sticky="we")

	tkgrid(sampleInfo,sampleLab,sampleValueText)
	tkgrid.configure(sampleInfo,sticky="w")
	tkgrid.configure(sampleLab,sticky="we")
	tkgrid.configure(sampleValueText,sticky="we")

	tkgrid(frameSources,sticky="nswe")

	done <- tclVar(0)   # tclVar() creates a Tcl variable

	OK.but <- tkbutton(dlg, text = "OK",
	    command = function() tclvalue(done) <- 1)
	tkgrid(OK.but)
	tkgrid.configure(OK.but,sticky="nswe")
	tkwait.variable(done)
        tkgrab.release(dlg)
        tkdestroy(dlg)
        tkfocus(tt)
    }

    plottingDlg <- function()
    {
	header="Plotting Options"
	dlg <- tktoplevel()
	tkwm.deiconify(dlg)
	tkgrab.set(dlg)
	tkfocus(dlg)
	title = "HRM Diversity Assay Analysis Tool"
	tkwm.title(dlg, title)
	#tkwm.geometry(dlg,"400x200+0+0")
	tkwm.geometry(dlg,"+0+0")

        inputFrame <- tkwidget(dlg, "labelframe",  text=header)
	tkgrid.columnconfigure(inputFrame,1,weight=1)
	tkgrid.columnconfigure(inputFrame,2,weight=1)
	tkgrid.columnconfigure(inputFrame,3,weight=1)
	tkgrid.columnconfigure(inputFrame,4,weight=1)
	tkgrid.columnconfigure(inputFrame,5,weight=1)
	tkgrid.columnconfigure(inputFrame,6,weight=1)

	plotInfo <- tklabel(inputFrame,text="i",font=fontForHelp)
	tkbind(plotInfo, "<Button-1>",plotHelp)
	rb1 <- tkradiobutton(inputFrame)
	rb2 <- tkradiobutton(inputFrame)
	rb3 <- tkradiobutton(inputFrame)

	tkconfigure(rb1,variable=plotTypeValue,value="PDF")
	tkconfigure(rb2,variable=plotTypeValue,value="PNG")
	tkconfigure(rb3,variable=plotTypeValue,value="None")

	opt1Lab=tklabel(inputFrame,text="PDF ")
	opt2Lab=tklabel(inputFrame,text="PNG ")
	opt3Lab=tklabel(inputFrame,text="None ")

	tkgrid(plotInfo,opt1Lab,rb1,opt2Lab,rb2,opt3Lab,rb3)
	tkgrid.configure(plotInfo,sticky="w")
	tkgrid.configure(opt1Lab,sticky="w")
	tkgrid.configure(rb1,sticky="w")
	tkgrid.configure(opt2Lab,sticky="w")
	tkgrid.configure(rb2,sticky="w")
	tkgrid.configure(opt3Lab,sticky="w")
	tkgrid.configure(rb3,sticky="w")

	#option for including linear fit line
	fluorSmoothLine <- tklabel(inputFrame,text="i",font=fontForHelp)
	tkbind(fluorSmoothLine, "<Button-1>",fluorSmoothHelp)
	fluorSmoothCB <- tkcheckbutton(inputFrame)
	tkconfigure(fluorSmoothCB,variable=fluorSmoothCBValue)
	fluorSmoothLab=tklabel(inputFrame,text="Show Fluorescence smoothing")

	tkgrid(fluorSmoothLine,fluorSmoothLab,fluorSmoothCB)

	tkgrid.configure(fluorSmoothLine,sticky="w")
	tkgrid.configure(fluorSmoothLab,sticky="w")
	tkgrid.configure(fluorSmoothCB,sticky="w")

	#option for including dfluorescence smoothing line
	dfluorSmoothLine <- tklabel(inputFrame,text="i",font=fontForHelp)
	tkbind(dfluorSmoothLine, "<Button-1>",dfluorSmoothHelp)
	dfluorSmoothCB <- tkcheckbutton(inputFrame)
	tkconfigure(dfluorSmoothCB,variable=dfluorSmoothCBValue)
	dfluorSmoothLab=tklabel(inputFrame,text="Show dFluorescence smoothing")

	tkgrid(dfluorSmoothLine,dfluorSmoothLab,dfluorSmoothCB)
	tkgrid.configure(dfluorSmoothLine,sticky="w")
	tkgrid.configure(dfluorSmoothLab,sticky="w")
	tkgrid.configure(dfluorSmoothCB,sticky="w")

	#option for overriding X-axis scale
	xrangeInfo <- tklabel(inputFrame,text="i",font=fontForHelp)
	tkbind(xrangeInfo, "<Button-1>",xrangeHelp)
	xrange <- tkcheckbutton(inputFrame)
	tkconfigure(xrange,variable=xrangeValue)
	xrangeLab = tklabel(inputFrame,text="Override X-axis range")

	xminLab <- tklabel(inputFrame,text="X min:")
	xminValueText <- tkentry(inputFrame,width="3",textvariable=xminValue)
	tkconfigure(xminValueText,textvariable=xminValue)
    
	xmaxLab <- tklabel(inputFrame,text="X max:")
	xmaxValueText <- tkentry(inputFrame,width="3",textvariable=xmaxValue)
	tkconfigure(xmaxValueText,textvariable=xmaxValue)
    
	tkgrid(xrangeInfo,xrangeLab,xrange,xminLab,xminValueText,xmaxLab,xmaxValueText)
	tkgrid.configure(xrangeInfo,sticky="w")
	tkgrid.configure(xrangeLab,sticky="w")
	tkgrid.configure(xrange,sticky="w")

	tkgrid.configure(xminLab,sticky="w")
	tkgrid.configure(xminValueText,sticky="e")
	tkgrid.configure(xmaxLab,sticky="w")
	tkgrid.configure(xmaxValueText,sticky="e")

	if(.Platform$OS.type == "unix") {
	    pdfVInfo <- tklabel(inputFrame,text="i",font=fontForHelp)
	    tkbind(pdfVInfo, "<Button-1>",pdfVHelp)

	    pdfVLab <- tklabel(inputFrame,text="PDF viewer program")
	    pdfVValueText <- tkentry(inputFrame,width="10",textvariable=pdfVName)
	    tkconfigure(pdfVValueText,textvariable=pdfVName)

	    pngVInfo <- tklabel(inputFrame,text="i",font=fontForHelp)
	    tkbind(pngVInfo, "<Button-1>",pngVHelp)

	    pngVLab <- tklabel(inputFrame,text="PNG viewer program")
	    pngVValueText <- tkentry(inputFrame,width="10",textvariable=pngVName)
	    tkconfigure(pngVValueText,textvariable=pngVName)

	    tkgrid(pdfVInfo,pdfVLab,pdfVValueText)
	    tkgrid.configure(pdfVInfo,sticky="w")
	    tkgrid.configure(pdfVLab,sticky="w")
	    tkgrid.configure(pdfVValueText,sticky="e")
	    tkgrid(pngVInfo,pngVLab,pngVValueText)
	    tkgrid.configure(pngVInfo,sticky="w")
	    tkgrid.configure(pngVLab,sticky="w")
	    tkgrid.configure(pngVValueText,sticky="e")
	}

	tkgrid(inputFrame,sticky="nswe")

	done <- tclVar(0)   # tclVar() creates a Tcl variable

	OK.but <- tkbutton(dlg, text = "OK",
	    command = function() tclvalue(done) <- 1)
	tkgrid(OK.but)
	tkgrid.configure(OK.but,sticky="nswe")
	tkwait.variable(done)
        tkgrab.release(dlg)
        tkdestroy(dlg)
        tkfocus(tt)
    }

    analysisDlg <- function()
    {
	header="Analysis Options"
	dlg <- tktoplevel()
	tkwm.deiconify(dlg)
	tkgrab.set(dlg)
	tkfocus(dlg)
	title = "HRM Diversity Assay Analysis Tool"
	tkwm.title(dlg, title)
	#tkwm.geometry(dlg,"400x400+0+0")
	tkwm.geometry(dlg,"+0+0")

        analysisFrame <- tkwidget(dlg, "labelframe",  text=header)
	tkgrid.columnconfigure(analysisFrame,1,weight=1)
	tkgrid.columnconfigure(analysisFrame,3,weight=1)

	#option for excluding early peaks
	earlyInfo <- tklabel(analysisFrame,text="i",font=fontForHelp)
	tkbind(earlyInfo, "<Button-1>",earlyHelp)
	epc <- tkcheckbutton(analysisFrame)
	tkconfigure(epc,variable=epcValue)
	earlyLab=tklabel(analysisFrame,text="Exclude early peaks")

	tkgrid(earlyInfo,earlyLab,epc)

	tkgrid.configure(earlyInfo,sticky="w")
	tkgrid.configure(earlyLab,sticky="w")
	tkgrid.configure(epc,sticky="w")

	earlyCutoffInfo <- tklabel(analysisFrame,text="i",font=fontForHelp)
	tkbind(earlyCutoffInfo, "<Button-1>",earlyCutoffHelp)

	earlyValueLabel <- tklabel(analysisFrame,text="Early peak cutoff (degrees C): ")
	earlyValueText <- tklabel(analysisFrame,text=as.character(tclvalue(earlyValue)))
	slider1 <- tkscale(analysisFrame, from=0, to=5,
		       showvalue=F, variable=earlyValue,
		       resolution=0.1, orient="horizontal")
	tkconfigure(earlyValueText,textvariable=earlyValue)
	tkgrid(earlyCutoffInfo,earlyValueLabel,earlyValueText,slider1)

	tkgrid.configure(earlyCutoffInfo,sticky="w")
	tkgrid.configure(earlyValueLabel,sticky="w")
	tkgrid.configure(earlyValueText,sticky="w")
	tkgrid.configure(slider1,sticky="we")


	#option for excluding late peaks
	lateInfo <- tklabel(analysisFrame,text="i",font=fontForHelp)
	tkbind(lateInfo, "<Button-1>",lateHelp)
	lpc <- tkcheckbutton(analysisFrame)
	tkconfigure(lpc,variable=lpcValue)
	lateLab=tklabel(analysisFrame,text="Exclude late peaks")

	tkgrid(lateInfo,lateLab,lpc)

	tkgrid.configure(lateInfo,sticky="w")
	tkgrid.configure(lateLab,sticky="w")
	tkgrid.configure(lpc,sticky="w")

	lateCutoffInfo <- tklabel(analysisFrame,text="i",font=fontForHelp)
	tkbind(lateCutoffInfo, "<Button-1>",lateCutoffHelp)

	lateValueLabel <- tklabel(analysisFrame,text="Late peak cutoff (degrees C): ")
	lateValueText <- tklabel(analysisFrame,text=as.character(tclvalue(lateValue)))
	tkconfigure(lateValueText,textvariable=lateValue)
	slider2 <- tkscale(analysisFrame, from=0, to=5,
		       showvalue=F, variable=lateValue,
		       resolution=0.1, orient="horizontal")

	tkgrid(lateCutoffInfo,lateValueLabel,lateValueText,slider2)

	tkgrid.configure(lateCutoffInfo,sticky="w")
	tkgrid.configure(lateValueLabel,sticky="w")
	tkgrid.configure(lateValueText,sticky="w")
	tkgrid.configure(slider2,sticky="we")

	#option for including shoulders
	shoulderInfo <- tklabel(analysisFrame,text="i",font=fontForHelp)
	tkbind(shoulderInfo, "<Button-1>",shoulderHelp)
	hht <- tkcheckbutton(analysisFrame)
	tkconfigure(hht,variable=hhtValue)
	shoulderLab=tklabel(analysisFrame,text="Include shoulders")

	tkgrid(shoulderInfo,shoulderLab,hht)

	tkgrid.configure(shoulderInfo,sticky="w")
	tkgrid.configure(shoulderLab,sticky="w")
	tkgrid.configure(hht,sticky="w")

	shoulderCutoffInfo <- tklabel(analysisFrame,text="i",font=fontForHelp)
	tkbind(shoulderCutoffInfo, "<Button-1>",shoulderCutoffHelp)

	shoulderValueLabel <- tklabel(analysisFrame,text="Shoulder height threshold (percent): ")
	shoulderValueText <- tklabel(analysisFrame,text=as.character(tclvalue(shoulderValue)))
	tkconfigure(shoulderValueText,textvariable=shoulderValue)
	slider3 <- tkscale(analysisFrame, from=0, to=100,
		       showvalue=F, variable=shoulderValue,
		       resolution=5, orient="horizontal")
	tkgrid(shoulderCutoffInfo,shoulderValueLabel,shoulderValueText,slider3)

	tkgrid.configure(shoulderCutoffInfo,sticky="w")
	tkgrid.configure(shoulderValueLabel,sticky="w")
	tkgrid.configure(shoulderValueText,sticky="w")
	tkgrid.configure(slider3,sticky="we")

	t1Info <- tklabel(analysisFrame,text="i",font=fontForHelp)
	tkbind(t1Info, "<Button-1>",t1Help)

	t1ValueLabel <- tklabel(analysisFrame,text="Theta 1 (degrees): ")
	t1ValueText <- tklabel(analysisFrame,text=as.character(tclvalue(t1Value)))
	tkconfigure(t1ValueText,textvariable=t1Value)

	slider4 <- tkscale(analysisFrame, from=0, to=90,
		       showvalue=F, variable=t1Value,
		       resolution=5, orient="horizontal")
	tkgrid(t1Info,t1ValueLabel,t1ValueText,slider4)

	tkgrid.configure(t1Info,sticky="w")
	tkgrid.configure(t1ValueLabel,sticky="w")
	tkgrid.configure(t1ValueText,sticky="w")
	tkgrid.configure(slider4,sticky="we")

	t1HoldInfo <- tklabel(analysisFrame,text="i",font=fontForHelp)
	tkbind(t1HoldInfo, "<Button-1>",t1HoldHelp)

	t1HoldValueLabel <- tklabel(analysisFrame,text="T1 hold interval (degrees C): ")
	t1HoldValueText <- tklabel(analysisFrame,text=as.character(tclvalue(t1HoldValue)))
	tkconfigure(t1HoldValueText,textvariable=t1HoldValue)

	slider5 <- tkscale(analysisFrame, from=0, to=3,
		       showvalue=F, variable=t1HoldValue,
		       resolution=0.1, orient="horizontal")

	tkgrid(t1HoldInfo,t1HoldValueLabel,t1HoldValueText,slider5)

	tkgrid.configure(t1HoldInfo,sticky="w")
	tkgrid.configure(t1HoldValueLabel,sticky="w")
	tkgrid.configure(t1HoldValueText,sticky="w")
	tkgrid.configure(slider5,sticky="we")

	t1SlopeWinInfo <- tklabel(analysisFrame,text="i",font=fontForHelp)
	tkbind(t1SlopeWinInfo, "<Button-1>",t1SlopeWinHelp)

	t1SlopeWinValueLabel <- tklabel(analysisFrame,text="T1 slope window (degrees C): ")
	t1SlopeWinValueText <- tklabel(analysisFrame,text=as.character(tclvalue(t1SlopeWinValue)))
	tkconfigure(t1SlopeWinValueText,textvariable=t1SlopeWinValue)

	slider6 <- tkscale(analysisFrame, from=0, to=2,
		       showvalue=F, variable=t1SlopeWinValue,
		       resolution=0.1, orient="horizontal")

	tkgrid(t1SlopeWinInfo,t1SlopeWinValueLabel,t1SlopeWinValueText,slider6)

	tkgrid.configure(t1SlopeWinInfo,sticky="w")
	tkgrid.configure(t1SlopeWinValueLabel,sticky="w")
	tkgrid.configure(t1SlopeWinValueText,sticky="w")
	tkgrid.configure(slider6,sticky="we")

	t2Info <- tklabel(analysisFrame,text="i",font=fontForHelp)
	tkbind(t2Info, "<Button-1>",t2Help)

	t2ValueLabel <- tklabel(analysisFrame,text="Theta 2 (degrees): ")
	t2ValueText <- tklabel(analysisFrame,text=as.character(tclvalue(t2Value)))
	tkconfigure(t2ValueText,textvariable=t2Value)
	slider7 <- tkscale(analysisFrame, from=0, to=90,
		       showvalue=F, variable=t2Value,
		       resolution=5, orient="horizontal")

	tkgrid(t2Info,t2ValueLabel,t2ValueText,slider7)
	tkgrid.configure(t2Info,sticky="w")
	tkgrid.configure(t2ValueLabel,sticky="w")
	tkgrid.configure(t2ValueText,sticky="w")
	tkgrid.configure(slider7,sticky="we")

	t2HoldInfo <- tklabel(analysisFrame,text="i",font=fontForHelp)
	tkbind(t2HoldInfo, "<Button-1>",t2HoldHelp)

	t2HoldValueLabel <- tklabel(analysisFrame,text="T2 hold interval (degrees C): ")
	t2HoldValueText <- tklabel(analysisFrame,text=as.character(tclvalue(t2HoldValue)))
	tkconfigure(t2HoldValueText,textvariable=t2HoldValue)

	slider8 <- tkscale(analysisFrame, from=0, to=3,
		       showvalue=F, variable=t2HoldValue,
		       resolution=0.1, orient="horizontal")

	tkgrid(t2HoldInfo,t2HoldValueLabel,t2HoldValueText,slider8)

	tkgrid.configure(t2HoldInfo,sticky="w")
	tkgrid.configure(t2HoldValueLabel,sticky="w")
	tkgrid.configure(t2HoldValueText,sticky="w")
	tkgrid.configure(slider8,sticky="we")

	t2SlopeWinInfo <- tklabel(analysisFrame,text="i",font=fontForHelp)
	tkbind(t2SlopeWinInfo, "<Button-1>",t2SlopeWinHelp)

	t2SlopeWinValueLabel <- tklabel(analysisFrame,text="T2 slope window (degrees C): ")
	t2SlopeWinValueText <- tklabel(analysisFrame,text=as.character(tclvalue(t2SlopeWinValue)))
	tkconfigure(t2SlopeWinValueText,textvariable=t2SlopeWinValue)

	slider9 <- tkscale(analysisFrame, from=0, to=2,
		       showvalue=F, variable=t2SlopeWinValue,
		       resolution=0.1, orient="horizontal")

	tkgrid(t2SlopeWinInfo,t2SlopeWinValueLabel,t2SlopeWinValueText,slider9)

	tkgrid.configure(t2SlopeWinInfo,sticky="w")
	tkgrid.configure(t2SlopeWinValueLabel,sticky="w")
	tkgrid.configure(t2SlopeWinValueText,sticky="w")
	tkgrid.configure(slider9,sticky="we")

	scoreInfo <- tklabel(analysisFrame,text="i",font=fontForHelp)
	tkbind(scoreInfo, "<Button-1>",scoreHelp)
	score.but <- tkbutton(analysisFrame,text="HRM score file (opt)",command=selectScoreFile)
	scoreValueText <- tkentry(analysisFrame,width="40",textvariable=scoreName)
	tkconfigure(scoreValueText,textvariable=scoreName)
	tkgrid(scoreInfo,score.but,scoreValueText)
	tkgrid.configure(scoreInfo,sticky="w")
	tkgrid.configure(score.but,sticky="we")
	tkgrid.configure(scoreValueText,sticky="we")

	#option for including Lasso model
	modelInfo <- tklabel(analysisFrame,text="i",font=fontForHelp)
	tkbind(modelInfo, "<Button-1>",modelHelp)
	lasso <- tkcheckbutton(analysisFrame)
	tkconfigure(lasso,variable=lassoValue)
	modelLab = tklabel(analysisFrame,text="Use Lasso regression")

	tkgrid(modelInfo,modelLab,lasso)
	tkgrid.configure(modelInfo,sticky="w")
	tkgrid.configure(modelLab,sticky="w")
	tkgrid.configure(lasso,sticky="w")

	modelFileInfo <- tklabel(analysisFrame,text="i",font=fontForHelp)
	tkbind(modelFileInfo, "<Button-1>",modelFileHelp)
	model.but <- tkbutton(analysisFrame,text="Lasso model file (opt)",command=selectModelFile)
	modelValueText <- tkentry(analysisFrame,width="40",textvariable=modelName)
	tkconfigure(modelValueText,textvariable=modelName)

	tkgrid(modelFileInfo,model.but,modelValueText)
	tkgrid.configure(modelFileInfo,sticky="w")
	tkgrid.configure(model.but,sticky="we")
	tkgrid.configure(modelValueText,sticky="we")


	tkgrid(analysisFrame)
	tkgrid.configure(analysisFrame,sticky="nswe")

	done <- tclVar(0)   # tclVar() creates a Tcl variable

	OK.but <- tkbutton(dlg, text = "OK",
	    command = function() tclvalue(done) <- 1)
	tkgrid(OK.but)
	tkgrid.configure(OK.but,sticky="nswe")
	tkwait.variable(done)
        tkgrab.release(dlg)
        tkdestroy(dlg)
        tkfocus(tt)
    }

    outputDlg <- function()
    {
	header="Output Options"
	dlg <- tktoplevel()
	tkwm.deiconify(dlg)
	tkgrab.set(dlg)
	tkfocus(dlg)
	title = "HRM Diversity Assay Analysis Tool"
	tkwm.title(dlg, title)
	#tkwm.geometry(dlg,"400x200+0+0")
	tkwm.geometry(dlg,"+0+0")

        frameOutput <- tkwidget(dlg, "labelframe",  text=header)
	tkgrid.columnconfigure(frameOutput,1,weight=1)
	tkgrid.columnconfigure(frameOutput,2,weight=1)

	outdirInfo <- tklabel(frameOutput,text="i",font=fontForHelp)
	outdir.but <- tkbutton(frameOutput,text="Output directory",command=selectOutputDir)
	outdirValueText <- tkentry(frameOutput,width="40",textvariable=outdirName)
	tkbind(outdirInfo, "<Button-1>",outdirHelp)

	tkgrid(outdirInfo,outdir.but,outdirValueText)
	tkgrid.configure(outdirInfo,sticky="w")
	tkgrid.configure(outdir.but,sticky="we")
	tkgrid.configure(outdirValueText,sticky="we")

	outfileInfo <- tklabel(frameOutput,text="i",font=fontForHelp)
	tkbind(outfileInfo, "<Button-1>",outfileHelp)
	outfileLab <- tklabel(frameOutput,text="Output file name (opt)")
	outfileValueText <- tkentry(frameOutput,width="40",textvariable=outfileName)
	tkconfigure(outfileValueText,textvariable=outfileName)

	tkgrid(outfileInfo,outfileLab,outfileValueText)
	tkgrid.configure(outfileInfo,sticky="w")
	tkgrid.configure(outfileLab,sticky="w")
	tkgrid.configure(outfileValueText,sticky="we")

	statsInfo <- tklabel(frameOutput,text="i",font=fontForHelp)
	tkbind(statsInfo, "<Button-1>",statsHelp)
	statsLab <- tklabel(frameOutput,text="Stats output file name (opt)")
	statsValueText <- tkentry(frameOutput,width="40",textvariable=statsName)
	tkconfigure(statsValueText,textvariable=statsName)

	tkgrid(statsInfo,statsLab,statsValueText)
	tkgrid.configure(statsInfo,sticky="w")
	tkgrid.configure(statsLab,sticky="w")
	tkgrid.configure(statsValueText,sticky="we")

	tkgrid(frameOutput,sticky="nswe")

	done <- tclVar(0)   # tclVar() creates a Tcl variable

	OK.but <- tkbutton(dlg, text = "OK",
	    command = function() tclvalue(done) <- 1)
	tkgrid(OK.but)
	tkgrid.configure(OK.but,sticky="nswe")
	tkwait.variable(done)
        tkgrab.release(dlg)
        tkdestroy(dlg)
        tkfocus(tt)
    }


    viewCommand <- function()
    {
	tkconfigure(tt,cursor="watch")
	tcl("update", "idletasks")

	outdirN <- as.character(tclvalue(outdirName))
	finf <- file.info(list.files(path = outdirN, pattern = "*.*"))
	newones <- finf[difftime(launchTime, finf[,"mtime"], units="secs")< 0,]
	print(rownames(newones))

	tkconfigure(tt,cursor="arrow")
	tcl("update", "idletasks")

	if (nrow(newones) > 0)
	{
	    if(.Platform$OS.type == "unix")
		header="Double click a file to view its content"
	    else
		header="Newly created files"

	    dlg <- tktoplevel()
	    tkwm.deiconify(dlg)
	    tkgrab.set(dlg)
	    tkfocus(dlg)
	    title = "HRM Diversity Assay Analysis Tool"
	    tkwm.title(dlg, title)
	    #tkwm.geometry(dlg,"400x200+0+0")
	    tkwm.geometry(dlg,"+0+0")

	    frameOutput <- tkwidget(dlg, "labelframe",  text=header)
	    tkgrid.columnconfigure(frameOutput,1,weight=1)
	    tkgrid.columnconfigure(frameOutput,2,weight=1)
	    tkgrid.configure(frameOutput,sticky="nswe")

	    maxw=max(nchar(rownames(newones)))
	    print(maxw)

	    if (nrow(newones) > 10)
	    {
		scr <- tkscrollbar(frameOutput, repeatinterval=5, command=function(...) tkyview(tl,...))
		tl<-tklistbox(frameOutput,height=10,width=maxw,selectmode="single", yscrollcommand=function(...) tkset(scr,...),background="white")
		tkgrid(tl,scr)
		tkgrid.configure(scr,rowspan=10,sticky="nsw")
		tkgrid.configure(tl,sticky="nswe")
	    }
	    else
	    {
		tl<-tklistbox(frameOutput,height=nrow(newones),width=maxw,selectmode="single",background="white")
		tkgrid(tl)
		tkgrid.configure(tl,sticky="nswe")
	    }

	    for (i in (1:nrow(newones)))
	    {
		tkinsert(tl,"end",rownames(newones)[i])
	    }
	    # setup the callbacks
	    if(.Platform$OS.type == "unix") {
		pdfV <- as.character(tclvalue(pdfVName))
		pngV <- as.character(tclvalue(pngVName))

		tkbind(tl, "<Double-Button-1>", function(...)
		{
		    fileN=rownames(newones)[as.numeric(tkcurselection(tl))+1]
		    suffix = gsub (pattern=".*\\.",replacement="",fileN)
		    if (suffix == "pdf")
			system(paste(pdfV,fileN),wait=FALSE)
		    else if (suffix == "png")
			system(paste(pngV,fileN),wait=FALSE)
		    else
			tkmessageBox(message=paste("No default viewer for file ",fileN,".",sep=""))
		})
	    }
	    done <- tclVar(0)   # tclVar() creates a Tcl variable

	    OK.but <- tkbutton(dlg, text = "Done",
		command = function() tclvalue(done) <- 1)
	    tkgrid(OK.but)
	    tkgrid.configure(OK.but,sticky="nswe")
	    tkwait.variable(done)
	    tkgrab.release(dlg)
	    tkdestroy(dlg)
	    tkfocus(tt)
	}
	else
	    tkmessageBox(message=paste("Found ",nrow(newones)," new files.",sep=""))
    }

    loadCommand <- function()
    {
	settingsFile<-tclvalue(tkgetOpenFile(filetypes="{{.txt Files} {.txt}} {{All files} *}"))
	if (!nchar(settingsFile))
	    return

	loadSettings(settingsFile)
    }

    saveCommand <- function()
    {
	settingsFile<-tclvalue(tkgetSaveFile(filetypes="{{.txt Files} {.txt}} {{All files} *}"))

	if (!nchar(settingsFile))
	    return

	write (c("sample",as.character(tclvalue(sampleName))),file=settingsFile,ncolumns=2,sep="=",append=FALSE)
	write (c("plottype",as.character(tclvalue(plotTypeValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("pdfViewer",as.character(tclvalue(pdfVName))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("pngViewer",as.character(tclvalue(pngVName))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("smoothFluor",as.character(tclvalue(fluorSmoothCBValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("smoothdFluor",as.character(tclvalue(dfluorSmoothCBValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("excludeEarly",as.character(tclvalue(epcValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("earlyCutoff",as.character(tclvalue(earlyValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("excludeLate",as.character(tclvalue(lpcValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("lateCutoff",as.character(tclvalue(lateValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("includeShoulders",as.character(tclvalue(hhtValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("shoulderHeight",as.character(tclvalue(shoulderValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("useLasso",as.character(tclvalue(lassoValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("lassoModel",as.character(tclvalue(modelName))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("theta1",as.character(tclvalue(t1Value))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("theta1Hold",as.character(tclvalue(t1HoldValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("theta1Window",as.character(tclvalue(t1SlopeWinValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("theta2",as.character(tclvalue(t2Value))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("theta2Hold",as.character(tclvalue(t2HoldValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("theta2Window",as.character(tclvalue(t2SlopeWinValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("useXrange",as.character(tclvalue(xrangeValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("minX",as.character(tclvalue(xminValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	write (c("maxX",as.character(tclvalue(xmaxValue))),file=settingsFile,ncolumns=2,sep="=",append=TRUE)
	tkmessageBox(title="Info",message=paste("Saved all settings to file ",settingsFile),icon="info",type="ok")
    }

    exitCommand <- function()
    {
	title = "HRM Diversity Assay Analysis Tool"
	ReturnVal <- tkmessageBox(title=title,message="Are you sure you want to exit?",icon="question",type="okcancel",default="ok")
	if (as.character(ReturnVal) == "ok")
	    tkdestroy(tt)
    }

    selectOutputDir <- function()
    {
	outdirN<-tclvalue(tkchooseDirectory())
	if (!nchar(outdirN))
	{
	    return
	}
	#else
	#    tkmessageBox(message=paste("The output directory selected was ",outdirN))

	tclvalue(outdirName)<<-outdirN
	setwd(outdirN)
    }
    selectDir <- function()
    {
	dirN<-tclvalue(tkchooseDirectory())
	if (!nchar(dirN))
	{
	    return
	}
	#else
	#    tkmessageBox(message=paste("The input directory selected was ",dirN))

	tclvalue(dirName)<<-dirN
	tclvalue(fileName)<<-""
    }

    selectFile <- function()
    {
	fileN<-tclvalue(tkgetOpenFile(filetypes="{{ABT Files} {.ABT}} {{All files} *}"))
	if (!nchar(fileN))
	{
	    return
	}
	#else
	#    tkmessageBox(message=paste("The file selected was ",fileN))

	tclvalue(fileName)<<-gsub (pattern=".*/",replacement="",fileN)
	if (fileN != "")
	    tclvalue(dirName)<<-gsub (pattern="/[^/]*$",replacement="",fileN)
    }

    selectScoreFile <- function()
    {
	scorefileN<-tclvalue(tkgetOpenFile(filetypes="{{csv Files} {.csv}} {{All files} *}"))
	if (!nchar(scorefileN))
	{
	    return
	}
	#else
	#    tkmessageBox(message=paste("The score file selected was ",scorefileN))

	tclvalue(scoreName)<<-scorefileN
    }

    selectModelFile <- function()
    {
	modelfileN<-tclvalue(tkgetOpenFile(filetypes="{{rda Files} {.rda}} {{All files} *}"))
	if (!nchar(modelfileN))
	{
	    return
	}
	#else
	#    tkmessageBox(message=paste("The model file selected was ",modelfileN))

	tclvalue(modelName)<<-modelfileN
    }

    helpDialog <- function(title,helpText)
    {
	ReturnVal <- tkmessageBox(title=title,message=helpText,icon="info",type="ok")
	return(ReturnVal)
    }
    outdirHelp <-function()
    {
	return(helpDialog("Output Directory","The \"Output directory\" field specifies the directory that should be used to store the output files generated during the analysis. This field also specifies the directory that the \"View New Files\" function should search in order to find and display the recently generated files."))
    }

    dirHelp <-function()
    {
	return(helpDialog("Input Directory","The \"Input directory\" field specifies the directory containing the ABT files to be analyzed.  Setting this clears the \"Input file\" field."))
    }

    fileHelp <-function()
    {
	return(helpDialog("Input File","The \"Input file\" field specifies a particular ABT file to be analyzed.  Setting this also sets the \"Input directory\"."))
    }

    sampleHelp <-function()
    {
	return(helpDialog("Sample Name","The \"Sample name\" field specifies a particular sample or subset of samples for analysis. If a complete sample designation is supplied, the analysis will be restricted to a single sample.  If a text string that is common to multiple samples is supplied, all samples whose names contain the supplied string will be subject to analysis."))
    }
    plotHelp <-function()
    {
	return(helpDialog("Plotting","The user may elect to generate summary plots in the PNG or PDF format.  PDF plots also present the T1, T2, and an HRM score for each sample.  The name assigned to each plot file is derived from the data file (for PDFs) or the sample name (for PNGs).  When subsets of the samples are chosen (using the \"Sample name\" field), each file will bear the sample name."))
    }
    xrangeHelp <-function()
    {
	return(helpDialog("X-axis Range","This option allows the user to override the default ranges for temperature used in the Fluorescence vs. Temperature and dFluorescence/dT vs. Temperature plots."))
    }

    pdfVHelp <-function()
    {
	return(helpDialog("PDF viewer program","The \"PDF viewer\" option specifies the program to launch to view a PDF file (ex. acroread).  The path must be included."))
    }

    pngVHelp <-function()
    {
	return(helpDialog("PNG viewer program","The \"PNG viewer\" option specifies the program to launch to view a PNG file (ex. display).  The path must be included."))
    }

    dfluorSmoothHelp <-function()
    {
	return(helpDialog("Show dFluorescence Smoothing","This option causes a local regression smoothing curve to be included when -dFluorescence/dT values are plotted against Temperature.  The smoothed curve is used in the analysis regardless of whether this line is displayed."))
    }

    fluorSmoothHelp <-function()
    {
	return(helpDialog("Show Fluorescence Smoothing","This option causes a local regression smoothing curve to be included when Fluorescence values are plotted against Temperature."))
    }

    earlyHelp <-function()
    {
	return(helpDialog("Exclude Early Peaks","This option causes early peaks preceding the main peak to be excluded from analysis.  Early peaks are excluded when they are separated from the principal peak by more than the number of degrees that were specified using the \"Early peak cutoff\" slider."))
    }

    lateHelp <-function()
    {
	return(helpDialog("Exclude Late Peaks","This option causes late peaks following the main peak to be excluded from analysis.  Late peaks are excluded when they are separated from the principal peak by more than the number of degrees that were specified using the \"Early peak cutoff\" slider."))
    }

    shoulderHelp <-function()
    {
	return(helpDialog("Include Shoulders","This option causes shoulders on the sides of the melting peaks to be included in the analysis. Shoulder inclusion thresholds may be set as a percentage of principal peak height using the \"Shoulder height threshold\" slider below."))
    }

    earlyCutoffHelp <-function()
    {
	return(helpDialog("Early Peak Cutoff","The slider can be used to set the cutoff number of degrees Celsius that must separate the early peak from the principal peak before the early peak is excluded."))
    }

    lateCutoffHelp <-function()
    {
	return(helpDialog("Late Peak Cutoff","The slider can be used to set the cutoff number of degrees Celsius that must separate the late peak from the principal peak before the late peak is excluded."))
    }

    shoulderCutoffHelp <-function()
    {
	return(helpDialog("Shoulder Height Threshold","The slider can be used to set the threshold as a percentage of the principal peak height that a shoulder must exceed before the shoulder is included in the HRM score calculation."))
    }

    t1Help <-function()
    {
	return(helpDialog("Theta 1","The slider can be used to set the number of degrees that the angle between the baseline and the tangent to the melting peak (-dFluorescence/dt vs. T curve) must exceed for a point to be identified as a T1 candidate temperature. The number of degrees Celsius over which this angle must be met or exceeded can be set using the \"T1 Hold Interval\" slider."))
    }

    t2Help <-function()
    {
	return(helpDialog("Theta 2","The slider can be used to set the number of degrees that the angle between the baseline and the tangent to the melting peak (-dFluorescence/dt vs. T curve) must exceed for a point to be identified as a T2 candidate temperature. The number of degrees Celsius over which this angle must be met or exceeded can be set using the \"T2 Hold Interval\" slider."))
    }

    t1HoldHelp <-function()
    {
	return(helpDialog("T1 Hold Interval","The slider can be used to set the number of degrees Celsius over which the angle specified above [using the \"Theta 1 (degrees)\" slider] must be met or exceeded for a point to be identified as a T1 candidate temperature."))
    }

    t2HoldHelp <-function()
    {
	return(helpDialog("T2 Hold Interval","The slider can be used to set the number of degrees Celsius over which the angle specified above [using the \"Theta 2 (degrees)\" slider] must be met or exceeded for a point to be identified as a T2 candidate temperature."))
    }

    t1SlopeWinHelp <-function()
    {
	return(helpDialog("T1 Slope Window","The slider can be used to set the number of degrees Celsius over which the slope is measured in an effort to determine the T1 candidate temperature."))
    }

    t2SlopeWinHelp <-function()
    {
	return(helpDialog("T2 Slope Window","The slider can be used to set the number of degrees Celsius over which the slope is measured in an effort to determine the T2 candidate temperature."))
    }

    statsHelp <-function()
    {
	return(helpDialog("Stats Output File Name","The \"Stats output file name\" field specifies a name for the text file that will contain the key sample characteristics (mean fluorescence, max fluorescence, etc.)."))
    }

    outfileHelp <-function()
    {
	return(helpDialog("Output File Name","The \"Output file name\" field specifies a name for the text file that will contain the analysis results (T1, T2, HRM score, etc.)."))
    }

    scoreHelp <-function()
    {
	return(helpDialog("DivMelt Score File","This option allows the user to specify a file of HRM scores to compare against on a sample by sample basis."))
    }

    modelHelp <-function()
    {
	return(helpDialog("Lasso Model Filtering","This option causes each sample to be evaluated against a Lasso model (included in the package) to determine if samples appear to be valid."))
    }


    modelFileHelp <-function()
    {
	return(helpDialog("Lasso Model File","This option allows the user to specify an alternate Lasso model.  The model must be stored as an R object file."))
    }


    header="File Menu"
    fileFrame <- tkwidget(tt, "labelframe",  text=header, borderwidth=4)

    inputOpts.button <- tkbutton(fileFrame, text = "Input Options...", command = inputDlg)
    tkgrid(inputOpts.button,sticky="we")

    outputOpts.button <- tkbutton(fileFrame, text = "Output Options...", command = outputDlg)
    tkgrid(outputOpts.button,sticky="we")

    tkgrid(fileFrame,sticky="we")
    tkgrid.columnconfigure(fileFrame,0,weight=1)
    tkgrid.rowconfigure(fileFrame,0,weight=1)

    header="Settings Menu"
    settingsFrame <- tkwidget(tt, "labelframe",  text=header, borderwidth=4)

    plottingOpts.button <- tkbutton(settingsFrame, text = "Plotting Options...", command = plottingDlg)
    tkgrid(plottingOpts.button,sticky="we")

    analysisOpts.button <- tkbutton(settingsFrame, text = "Analysis Options...", command = analysisDlg)
    tkgrid(analysisOpts.button,sticky="we")

    Load.button <- tkbutton(settingsFrame,text="Load Settings...",command=loadCommand)
    tkgrid(Load.button,sticky="we")

    Save.button <- tkbutton(settingsFrame,text="Save Settings...",command=saveCommand)
    tkgrid(Save.button,sticky="we")

    Show.button <- tkbutton(settingsFrame,text="Show Settings...",command=showCommand)
    tkgrid(Show.button,sticky="we")

    tkgrid(settingsFrame,sticky="we")
    tkgrid.columnconfigure(settingsFrame,0,weight=1)
    tkgrid.rowconfigure(settingsFrame,0,weight=1)

    header="Main Menu"
    mainFrame <- tkwidget(tt, "labelframe",  text=header, borderwidth=4)

    Go.button <- tkbutton(mainFrame,text="Run Analysis...",command=runCommand)
    tkgrid(Go.button,sticky="we")

    View.button <- tkbutton(mainFrame,text="View New Files...",command=viewCommand)
    tkgrid(View.button,sticky="we")

    Exit.button <- tkbutton(mainFrame,text="Exit",command=exitCommand)
    tkgrid(Exit.button,sticky="we")

    tkgrid(mainFrame,sticky="we")
    tkgrid.columnconfigure(mainFrame,0,weight=1)
    tkgrid.rowconfigure(mainFrame,0,weight=1)


    tkgrid.columnconfigure(tt,0,weight=1)
    tkgrid.rowconfigure(tt,0,weight=1)

    tkfocus(tt)
}
