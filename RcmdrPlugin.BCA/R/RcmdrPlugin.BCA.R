# Created on 26Sep2010 by Dan Putler
# Last modified on 02Sep2014 by Dan Putler

if (getRversion() >= "2.15.1")  
    globalVariables(c('.grab.focus', 'boxplotsVariable', 'buttonsFrame',
	'criteriaFrame', 'criteriaVariable', 'diagonalFrame', 'diagonalVariable',
	'directionVariable', 'distanceTypeFrame', 'distanceTypeVariable', 
	'entry1', 'formulaFrame', 'generalizedLinearModel', 'GenerializedLinearModel',
	'groupsFrame', 'hclustSummary', 'hierarchicalCluster', 'identifyVariable',
	'jitterXVariable', 'jitterYVariable', 'kmeansClustering', 'leafVariable',
	'lhsEntry', 'lhsVariable', 'logXVariable', 'logYVariable', 'lsLineVariable',
	'methodFrame', 'methodVariable', 'onHelp', 'optionsFrame', 'outerOperatorsFrame',
	'psprintf', 'relableFactor', 'rhsVariable', 'scatterPlot', 'scatterPlotMatrix',
	'smoothLineVariable', 'spreadVariable', 'subButtonsFrame', 'subdialog',
	'subsetFrame', 'subsetVariable', 'top', 'typeFrame', 'typeVariable', 'xBox'))

# .onAttach function to place .subset, .targetLevel, and .trueResp into the
# Rcmdr environment. Hopefully this works. This is done to deal with R2.14.0
.onAttach <- function(libname, pkgname){
    if (!interactive()) return()
    Rcmdr <- options()$Rcmdr
    plugins <- Rcmdr$plugins
    if ((!pkgname %in% plugins) && !getRcmdr("autoRestart")) {
        Rcmdr$plugins <- c(plugins, pkgname)
        options(Rcmdr=Rcmdr)
        closeCommander(ask=FALSE, ask.save=TRUE)
        putRcmdr(".subset", gettextRcmdr("<all valid cases>"))
        putRcmdr(".Sample", NULL)
        putRcmdr(".targetLevel", NULL)
        putRcmdr(".trueResp", NULL)
        Commander()
        }
    palette(c("#0072B2","#D55E00","#000000","#F0E442","#2B9F78","#E69F00","#5684E9","#CC79A7"))
	invisible(dev.off())
    }

# The subsetBox with the subset expression as an Rcmdr environment variable
subsetBoxBCA <- defmacro(window=top, model=FALSE, 
    expr={
            subsetVariable <- if (model){
                if (currentModel && currentFields$subset != "") 
                    tclVar(currentFields$subset) else tclVar(getRcmdr(".subset"))
                }
            else tclVar(getRcmdr(".subset"))
            subsetFrame <- tkframe(window)
            subsetEntry <- tkentry(subsetFrame, width="20", textvariable=subsetVariable)
            subsetScroll <- tkscrollbar(subsetFrame, orient="horizontal",
                repeatinterval=5, command=function(...) tkxview(subsetEntry, ...))
            tkconfigure(subsetEntry, xscrollcommand=function(...) tkset(subsetScroll, ...))
            tkgrid(tklabel(subsetFrame, text=gettextRcmdr("Subset expression"), fg="blue"), sticky="w")
            tkgrid(subsetEntry, sticky="w")
            tkgrid(subsetScroll, sticky="ew")
            })


# DBF IMPORT/EXPORT FUNCTIONS. The data input options should expand to
# allow for other DBMS options
importDBF <- function() {
    initializeDialog(title=gettextRcmdr("Import dBase (dbf) Data Set"))
    dsname <- tclVar(gettextRcmdr("Dataset"))
    entryDsname <- tkentry(top, width="20", textvariable=dsname)
    asFactor <- tclVar("0")
    asFactorCheckBox <- tkcheckbutton(top, variable=asFactor)
    onOK <- function(){
        closeDialog()
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == ""){
            errorCondition(recall=importDBF,
                message=gettextRcmdr("You must enter the name of a data set."))
                return()
                }
        if (!is.valid.name(dsnameValue)){
            errorCondition(recall=importDBF, message=
              sprintf(gettextRcmdr('"%s" is not a valid name.'), dsnameValue))
	      return()
            }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
                importDBF()
                return()
                }
            }
        file <- tclvalue(tkgetOpenFile(
            filetypes=gettextRcmdr('{"dBase datasets" {".dbf" ".DBF"}} {"All Files" {"*"}}')))
        if (file == "") {
            tkfocus(CommanderWindow())
            return()
            }
        factor <- tclvalue(asFactor) == "1"
        command <- paste('read.dbf("', file,'", as.is=', factor,')', sep="")
		## The command below is not for this case, but it the way to avoid the use of an assign() call to .GlobalEnv
 	    doItAndPrint(paste(dsnameValue, "<-", command))
		activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="read.dbf")
    tkgrid(tklabel(top, text=gettextRcmdr("Enter name for data set:")), entryDsname, sticky="w")
    tkgrid(tklabel(top, text=gettextRcmdr("Do not convert character\nvariables to factors"), justify="left"), 
        asFactorCheckBox, sticky="w")
    tkgrid(buttonsFrame, columnspan="2", sticky="w")
    tkgrid.configure(entryDsname, sticky="w")
    tkgrid.configure(asFactorCheckBox, sticky="w")
    dialogSuffix(rows=4, columns=2, focus=entryDsname)
    } 

exportDBF <- function() {
    dsname <- activeDataSet()
    saveFile <- tclvalue(tkgetSaveFile(filetypes=gettextRcmdr('{"DBF Files" {".dbf" ".DBF"}} {"All Files" {"*"}}'),
      defaultextension="dbf", initialfile=paste(dsname, ".dbf", sep="")))
    if (saveFile == "") {
        tkfocus(CommanderWindow())
        return()
        }
    command <- paste("write.dbf(", dsname, ', "', saveFile, '")', sep="")
    justDoIt(command)
    logger(command)
    Message(paste(gettextRcmdr("Active dataset exported to file"), saveFile), type="note")
    tkfocus(CommanderWindow())
    }

# DATA MENU
createSamples <- function(){
	.activeDataSet <- ActiveDataSet()
    initializeDialog(title=gettextRcmdr("Create Samples"))
    optionsFrame <- tkframe(top)
    estPercent <- tclVar("34")
    estPctSlider <- tkscale(optionsFrame, from=20, to=80, showvalue=TRUE,
      variable=estPercent, resolution=1, orient="horizontal")
    valPercent <- tclVar("33")
    valPctSlider <- tkscale(optionsFrame, from=10, to=70, showvalue=TRUE,
      variable=valPercent, resolution=1, orient="horizontal")
    seedVal <- tclVar("1")
    seedField <- tkentry(optionsFrame, width="3", textvariable=seedVal)
    assignName <- tclVar("Sample")
    assignField <- tkentry(optionsFrame, width="15",
      textvariable=assignName)
    onOK <- function(){
        prctEst <- paste("0.", tclvalue(estPercent), sep="")
        prctVal <- paste("0.", tclvalue(valPercent), sep="")
        if((as.numeric(prctEst) + as.numeric(prctVal)) > 1) {
            errorCondition(recall=createSamples, 
              message=gettextRcmdr("The estimation and validation samples exceed 100%."))
            return()
            }
        randSeed <- trim.blanks(tclvalue(seedVal))
        assignVar <- trim.blanks(tclvalue(assignName))
        if(!is.valid.name(assignVar)) {
            errorCondition(recall=createSamples, message=
	      sprintf(gettextRcmdr('"%s" is not a valid name.'), assignVar))
            return()
            }
        closeDialog()
        command <- paste(.activeDataSet, "$", assignVar," <- create.samples(",
          .activeDataSet, ", est = ", prctEst, ", val = ", prctVal,
          ", rand.seed = ", randSeed, ")", sep="")
        justDoIt(command)
        logger(command)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="ls")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Estimation sample percent:")),
      estPctSlider, sticky="sw")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Validation sample percent:")),
      valPctSlider, sticky="sw")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Random seed: ")),
      seedField, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Assignment variable: ")),
      assignField, sticky="w")
    tkgrid(optionsFrame, sticky="nw")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=2, columns=2)
    }

variableSummary <- function() {
    command <- paste("variable.summary(", ActiveDataSet(), ")", sep="")
    doItAndPrint(command)
    }

relabelFactor <- function(){
	.activeDataSet <- ActiveDataSet()
    initializeDialog(title=gettextRcmdr("Relabel a Factor"))
    variableBox <- variableListBox(top, Factors(), title=gettextRcmdr("Factor (pick one)"))
    factorName <- tclVar(gettextRcmdr("<same as original>"))
    factorNameField <- tkentry(top, width="20", textvariable=factorName)
    onOK <- function(){
        variable <- getSelection(variableBox)
        closeDialog()
        if (length(variable) == 0) {
            errorCondition(recall=relabelFactor, 
              message=gettextRcmdr("You must select a variable."))
            return()
            }
        name <- trim.blanks(tclvalue(factorName))
        if (name == gettextRcmdr("<same as original>")) name <- variable
        if (!is.valid.name(name)){
            errorCondition(recall=relabelFactor, message=
	     sprintf(gettextRcmdr('"%s" is not a valid name.'), name))
	     return()
            }
        if (is.element(name, Variables())) {


            if ("no" == tclvalue(checkReplace(name))){
                if (.grab.focus) tkgrab.release(top)
                tkdestroy(top)
                relabelFactor()
                return()
                }
            }
        old.labels <- eval(parse(text=paste("levels(", .activeDataSet, 
          "$", variable, ")", 
            sep="")), envir=.GlobalEnv)
        nvalues <- length(old.labels)
        if (nvalues > 30) {
            errorCondition(recall=relabelFactor,
                message=psprintf(gettextRcmdr("Number of levels (%d) is too large."),
                nvalues))
            return()
            }
        initializeDialog(subdialog, title=gettextRcmdr("New Labels"))
        labels <- rep("", nvalues)
        onOKsub <- function() {
            closeDialog(subdialog)
            for (i in 1:nvalues){
                labels[i] <-
      eval(parse(text=paste("tclvalue(levelLabel", i, ")", sep="")))
                }
            if (any(labels == "")){
                errorCondition(recall=relabelFactor,
                  message=gettextRcmdr("All factor levels must be given a label."))
                relableFactor()
                return()
                }
            command <- paste("relabel.factor(", .activeDataSet, "$", variable,
              ", new.labels=c(", paste(paste("'", labels, "'", sep=""),
              collapse=","), "))", sep="")
            justDoIt(paste(.activeDataSet, "$", name, " <- ", command, sep=""))
            logger(paste(.activeDataSet,"$", name," <- ", command, sep=""))
            activeDataSet(.activeDataSet)
            tkfocus(CommanderWindow())
            }
        subOKCancelHelp()
        tkgrid(tklabel(subdialog, text=gettextRcmdr("Old Labels"), fg="blue"),
            tklabel(subdialog, text=gettextRcmdr("New Labels"), fg="blue"), sticky="w")
        for (i in 1:nvalues){
            valVar <- paste("levelLabel", i, sep="")
            assign(valVar, tclVar(""))
            assign(paste("entry", i, sep=""), tkentry(subdialog, width="12", 
                textvariable=eval(parse(text=valVar))))
            tkgrid(tklabel(subdialog, text=old.labels[i]),
              eval(parse(text=paste("entry", i, sep=""))), sticky="w")
            }
        tkgrid(subButtonsFrame, sticky="w", columnspan=2)
        dialogSuffix(subdialog, focus=entry1, rows=nvalues+1, columns=2)
        }
    OKCancelHelp(helpSubject="relabel.factor")
    tkgrid(getFrame(variableBox), sticky="nw")
    tkgrid(tklabel(top, text=gettextRcmdr("Name for factor")), sticky="w")
    tkgrid(factorNameField, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1, preventGrabFocus=TRUE)
    }

## Clustering
# kmeans menu corrected for the subsetBox
# This function needs to be worked on. It currently does things that the R
# kmeans function already does
kmeansClusteringBCA <- function(){
	.activeDataSet <- ActiveDataSet()
    initializeDialog(title=gettextRcmdr("KMeans Clustering"))
    dataFrame <- tkframe(top)
    xBox <- variableListBox(dataFrame, Numeric(), selectmode="multiple",
      title=gettextRcmdr("Variables (pick one or more)"))
    subsetBoxBCA(dataFrame)
    assignFrame <- tkframe(dataFrame)
    assignClusters <- tclVar("0")
    assignCB <- tkcheckbutton(assignFrame)
    tkconfigure(assignCB, variable=assignClusters)
    assignName <- tclVar("KMeans")
    assignField <- tkentry(assignFrame, width="15",
      textvariable=assignName)
    optionsFrame <- tkframe(top)
    clusterNumber <- tclVar("2")
    clusterNumSlider <- tkscale(optionsFrame, from=2, to=10, showvalue=TRUE,
      variable=clusterNumber, resolution=1, orient="horizontal")
    seedNumber <- tclVar("10")
    seedNumSlider <- tkscale(optionsFrame, from=1, to=20, showvalue=TRUE,
      variable=seedNumber, resolution=1, orient="horizontal")
    iterNumber <- tclVar("10")
    iterNumSlider <- tkscale(optionsFrame, from=5, to=30, showvalue=TRUE,
      variable=iterNumber, resolution=5, orient="horizontal")
    summaryClusters <- tclVar("1")
    summaryCB <- tkcheckbutton(optionsFrame)
    tkconfigure(summaryCB, variable=summaryClusters)
    plotCls2d <- tclVar("1")
    plot2dCB <- tkcheckbutton(optionsFrame)
    tkconfigure(plot2dCB, variable=plotCls2d)
    plotCls3d <- tclVar("0")
    plot3dCB <- tkcheckbutton(optionsFrame)
    tkconfigure(plot3dCB, variable=plotCls3d)
    plotFrame <- tkframe(optionsFrame)
    plotPts <- tclVar("1")
    plotPtsCB <- tkcheckbutton(plotFrame)
    tkconfigure(plotPtsCB, variable=plotPts)
    plotCnt <- tclVar("0")
    plotCntCB <- tkcheckbutton(plotFrame)
    tkconfigure(plotCntCB, variable=plotCnt)
    onOK <- function(){
        x <- getSelection(xBox)
        nvar <- length(x)
        subset <- trim.blanks(tclvalue(subsetVariable))
        nClusters <- tclvalue(clusterNumber)
        seeds <- tclvalue(seedNumber)
        iters <- tclvalue(iterNumber)
        clusterSummary <- tclvalue(summaryClusters)
        clsPlot2d <- tclvalue(plotCls2d)
        clsPlot3d <- tclvalue(plotCls3d)
        ptsPlot <- tclvalue(plotPts)
        cntPlot <- tclvalue(plotCnt)
        clusterAssign <- tclvalue(assignClusters)
        clusterVariable <- trim.blanks(tclvalue(assignName))
        if (clusterAssign == "1" & !is.valid.name(clusterVariable)){
            errorCondition(recall=kmeansClustering, message=
              sprintf(gettextRcmdr('"%s" is not a valid name.'),
	      clusterVariable))
	      return()
            }
        closeDialog()
        if (clusterAssign == "1"){
           if (is.element(clusterVariable, Variables())) {
                if ("no" == tclvalue(checkReplace(clusterVariable))){
                    kmeansClustering()
                    return()
                    }
                }
           } 
        if (length(x)==0) {
            errorCondition(recall=kmeansClustering, 
              message==gettextRcmdr("No variables selected."))
            return()
            }
        varFormula <- paste(x, collapse=" + ")
        vars <- paste(x, collapse=",", sep="")
        .activeDataSet <- ActiveDataSet()
        dset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) .activeDataSet
          else {paste(.activeDataSet, "[", .activeDataSet, "$", subset, ", ]",
            sep="")}
        xmat <- paste("model.matrix(~-1 + ", varFormula, ", ", dset, ")",
          sep="")
        command <- paste("KMeans(", xmat, ", centers = ", nClusters,
          ", iter.max = ", iters, ", num.seeds = ", seeds, ")", sep="")
 	    doItAndPrint(paste(".cluster", "<-", command))
        if (clusterSummary == "1") {
            doItAndPrint(paste(".cluster$size # Cluster Sizes"))
            doItAndPrint(paste(".cluster$centers # Cluster Centroids"))
            doItAndPrint(paste(
              ".cluster$withinss # Within Cluster Sum of Squares"))
            doItAndPrint(paste(
              ".cluster$tot.withinss # Total Within Sum of Squares"))
            doItAndPrint(paste(
              ".cluster$betweenss # Between Cluster Sum of Squares"))
            }
	if((clsPlot2d=="1" | clsPlot3d=="1") & ptsPlot=="0" & cntPlot=="0") {
	   Message(gettextRcmdr("Bi-plot(s) contains no points or centroids."),
	     type="warning")
	   }
	ptsOpt <- "FALSE"
	if(ptsPlot == "1") ptsOpt <- "TRUE"
	cntOpt <- "FALSE"
	if(cntPlot == "1") cntOpt <- "TRUE"
        if (clsPlot2d == "1") {
            plotCmd2d <- paste("bpCent(prcomp(", xmat, 
              "), .cluster$cluster, data.pts = ", ptsOpt, ", centroids = ",
	      cntOpt, ", xlabs = as.character(.cluster$cluster))", sep="")
           justDoIt(plotCmd2d)
           logger(plotCmd2d)
           }
        if (clsPlot3d == "1") {
			if (requireNamespace("rgl", quietly = TRUE)) {
				plotCmd3d <- paste("bpCent3d(prcomp(", xmat, 
				"), .cluster$cluster, data.pts = ", ptsOpt, ", centroids = ",
				cntOpt, ", xlabs = as.character(.cluster$cluster))", sep="")
				justDoIt(plotCmd3d)
				logger(plotCmd3d)
				.Tcl("update")
				activateMenus()
				rgl::rgl.bringtotop()
			} else {
				warning("The rgl package is not available so a three dimensional bi-plot was not created")
			}
		}
        if (clusterAssign == "1") {
            assignCommand <- paste(.activeDataSet, "$", clusterVariable,
              " <- assignCluster(", xmat, ", ", .activeDataSet,
              ", .cluster$cluster)", sep="")
            justDoIt(assignCommand)
            logger(assignCommand)
            activeDataSet(.activeDataSet)
            }
        justDoIt(paste("remove(.cluster)"))
        logger(paste("remove(.cluster)"))
        activateMenus()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="KMeans")
    tkgrid(getFrame(xBox), sticky="nw")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(tklabel(assignFrame, 
      text=gettextRcmdr("Assign clusters to\nthe data set         ")),
      assignCB, sticky="w")
    tkgrid(tklabel(assignFrame, text=gettextRcmdr("Assignment variable: ")),
      assignField, sticky="w")
    tkgrid(assignFrame, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Number of clusters:")),
      clusterNumSlider, sticky="sw")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Number of starting seeds:")),
      seedNumSlider, sticky="sw")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Maximum iterations:")),
      iterNumSlider, sticky="sw")
    tkgrid(tklabel(optionsFrame, 
      text=gettextRcmdr("Print cluster summary")), summaryCB, sticky="w")
    tkgrid(tklabel(optionsFrame, 
      text=gettextRcmdr("2D bi-plot of clusters")), plot2dCB, sticky="w")
    tkgrid(tklabel(optionsFrame, 
      text=gettextRcmdr("3D bi-plot of clusters")), plot3dCB, sticky="w")
    tkgrid(tklabel(plotFrame, text=gettextRcmdr("Plot points")), plotPtsCB,
      tklabel(plotFrame, text=gettextRcmdr("Plot centroids")), plotCntCB,
      sticky="w")
    tkgrid(plotFrame, columnspan=2, sticky="w")
    tkgrid(dataFrame, tklabel(top, text="  "), optionsFrame,
        sticky="nw")
    tkgrid(buttonsFrame, columnspan=3, sticky="w")
    dialogSuffix(rows=3, columns=3)
    }

# The SD Index plot functions for k-means
SDIndexPlot <- function(){
	.activeDataSet <- ActiveDataSet()
    initializeDialog(title=gettextRcmdr("SD Index Plot"))
    dataFrame <- tkframe(top)
    xBox <- variableListBox(dataFrame, Numeric(), selectmode="multiple",
      title=gettextRcmdr("Variables (pick one or more)"))
    subsetBoxBCA(dataFrame)
    optionsFrame <- tkframe(top)
    minCluster <- tclVar("2")
    minClusterSlider <- tkscale(optionsFrame, from=2, to=8, showvalue=TRUE,
      variable=minCluster, resolution=1, orient="horizontal")
    maxCluster <- tclVar("4")
    maxClusterSlider <- tkscale(optionsFrame, from=4, to=10, showvalue=TRUE,
      variable=maxCluster, resolution=1, orient="horizontal")
    iterNumber <- tclVar("10")
    iterNumSlider <- tkscale(optionsFrame, from=5, to=30, showvalue=TRUE,
      variable=iterNumber, resolution=5, orient="horizontal")
    seedNumber <- tclVar("10")
    seedNumSlider <- tkscale(optionsFrame, from=1, to=20, showvalue=TRUE,
      variable=seedNumber, resolution=1, orient="horizontal")
    .activeDataset <- ActiveDataSet()
     onOK <- function(){
        x <- getSelection(xBox)
        nvar <- length(x)
        subset <- tclvalue(subsetVariable)
        minClus <- tclvalue(minCluster)
        maxClus <- tclvalue(maxCluster)
        seeds <- tclvalue(seedNumber)
        iters <- tclvalue(iterNumber)
        closeDialog()
        if (length(x)==0) {
            errorCondition(recall=SDIndexPlot, 
              message=gettextRcmdr("No variables selected."))
            return()
            }
        if (as.numeric(minClus) > as.numeric(maxClus)) {
            errorCondition(recall=SDIndexPlot, message=gettextRcmdr(
              "The minimum number of clusters selected exceeds the maximum."))
            return()
            }
        varFormula <- paste(x, collapse=" + ")
        vars <- paste(x, collapse=",", sep="")
        dset <- if (subset == gettextRcmdr("<all valid cases>")) .activeDataSet
          else {paste(.activeDataSet, "[", .activeDataSet, "$", subset, ", ]",
            sep="")}
        xmat <- paste("model.matrix(~-1 + ", varFormula, ", ", dset, ")",
          sep="")
        command <- paste("SDIndex(", xmat, ", minClus = ", minClus,
           ", maxClus = ", maxClus, ", iter.max = ", iters, 
           ", num.seeds=", seeds, ")", sep="")
        justDoIt(command)
        logger(command)
        activateMenus()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="SDIndex")
    tkgrid(getFrame(xBox), sticky="nw")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Minimum number of clusters:")),
      minClusterSlider, sticky="sw")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Maximum number of clusters:")),
      maxClusterSlider, sticky="sw")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Maximum iterations:")),
      iterNumSlider, sticky="sw")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Number of random seeds:")),
      seedNumSlider, sticky="sw")
    tkgrid(dataFrame, tklabel(top, text="  "), optionsFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=3, sticky="w")
    dialogSuffix(rows=3, columns=3)
    }

# K-Centroids Clustering
kcentroidsClustering <- function(){
	.activeDataSet <- ActiveDataSet()
    initializeDialog(title=gettextRcmdr("K-Centroids Clustering"))
    dataFrame <- tkframe(top)
    xBox <- variableListBox(dataFrame, Numeric(), selectmode="multiple",
      title=gettextRcmdr("Variables (pick one or more)"))
    subsetBoxBCA(dataFrame)
    assignFrame <- tkframe(dataFrame)
    assignClusters <- tclVar("0")
    assignCB <- tkcheckbutton(assignFrame)
    tkconfigure(assignCB, variable=assignClusters)
    assignName <- tclVar("KCentroids")
    assignField <- tkentry(assignFrame, width="15",
      textvariable=assignName)
    radioButtons(name="method", buttons=c("kmn", "kmd", "neuralgas"),
      labels=gettextRcmdr(c("K-Means", "K-Medians", "Neural Gas")),
      title=gettextRcmdr("Clustering Method"))
    optionsFrame <- tkframe(top)
    clusterNumber <- tclVar("2")
    clusterNumSlider <- tkscale(optionsFrame, from=2, to=10, showvalue=TRUE,
      variable=clusterNumber, resolution=1, orient="horizontal")
    seedNumber <- tclVar("10")
    seedNumSlider <- tkscale(optionsFrame, from=1, to=20, showvalue=TRUE,
      variable=seedNumber, resolution=1, orient="horizontal")
    #iterNumber <- tclVar("10")
    #iterNumSlider <- tkscale(optionsFrame, from=5, to=30, showvalue=TRUE,
    #  variable=iterNumber, resolution=5, orient="horizontal")
    summaryClusters <- tclVar("1")
    summaryCB <- tkcheckbutton(optionsFrame)
    tkconfigure(summaryCB, variable=summaryClusters)
    plotCls2d <- tclVar("1")
    plot2dCB <- tkcheckbutton(optionsFrame)
    tkconfigure(plot2dCB, variable=plotCls2d)
    plotCls3d <- tclVar("0")
    plot3dCB <- tkcheckbutton(optionsFrame)
    tkconfigure(plot3dCB, variable=plotCls3d)
    plotFrame <- tkframe(optionsFrame)
    plotPts <- tclVar("1")
    plotPtsCB <- tkcheckbutton(plotFrame)
    tkconfigure(plotPtsCB, variable=plotPts)
    plotCnt <- tclVar("0")
    plotCntCB <- tkcheckbutton(plotFrame)
    tkconfigure(plotCntCB, variable=plotCnt)
    onOK <- function(){
        x <- getSelection(xBox)
        nvar <- length(x)
        subset <- trim.blanks(tclvalue(subsetVariable))
        clusMethod <- tclvalue(methodVariable)
        nClusters <- tclvalue(clusterNumber)
        seeds <- tclvalue(seedNumber)
        #iters <- tclvalue(iterNumber)
        clusterSummary <- tclvalue(summaryClusters)
        clsPlot2d <- tclvalue(plotCls2d)
        clsPlot3d <- tclvalue(plotCls3d)
        ptsPlot <- tclvalue(plotPts)
        cntPlot <- tclvalue(plotCnt)
        clusterAssign <- tclvalue(assignClusters)
        clusterVariable <- trim.blanks(tclvalue(assignName))
        if (clusterAssign == "1" & !is.valid.name(clusterVariable)){
            errorCondition(recall=kcentroidsClustering, message=
              sprintf(gettextRcmdr('"%s" is not a valid name.'),
	      clusterVariable))
	      return()
            }
        closeDialog()
        if (clusterAssign == "1"){
           if (is.element(clusterVariable, Variables())) {
                if ("no" == tclvalue(checkReplace(clusterVariable))){
                    kcentroidsClustering()
                    return()
                    }
                }
           } 
        if (length(x)==0) {
            errorCondition(recall=kcentroidsClustering, 
              message==gettextRcmdr("No variables selected."))
            return()
            }
        varFormula <- paste(x, collapse=" + ")
        vars <- paste(x, collapse=",", sep="")
        .activeDataSet <- ActiveDataSet()
        dset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) .activeDataSet
          else {paste(.activeDataSet, "[", .activeDataSet, "$", subset, ", ]",
            sep="")}
        xmat <- paste("model.matrix(~-1 + ", varFormula, ", ", dset, ")",
          sep="")
        details <- 'FUN = kcca, family = kccaFamily("kmeans")'
        if(clusMethod == "kmd") {
            details <- 'FUN = kcca, family = kccaFamily("kmedians")'
            }
        else if(clusMethod == "neuralgas") {
            details <- 'FUN = cclust, dist = "euclidean", method = "neuralgas"'
            }
        command <- paste("stepFlexclust(", xmat, ", k = ", nClusters,
          ", nrep = ", seeds, ", ", details, ")", sep="")
 	    doItAndPrint(paste(".cluster", "<-", command))
        if (clusterSummary == "1") {
            doItAndPrint("summary(.cluster)")
            doItAndPrint("print(.cluster@centers) # Cluster Centroids")
            }
	if((clsPlot2d=="1" | clsPlot3d=="1") & ptsPlot=="0" & cntPlot=="0") {
	   Message(gettextRcmdr("Bi-plot(s) contains no points or centroids."),
	     type="warning")
	   }
	ptsOpt <- "FALSE"
	if(ptsPlot == "1") ptsOpt <- "TRUE"
	cntOpt <- "FALSE"
	if(cntPlot == "1") cntOpt <- "TRUE"
        if (clsPlot2d == "1") {
            plotCmd2d <- paste("bpCent(prcomp(", xmat, 
              "), clusters(.cluster), data.pts = ", ptsOpt, ", centroids = ",
	          cntOpt, ", xlabs = as.character(clusters(.cluster)))", sep="")
           justDoIt(plotCmd2d)
           logger(plotCmd2d)
           }
        if (clsPlot3d == "1") {
            plotCmd3d <- paste("bpCent3d(prcomp(", xmat, 
              "), clusters(.cluster), data.pts = ", ptsOpt, ", centroids = ",
	          cntOpt, ", xlabs = as.character(clusters(.cluster)))", sep="")
           justDoIt(plotCmd3d)
           logger(plotCmd3d)
           .Tcl("update")
           activateMenus()
           rgl::rgl.bringtotop()
           }
        if (clusterAssign == "1") {
#            assignCommand <- paste(.activeDataSet, "$", clusterVariable,
#              " <- assignCluster(", xmat, ", ", .activeDataSet,
#              ", clusters(.cluster))", sep="")
            assignCommand <- paste(.activeDataSet, "$", clusterVariable,
              " <- as.factor(clusters(.cluster))", sep="")
            justDoIt(assignCommand)
            logger(assignCommand)
            activeDataSet(.activeDataSet)
            }
        justDoIt(paste("remove(.cluster)"))
        logger(paste("remove(.cluster)"))
        activateMenus()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="kcca")
    tkgrid(getFrame(xBox), sticky="nw")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(tklabel(assignFrame, 
      text=gettextRcmdr("Assign clusters to\nthe data set         ")),
      assignCB, sticky="w")
    tkgrid(tklabel(assignFrame, text=gettextRcmdr("Assignment variable: ")),
      assignField, sticky="w")
    tkgrid(assignFrame, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Number of clusters:")),
      clusterNumSlider, sticky="sw")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Number of starting seeds:")),
      seedNumSlider, sticky="sw")
    #tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Maximum iterations:")),
    #  iterNumSlider, sticky="sw")
    tkgrid(tklabel(optionsFrame, 
      text=gettextRcmdr("Print cluster summary")), summaryCB, sticky="w")
    tkgrid(tklabel(optionsFrame, 
      text=gettextRcmdr("2D bi-plot of clusters")), plot2dCB, sticky="w")
    tkgrid(tklabel(optionsFrame, 
      text=gettextRcmdr("3D bi-plot of clusters")), plot3dCB, sticky="w")
    tkgrid(tklabel(plotFrame, text=gettextRcmdr("Plot points")), plotPtsCB,
      tklabel(plotFrame, text=gettextRcmdr("Plot centroids")), plotCntCB,
      sticky="w")
    tkgrid(plotFrame, columnspan=2, sticky="w")
    tkgrid(dataFrame, tklabel(top, text="  "), methodFrame, optionsFrame,
        sticky="nw")
    tkgrid(buttonsFrame, columnspan=3, sticky="w")
    dialogSuffix(rows=3, columns=4)
    }

# The bootstrap replicate diagnostic testing for K-Centroids
bootDiagnostics <- function(){
    initializeDialog(title=gettextRcmdr("K-Centroids Clustering Diagnostics"))
    dataFrame <- tkframe(top)
    xBox <- variableListBox(dataFrame, Numeric(), selectmode="multiple",
      title=gettextRcmdr("Variables (pick one or more)"))
    subsetBoxBCA(dataFrame)
    radioButtons(name="method", buttons=c("kmn", "kmd", "neuralgas"),
      labels=gettextRcmdr(c("K-Means", "K-Medians", "Neural Gas")),
      title=gettextRcmdr("Clustering Method"))
    optionsFrame <- tkframe(top)
    minCluster <- tclVar("2")
    minClusterSlider <- tkscale(optionsFrame, from=2, to=8, showvalue=TRUE,
      variable=minCluster, resolution=1, orient="horizontal")
    maxCluster <- tclVar("4")
    maxClusterSlider <- tkscale(optionsFrame, from=3, to=10, showvalue=TRUE,
      variable=maxCluster, resolution=1, orient="horizontal")
    bootRep <- tclVar("100")
    bootRepSlider <- tkscale(optionsFrame, from=50, to=200, showvalue=TRUE,
      variable=bootRep, resolution=10, orient="horizontal")
    seedNumber <- tclVar("3")
    seedNumSlider <- tkscale(optionsFrame, from=1, to=5, showvalue=TRUE,
      variable=seedNumber, resolution=1, orient="horizontal")
    .activeDataSet <- ActiveDataSet()
     onOK <- function(){
        x <- getSelection(xBox)
        nvar <- length(x)
        subset <- tclvalue(subsetVariable)
        clusMethod <- tclvalue(methodVariable)
        minClus <- tclvalue(minCluster)
        maxClus <- tclvalue(maxCluster)
        bReps <- tclvalue(bootRep)
        seeds <- tclvalue(seedNumber)
        closeDialog()
        if (length(x)==0) {
            errorCondition(recall=bootDiagnostics, 
              message=gettextRcmdr("No variables selected."))
            return()
            }
        if (as.numeric(minClus) > as.numeric(maxClus)) {
            errorCondition(recall=SDIndexPlot, message=gettextRcmdr(
              "The minimum number of clusters selected exceeds the maximum."))
            return()
            }
        varFormula <- paste(x, collapse=" + ")
        vars <- paste(x, collapse=",", sep="")
        dset <- if (subset == gettextRcmdr("<all valid cases>")) .activeDataSet
          else {paste(.activeDataSet, "[", .activeDataSet, "$", subset, ", ]",
            sep="")}
        xmat <- paste("model.matrix(~-1 + ", varFormula, ", ", dset, ")",
          sep="")
        details <- 'FUN = cclust, dist = "euclidean", method = "kmeans"'
        if(clusMethod == "kmd") {
            details <- 'FUN = kcca, family = kccaFamily("kmedians")'
            }
        else if(clusMethod == "neuralgas") {
            details <- 'FUN = cclust, dist = "euclidean", method = "neuralgas"'
            }
        command <- paste("bootCVD(", xmat, ", k = ", minClus, ":", maxClus,
            ", nboot = ", bReps, ", nrep = ", seeds, ', method = "', clusMethod,
            '", col1 = 1, col2 = 2, dsname = "', .activeDataSet, '")', sep="")
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="bootCVD")
    tkgrid(getFrame(xBox), sticky="nw")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Minimum number of clusters:")),
      minClusterSlider, sticky="sw")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Maximum number of clusters:")),
      maxClusterSlider, sticky="sw")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Bootstrap replicates:")),
      bootRepSlider, sticky="sw")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Number of random seeds:")),
      seedNumSlider, sticky="sw")
    tkgrid(dataFrame, tklabel(top, text="  "), methodFrame, optionsFrame, sticky="nw")
    tkgrid(buttonsFrame, columnspan=3, sticky="w")
    dialogSuffix(rows=3, columns=4)
    }

# List hierachical clustering solutions
listHclustSolutions <- function(envir=.GlobalEnv, ...) {
    objects <- ls(envir=envir, ...)
    if (length(objects) == 0) NULL
    else objects[sapply(objects,
        function(.x) "hclust" == class(get(.x, envir=envir))[1]) ]
#        function(.x) "hclust" == (class(eval(parse(text=.x), envir=envir))[1]))]
    }

# Hierarchical clustering dialog box with the subsetBox correction
hCluster <- function(){
    solutionNumber=length(listHclustSolutions())
    initializeDialog(title=gettextRcmdr("Hierarchical Clustering"))
    solutionFrame <- tkframe(top)
    solutionName <- tclVar(paste("HClust.", (solutionNumber+1),
        sep=""))
    solutionField <- ttkentry(solutionFrame, width="20",
      textvariable=solutionName)
    dataFrame <- tkframe(top)
    xBox <- variableListBox(dataFrame, Numeric(), selectmode="multiple",
      title=gettextRcmdr("Variables (pick one or more)"))
    subsetBoxBCA(dataFrame)
    radioButtons(name="method",
      buttons=c("ward", "single", "complete","average", "mcquitty", "median",
      "centroid"), labels=gettextRcmdr(c("Ward's Method", "Single Linkage",
      "Complete Linkage", "Average Linkage", "McQuitty's Method",
      "Median Linkage", "Centroid Linkage")), title=gettextRcmdr("Clustering Method"))
    optionsFrame <- tkframe(top)
    radioButtons(optionsFrame, name="distanceType", buttons=c("euc", "euc2",
      "city", "none"), labels=gettextRcmdr(c("Euclidean", "Squared-Euclidian", 
      "Manhattan (City Block)", "No Transformation")), initialValue="euc2", 
      title=gettextRcmdr("Distance Measure"))
    checkFrame <- tkframe(optionsFrame)
    plotDendro <- tclVar("1")
    plotCB <- tkcheckbutton(checkFrame)
    tkconfigure(plotCB, variable=plotDendro)
    onOK <- function(){
        x <- getSelection(xBox)
        nvar <- length(x)
        clusMethod <- tclvalue(methodVariable)
        distance <- tclvalue(distanceTypeVariable)
        subset <- trim.blanks(tclvalue(subsetVariable))
        dendro <- tclvalue(plotDendro)
        solution <- trim.blanks(tclvalue(solutionName))
        if (length(x)==0) {
            errorCondition(recall=hierarchicalCluster, 
              message=gettextRcmdr("No variables selected."))
            return()
            }
        closeDialog()
        varFormula <- paste(x, collapse="+")
        vars <- paste(x, collapse=",", sep="")
        .activeDataSet <- ActiveDataSet()
        dset <- if (subset == gettextRcmdr("<all valid cases>")) .activeDataSet
          else {paste(.activeDataSet, "[", .activeDataSet, "$", subset, ", ]",
            sep="")}
        xmat <- paste("model.matrix(~-1 + ", varFormula, ", ", dset, ")",
          sep="")
        if(distance=="euc") {
            dx <- paste("dist(", xmat, ")", sep="")
            distlab <- "euclidian"
        }
        else if(distance=="euc2") {
            dx <- paste("dist(", xmat, ")^2", sep="")
            distlab <- "squared-euclidian"
        }
        else if(distance=="city") {
            dx <- paste("dist(", xmat, ", method= ", '"manhattan"', ")",
                sep="")
            distlab <- "city-block"
        }
        else {
            dx <- xmat
            distlab <- "untransformed"
        }
        command <- paste("hclust(", dx, " , method= ", '"', clusMethod, '"',
          ")", sep="")
  	    doItAndPrint(paste(solution, "<-", command))
       if (dendro == "1") {
            justDoIt(paste("plot(", solution, ", main= ",'"',
              "Cluster Dendrogram for Solution ", solution, '"', ", xlab= ",
              '"',"Observation Number in Data Set ", dset, '"',
               ", sub=", '"', "Method=", clusMethod,
              "; Distance=", distlab, '"', ")", sep=""))
            logger(paste("plot(", solution, ", main= ",'"',
              "Cluster Dendrogram for Solution ", solution, '"', ", xlab= ",
              '"',"Observation Number in Data Set ", dset, '"',
               ", sub=", '"', "Method=", clusMethod,
              "; Distance=", distlab, '"', ")",
              sep=""))
            }
         activateMenus()
         tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="hclust")
    tkgrid(solutionField, sticky="w")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Clustering solution name:")),
      solutionFrame, sticky="w")
    tkgrid(getFrame(xBox), sticky="nw")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(distanceTypeFrame, sticky="w")
    tkgrid(labelRcmdr(checkFrame, text="  "), sticky="w")
    tkgrid(labelRcmdr(checkFrame, text=gettextRcmdr("Plot Dendrogram  ")), plotCB,
      sticky="w")
    tkgrid(checkFrame, sticky="w")
    tkgrid(dataFrame, methodFrame, optionsFrame, sticky="nw")
    tkgrid(buttonsFrame, columnspan=3, sticky="w")
    dialogSuffix(rows=3, columns=3)
    }

hclustSummaryBCA <- function(){
    parseDataSet <- function(x) {
        y <- eval(parse(text=paste(x, "$call", sep="")))
        string1 <- unlist(strsplit(as.character(y)[2], "\\("))
        string2 <- unlist(strsplit(string1[3], ","))
        if(length(grep("\\[", string2[2])) == 0) {
            out <- gsub(")", "", gsub(" ", "", gsub("\\^2", "",
                string2[2])))
            }
        else {
            string3 <- unlist(strsplit(string2[2], "\\["))
            out <- gsub(" ", "", string3[1])
            }
        return(out)
        }
    hclustObjects <- listHclustSolutions()
    testDataSet <- tapply(hclustObjects, as.factor(1:length(hclustObjects)),
      parseDataSet)
    .activeDataSet <- ActiveDataSet()
    validHclust <- hclustObjects[testDataSet==.activeDataSet]
    initializeDialog(
      title=gettextRcmdr("Hierarchical Cluster Summary"))
    hclustBox <- variableListBox(top, validHclust, selectmode="single",
      title=gettextRcmdr("Select One Clustering Solution"))
    optionsFrame <- tkframe(top)
    clusterNumber <- tclVar("2")
    slider <- tkscale(optionsFrame, from=2, to=10, showvalue=TRUE,
      variable=clusterNumber, resolution=1, orient="horizontal")
    summaryClusters <- tclVar("1")
    summaryCB <- tkcheckbutton(optionsFrame)
    tkconfigure(summaryCB, variable=summaryClusters)
    plotCls2d <- tclVar("1")
    plot2dCB <- tkcheckbutton(optionsFrame)
    tkconfigure(plot2dCB, variable=plotCls2d)
    plotCls3d <- tclVar("0")
    plot3dCB <- tkcheckbutton(optionsFrame)
    tkconfigure(plot3dCB, variable=plotCls3d)
    plotFrame <- tkframe(optionsFrame)
    plotPts <- tclVar("1")
    plotPtsCB <- tkcheckbutton(plotFrame)
    tkconfigure(plotPtsCB, variable=plotPts)
    plotCnt <- tclVar("0")



    plotCntCB <- tkcheckbutton(plotFrame)
    tkconfigure(plotCntCB, variable=plotCnt)
    if(length(hclustObjects)==0) {
        errorCondition(recall=return,
          message=gettextRcmdr("There are no hierachical clustering solutions"))
        }
    if(length(validHclust)==0) {
        errorCondition(recall=return, message=
     gettextRcmdr("No hierachical clustering solutions are associated with this data set."))
        }
   onOK <- function(){
        solution <- getSelection(hclustBox)
        if(length(solution)==0) {
          errorCondition(recall=hclustSummary,
            message=gettextRcmdr("A clustering solution has not been selected."))
          return()
            }
        clusters <- as.numeric(tclvalue(clusterNumber))
        clusterVar <- paste("cutree(", solution, ", k = ", clusters, ")",
          sep="")
        clusterSummary <- tclvalue(summaryClusters)
        clsPlot2d <- tclvalue(plotCls2d)
        clsPlot3d <- tclvalue(plotCls3d)
        ptsPlot <- tclvalue(plotPts)
        cntPlot <- tclvalue(plotCnt)
        hclustCall <- eval(parse(text=paste(solution,"$call",sep="")))
        string1 <- unlist(strsplit(as.character(hclustCall)[2], "\\("))
        string2 <- unlist(strsplit(string1[3], ","))
        form.vars <- string2[1]
        closeDialog()
        if(length(grep("\\[", string2[2])) == 0) {
            xmat <- paste("model.matrix(", form.vars, ", ", .activeDataSet, ")",
              sep="")
            }
        else {
            string3 <- unlist(strsplit(string2[2], "\\["))
            xmat <- paste("model.matrix(", form.vars, ", ", .activeDataSet, "[",
              string3[2], ", ]",")", sep="")
            }
        if (clusterSummary == "1") {
            doItAndPrint(paste("summary(as.factor(", clusterVar,
              ")) # Cluster Sizes", sep=""))
			#centroidsCommand <- paste("apply(", xmat, ", 2, function(x) tapply(x, as.factor(", clusterVar, "), mean)) # Cluster Centroids", sep = "")
            centroidsCommand <- paste("by(", xmat, ", as.factor(", clusterVar,
              "), colMeans) # Cluster Centroids", sep="")
            doItAndPrint(centroidsCommand)
            }
	if((clsPlot2d=="1" | clsPlot3d=="1") & ptsPlot=="0" & cntPlot=="0") {
	   Message(gettextRcmdr("Bi-plot(s) contains no points or centroids."),
	     type="warning")
	   }
	ptsOpt <- "FALSE"
	if(ptsPlot == "1") ptsOpt <- "TRUE"
	cntOpt <- "FALSE"
	if(cntPlot == "1") cntOpt <- "TRUE"
        if (clsPlot2d == "1") {
            plotCmd2d <- paste("bpCent(prcomp(", xmat, 
              "), ", clusterVar, ", data.pts = ", ptsOpt, ", centroids = ",
	      cntOpt, ", xlabs = as.character(", clusterVar, "))", sep="")
           justDoIt(plotCmd2d)
           logger(plotCmd2d)
           }
        if (clsPlot3d == "1") {
            plotCmd3d <- paste("bpCent3d(prcomp(", xmat, 
              "), ", clusterVar, ", data.pts = ", ptsOpt, ", centroids = ",
	      cntOpt, ", xlabs = as.character(", clusterVar, "))", sep="")
           justDoIt(plotCmd3d)
           logger(plotCmd3d)
          .Tcl("update")
           activateMenus()
           rgl::rgl.bringtotop()
           }
        activateMenus()
        tkfocus(CommanderWindow())
        } 
    OKCancelHelp(helpSubject="biplot")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Number of clusters:")), slider,
      sticky="sw")
    tkgrid(tklabel(optionsFrame, 
      text=gettextRcmdr("Print cluster summary")), summaryCB, sticky="w")
    tkgrid(tklabel(optionsFrame, 
      text=gettextRcmdr("2D bi-plot of clusters")), plot2dCB, sticky="w")
    tkgrid(tklabel(optionsFrame, 
      text=gettextRcmdr("3D bi-plot of clusters")), plot3dCB, sticky="w")
    tkgrid(tklabel(plotFrame, text=gettextRcmdr("Plot points")), plotPtsCB,
      tklabel(plotFrame, text=gettextRcmdr("Plot centroids")), plotCntCB,
      sticky="w")
    tkgrid(plotFrame, columnspan=2, sticky="w")
    tkgrid(getFrame(hclustBox), optionsFrame, sticky="nw")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=2, columns=3)
    }

# STATISTICS-MODELS MENU
# GLM changes to add the McFadden R^2 value to the model summary
generalizedLinearModelBCA <- function(){
    families <- c("gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian",
        "quasibinomial", "quasipoisson")
    links <- c("identity", "inverse", "log", "logit", "probit",
        "cloglog", "sqrt", "1/mu^2")
    availableLinks <- matrix(c(
        TRUE,  TRUE,  TRUE,  FALSE, FALSE, FALSE, FALSE, FALSE,
        FALSE, FALSE, FALSE, TRUE,  TRUE,  TRUE,  FALSE, FALSE,
        TRUE,  FALSE, TRUE,  FALSE, FALSE, FALSE, TRUE,  FALSE,
        TRUE,  TRUE,  TRUE,  FALSE, FALSE, FALSE, FALSE, FALSE,
        TRUE,  TRUE,  TRUE,  FALSE, FALSE, FALSE, FALSE, TRUE,
        FALSE, FALSE, FALSE, TRUE,  TRUE,  TRUE,  FALSE, FALSE,
        TRUE,  FALSE, TRUE,  FALSE, FALSE, FALSE, TRUE,  FALSE),
        7, 8, byrow=TRUE)
    rownames(availableLinks) <- families
    colnames(availableLinks) <- links
    canonicalLinks <- c("identity", "logit", "log", "inverse", "1/mu^2", "logit", "log")
    names(canonicalLinks) <- families
    initializeDialog(title=gettextRcmdr("Generalized Linear Model"))
    .activeModel <- ActiveModel()
    currentModel <- if (!is.null(.activeModel))
        class(get(.activeModel, envir=.GlobalEnv))[1] == "glm"
#        eval(parse(text=paste("class(", .activeModel, ")[1] == 'glm'", sep="")),
#            envir=.GlobalEnv)
        else FALSE
    if (currentModel) {
        currentFields <- formulaFields(get(.activeModel, envir=.GlobalEnv), glm=TRUE)
#        currentFields <- formulaFields(eval(parse(text=.activeModel),
#            envir=.GlobalEnv), glm=TRUE)
        if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
        }
    modelFormula()
    UpdateModelNumber()
    modelName <- tclVar(paste("GLM.", getRcmdr("modelNumber"), sep=""))
    modelFrame <- tkframe(top)
    model <- ttkentry(modelFrame, width="20", textvariable=modelName)
    linkFamilyFrame <- tkframe(top)
    familyFrame <- tkframe(linkFamilyFrame)
    familyBox <- tklistbox(familyFrame, height="4", exportselection="FALSE",
        selectmode="single", background="white")
    familyScroll <- ttkscrollbar(familyFrame,
        command=function(...) tkyview(familyBox, ...))
    tkconfigure(familyBox, yscrollcommand=function(...) tkset(familyScroll, ...))
    for (fam in families) tkinsert(familyBox, "end", fam)
    linkFrame <- tkframe(linkFamilyFrame)
    linkBox <- tklistbox(linkFrame, height="4", exportselection="FALSE",
        selectmode="single", background="white")
    subsetBox(model=TRUE)
    onFamilySelect <- function(){
        family <- families[as.numeric(tkcurselection(familyBox)) + 1]
        availLinks <- links[availableLinks[family,]]
        tkdelete(linkBox, "0", "end")
        for (lnk in availLinks) tkinsert(linkBox, "end", lnk)
        canLink <- canonicalLinks[family]
        tkconfigure(linkBox, height=length(availLinks))
        tkselection.set(linkBox, which(canLink == availLinks) - 1)
        }
    onOK <- function(){
        check.empty <- gsub(" ", "", tclvalue(lhsVariable))
        if ("" == check.empty) {
            errorCondition(recall=generalizedLinearModel, model=TRUE, message=gettextRcmdr("Left-hand side of model empty."))
            return()
            }
        check.empty <- gsub(" ", "", tclvalue(rhsVariable))
        if ("" == check.empty) {
            errorCondition(recall=generalizedLinearModel, model=TRUE, message=gettextRcmdr("Right-hand side of model empty."))
            return()
            }
        modelValue <- trim.blanks(tclvalue(modelName))
        if (!is.valid.name(modelValue)){
            errorCondition(recall=generalizedLinearModel, model=TRUE, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue))
            return()
            }
        if (is.element(modelValue, listGeneralizedLinearModels())) {
            if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
                UpdateModelNumber(-1)
                closeDialog()
                generalizedLinearModel()
                return()
                }
            }
        formula <- paste(tclvalue(lhsVariable), tclvalue(rhsVariable), sep=" ~ ")
        family <- families[as.numeric(tkcurselection(familyBox)) + 1]
        availLinks <- links[availableLinks[family,]]
        link <- availLinks[as.numeric(tkcurselection(linkBox)) + 1]
        subset <- tclvalue(subsetVariable)
        closeDialog()
        # The five lines below added to address the "=" vs "==" problem
        if(is.equality.prob(trim.blanks(subset))) {
            errorCondition(recall=GenerializedLinearModel, 
            message=gettextRcmdr('Use "==" not "=" for equalties in subsetting.'))
            return()
            }
        if (trim.blanks(subset) == gettextRcmdr("<all valid cases>") || trim.blanks(subset) == ""){
            subset <- ""
            putRcmdr("modelWithSubset", FALSE)
            }
        else{
            subset <- paste(", subset=", subset, sep="")
            putRcmdr("modelWithSubset", TRUE)
            }
        command <- paste("glm(", formula, ", family=", family, "(", link,
            "), data=", ActiveDataSet(), subset, ")", sep="")
 	    doItAndPrint(paste(modelValue, "<-", command))
        doItAndPrint(paste("summary(", modelValue, ")", sep=""))
        # The McFadden R^2 computation added by Dan Putler 29Sep06
        doItAndPrint(paste("1 - (", modelValue, "$deviance/",
          modelValue, "$null.deviance) # McFadden R2", sep=""))
        activeModel(modelValue)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="generalizedLinearModel")
    helpButton <- buttonRcmdr(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(getFrame(xBox), sticky="w")
    tkgrid(outerOperatorsFrame, sticky="w")
    tkgrid(formulaFrame, sticky="w")
    tkgrid(subsetFrame, sticky="w")
        tkgrid(labelRcmdr(linkFamilyFrame, text=gettextRcmdr("Family (double-click to select)"), fg="blue"),
        labelRcmdr(linkFamilyFrame, text="   "), labelRcmdr(linkFamilyFrame, text=gettextRcmdr("Link function"), fg="blue"), sticky="w")
    tkgrid(familyBox, familyScroll, sticky="nw")
    tkgrid(linkBox, sticky="nw")
    tkgrid(familyFrame, labelRcmdr(linkFamilyFrame, text="   "), linkFrame, sticky="nw")
    tkgrid(linkFamilyFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(familyScroll, sticky="ns")
    fam <- if (currentModel) which(currentFields$family == families) - 1
        else 1
    tkselection.set(familyBox, fam)
    availLinks <- links[availableLinks[fam + 1,]]
    for (lnk in availLinks) tkinsert(linkBox, "end", lnk)
    tkconfigure(linkBox, height=length(availLinks))
    lnk <- if (currentModel) which(currentFields$link == availLinks) - 1
            else 0
    tkselection.set(linkBox, lnk)
    tkbind(familyBox, "<Double-ButtonPress-1>", onFamilySelect)
    dialogSuffix(rows=7, columns=1, focus=lhsEntry, preventDoubleClick=TRUE)
    }

# Helper functions for the added nnet and rpart models
listNnetModels <- function(envir=.GlobalEnv, ...) {
    objects <- ls(envir=envir, ...)
    if (length(objects) == 0) NULL
    else objects[sapply(objects, 
        function(.x) "nnet.formula" == 
          (class(eval(parse(text=.x), envir=envir))[1]))]
    }

nnetP <- function() length(listNnetModels()) > 0

listRpartModels <- function(envir=.GlobalEnv, ...) {
    objects <- ls(envir=envir, ...)
    if (length(objects) == 0) NULL
    else objects[sapply(objects, 
        function(.x) "rpart" == (class(eval(parse(text=.x), envir=envir))[1]))]
    }

rpartP <- function() length(listRpartModels()) > 0

activeRpartP <- function() activeModelP() && eval(parse(text=paste("class(",
  ActiveModel(), ")[1] == 'rpart'")))

# The neural network functions
nnetModel <- function(){
    initializeDialog(title=gettextRcmdr("Neural Network Model"))
    .activeModel <- ActiveModel()
    .activeDataSet <- ActiveDataSet()
    currentModel <- if (!is.null(.activeModel))
        class(get(.activeModel, envir=.GlobalEnv))[1] == "nnet.formula"
#        eval(parse(text=paste("class(", .activeModel, ")[1] == 'nnet.formula'",
#          sep="")), envir=.GlobalEnv)
        else FALSE
    if (currentModel) {
        currentFields <- formulaFields(get(.activeModel, envir=.GlobalEnv))
#        currentFields <- formulaFields(eval(parse(text=.activeModel),
#            envir=.GlobalEnv))
        if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
        }
    UpdateModelNumber()
    modelName <- tclVar(paste("NNET.", getRcmdr("modelNumber"), sep=""))
    modelFrame <- tkframe(top)
    model <- ttkentry(modelFrame, width="20", textvariable=modelName)
#    model <- tkentry(modelFrame, width="20", textvariable=modelName)
#    subsetBox(model=TRUE)
    subsetBoxBCA(model=TRUE)
    optionsFrame <- tkframe(top)
    sizeNumber <- tclVar("1")
    sizeNumSlider <- tkscale(optionsFrame, from=1, to=10, showvalue=TRUE,
      variable=sizeNumber, resolution=1, orient="horizontal")
    decayValue <- tclVar("0")
    decayValSlider <- tkscale(optionsFrame, from=0, to=1, showvalue=TRUE,
      variable=decayValue, resolution=0.05, orient="horizontal")
    onOK <- function(){
        check.empty <- gsub(" ", "", tclvalue(lhsVariable))
        modelValue <- trim.blanks(tclvalue(modelName))
        closeDialog()
        if (!is.valid.name(modelValue)){
            errorCondition(recall=nnetModel, message=
	      sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue),
	      model=TRUE)
            return()
            }
        subset <- tclvalue(subsetVariable)
        if(is.equality.prob(trim.blanks(subset))) {
            errorCondition(recall=nnetModel, 
            message=gettextRcmdr('Use "==" not "=" for equalties in subsetting.'))
            return()
            }
        if (trim.blanks(subset) == gettextRcmdr("<all valid cases>") || trim.blanks(subset) == ""){
            subset <- ""
            putRcmdr("modelWithSubset", FALSE)
            }
        else{
            putRcmdr(".subset", subset)
            subset <- paste(", subset='",subset, "'", sep="")
#            subset <- paste(", subset=", subset, sep="")
            putRcmdr("modelWithSubset", TRUE)
            }
        check.empty <- gsub(" ", "", tclvalue(lhsVariable))
        if ("" == check.empty) {
            errorCondition(recall=nnetModel, message=gettextRcmdr("Left-hand side of model empty."), model=TRUE)
            return()
            }
        check.empty <- gsub(" ", "", tclvalue(rhsVariable))
        if ("" == check.empty) {
            errorCondition(recall=nnetModel, message=gettextRcmdr("Right-hand side of model empty."), model=TRUE)
            return()
            }
        if (!is.factor(eval(parse(text=tclvalue(lhsVariable)), envir=eval(parse(text=.activeDataSet), envir=.GlobalEnv)))){
            errorCondition(recall=nnetModel, message=gettextRcmdr("Response variable must be a factor"))
            return()
            }
        if (is.element(modelValue, listNnetModels())) {
            if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
                UpdateModelNumber(-1)
                nnetModel()
                return()
                }
            }
        formula <- paste(trim.blanks(tclvalue(lhsVariable)),
          trim.blanks(tclvalue(rhsVariable)), sep=" ~ ")
        sizeNum <- tclvalue(sizeNumber)
        decayVal <- tclvalue(decayValue)
#        command <- paste("Nnet(", formula, ", data=", .activeDataSet, 
#        command <- paste("nnet(", formula, ", data=", .activeDataSet, 
#          ", decay=", decayVal, ", size=", sizeNum, subset, ")", sep="")
        if(subset=="") {
            command <- paste("Nnet(", formula, ", data=", .activeDataSet, 
              ", decay=", decayVal, ", size=", sizeNum, ")", sep="")
            }
        else {
            command <- paste("Nnet(", formula, ", data=", .activeDataSet, 
              ", decay=", decayVal, ", size=", sizeNum, subset, ")", sep="")
            }
 	    doItAndPrint(paste(modelValue, "<-", command))
        doItAndPrint(paste(modelValue,"$value # Final Objective Function Value",
          sep=""))
        doItAndPrint(paste("summary(", modelValue, ")", sep=""))
        activeModel(modelValue)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="Nnet")
# This is new and seems to control the styling of the help button
    helpButton <- buttonRcmdr(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    modelFormula()
    subsetBox(model=TRUE)
    tkgrid(getFrame(xBox), sticky="w")

    tkgrid(formulaFrame, sticky="w")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Hidden layer units:")),
        sizeNumSlider, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Decay weight parameter:")),
        decayValSlider, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1, focus=lhsEntry, preventDoubleClick=TRUE)
    }
# These are no longer in use    
#Printcp <- function(rpartObj) {
#    printcp(rpartObj)
#    invisible()
#    }

#Print <- function(rpartObj) {
#    print(rpartObj)
#    invisible()
#    }

# The rpart functions
rpartModel <- function(){
    initializeDialog(title=gettextRcmdr("rpart Tree"))
    .activeModel <- ActiveModel()
    .activeDataSet <- ActiveDataSet()
    currentModel <- if (!is.null(.activeModel))
      class(get(.activeModel, envir=.GlobalEnv))[1] == "rpart"
#      eval(parse(text=paste("class(", .activeModel, ")[1] == 'rpart'",
#        sep="")), envir=.GlobalEnv)
      else FALSE
    if (currentModel) {
        currentFields <- formulaFields(get(.activeModel, envir=.GlobalEnv))
#        currentFields <- formulaFields(eval(parse(text=.activeModel),
#            envir=.GlobalEnv))
        if (currentFields$data != .activeDataSet) currentModel <- FALSE
        }
    UpdateModelNumber()
    modelName <- tclVar(paste("RPART.", getRcmdr("modelNumber"), sep=""))
    modelFrame <- tkframe(top)
#    model <- tkentry(modelFrame, width="20", textvariable=modelName)
    model <- ttkentry(modelFrame, width="20", textvariable=modelName)
#    subsetBox(model=TRUE)
    subsetBoxBCA(model=TRUE)
    optionsFrame <- tkframe(top)
    complexParam <- tclVar("0.01")
    comParField <- tkentry(optionsFrame, width="6", textvariable=complexParam)
    printTree <- tclVar("1")
    printTreeCB <- tkcheckbutton(optionsFrame)
    tkconfigure(printTreeCB, variable=printTree)
    printPrune <- tclVar("1")
    printCB <- tkcheckbutton(optionsFrame)
    tkconfigure(printCB, variable=printPrune)
    plotPrune <- tclVar("1")
    plotCB <- tkcheckbutton(optionsFrame)
    tkconfigure(plotCB, variable=plotPrune)
    onOK <- function(){
        check.empty <- gsub(" ", "", tclvalue(lhsVariable))
        closeDialog()
        if ("" == check.empty) {
            errorCondition(recall=rpartModel, model=TRUE,
              message=gettextRcmdr("Left-hand side of model empty."))
            return()
            }
        check.empty <- gsub(" ", "", tclvalue(rhsVariable))
        if ("" == check.empty) {
            errorCondition(recall=rpartModel, model=TRUE,
              message=gettextRcmdr("Right-hand side of model empty."))
            return()
            }
        modelValue <- trim.blanks(tclvalue(modelName))
        if (!is.valid.name(modelValue)){
            errorCondition(recall=rpartModel, model=TRUE,
              message=sprintf(gettextRcmdr('"%s" is not a valid name.'),
	      modelValue))
            return()
            }
        subset <- tclvalue(subsetVariable)
        if(is.equality.prob(trim.blanks(subset))) {
            errorCondition(recall=rpartModel, 
            message=gettextRcmdr('Use "==" not "=" for equalties in subsetting.'))
            return()
            }
        if (is.element(modelValue, listRpartModels())) {
            if ("no" == tclvalue(checkReplace(modelValue,
              type=gettextRcmdr("Model")))){
                UpdateModelNumber(-1)
                rpartModel()
                return()
                }
            }
        formula <- paste(trim.blanks(tclvalue(lhsVariable)),
          trim.blanks(tclvalue(rhsVariable)), sep=" ~ ")
        if (trim.blanks(subset) == 
          gettextRcmdr("<all valid cases>") || trim.blanks(subset) == ""){
            subset <- ""
            putRcmdr("modelWithSubset", FALSE)
            }
        else{
            putRcmdr(".subset", subset)
            subset <- paste(", subset=", subset, sep="")
            putRcmdr("modelWithSubset", TRUE)
            }
        complexPar <- tclvalue(complexParam)
        treePrint <- tclvalue(printTree)
        prunePrint <- tclvalue(printPrune)
        prunePlot <- tclvalue(plotPrune)
        command <- paste("rpart(", formula, ", data=", .activeDataSet, 
          ", cp=" , complexPar, subset,")", sep="")
#        if(subset=="") {
#            command <- paste("rpart(", formula, ", data=", .activeDataSet, 
#              ", cp=" , complexPar, ")", sep="")
#            }
#        else {
#            command <- paste("rpart(", formula, ", data=", .activeDataSet, 
#              subset, ", cp=" , complexPar, ")", sep="")
#            }
 	    doItAndPrint(paste(modelValue, "<-", command))
        if(prunePlot == "1") {
            logger(paste("plotcp(", modelValue, ")", sep=""))
            justDoIt(paste("plotcp(", modelValue, ")", sep=""))
            }
        if(prunePrint == "1") {
#            doItAndPrint(paste("Printcp(", modelValue, ") # Pruning Table",
#              sep=""))
            doItAndPrint(paste("printcp(", modelValue, ") # Pruning Table",
              sep=""))
        }
        if(treePrint == "1") {
#            doItAndPrint(paste("Print(", modelValue, ") # Tree Leaf Summary",
#              sep=""))
#            doItAndPrint(paste("print(", modelValue, ") # Tree Leaf Summary",
#              sep=""))
            doItAndPrint(paste(modelValue, " # Tree Leaf Summary", sep=""))
        }
        activeModel(modelValue)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="rpart")
    helpButton <- buttonRcmdr(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(modelFrame, text=gettextRcmdr("Enter name for model:")),
      model, sticky="w")
    tkgrid(modelFrame, sticky="w")    
    modelFormula()
    tkgrid(getFrame(xBox), sticky="w")
    tkgrid(formulaFrame, sticky="w")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(tklabel(optionsFrame, text=" "), sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Complexity parameter:")),
      comParField, sticky="w")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Print pruning table:")),
      printCB, tklabel(optionsFrame, text=gettextRcmdr("Plot pruning table:")),
      plotCB,tklabel(optionsFrame, text=gettextRcmdr(
      "Print tree leaf summary:")), printTreeCB, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1, focus=lhsEntry, preventDoubleClick=TRUE)
    }

rpartPlot <- function(){
    rpartMod <- ActiveModel()
    initializeDialog(title=gettextRcmdr("Plot an rpart Tree"))
    leafFrame <- tkframe(top)
    radioButtons(name="leaf", buttons=c("counts", "proportions"), 
      labels=gettextRcmdr(c("Counts", "Proportions")), title=gettextRcmdr("Leaf Summary"))
    optionsFrame <- tkframe(top)
    uniform <- tclVar("1")
    uniformCB <- tkcheckbutton(optionsFrame)
    tkconfigure(uniformCB, variable=uniform)
    onOK <- function(){
        leafSum <- tclvalue(leafVariable)
        uniformOp <- tclvalue(uniform)
        closeDialog()
        if(uniformOp == 1) uniformStr <- ", uniform = TRUE, fallen.leaves = FALSE)"
        else uniformStr <- ", uniform = FALSE, fallen.leaves = TRUE)"
        if(leafSum == "counts") leafStr <- "extra = 2"
        else leafStr <- "extra = 4"
        plotCom <- paste("rpart.plot(", rpartMod, ", type = 0, ", leafStr, uniformStr, sep="")
        logger(plotCom)
        justDoIt(plotCom)
        activateMenus()
        tkfocus(CommanderWindow())
        } 
    OKCancelHelp(helpSubject="rpart.plot")
    tkgrid(leafFrame, sticky="w")
    tkgrid(tklabel(optionsFrame, 
      text=gettextRcmdr("Uniform branch distances")), uniformCB, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=4, columns=2)
    }

# MODELS MENU
summarizeModelBCA <- function(){
    .activeModel <- ActiveModel()
    if (!checkMethod("summary", .activeModel)) return()
    modClass <- eval(parse(text=paste("class(",.activeModel,")[1]", sep="")))
    if(modClass == "rpart") {
        doItAndPrint(paste("summary(", .activeModel, ")", sep=""))
    }
    else doItAndPrint(paste("summary(", .activeModel, ", cor=FALSE)", sep=""))
    if(modClass == "glm") {
        doItAndPrint(paste("1 - (", .activeModel, "$deviance/",
          .activeModel, "$null.deviance) # McFadden R^2", sep=""))
        }
    }

stepwiseBCA <- function(){
    initializeDialog(title=gettextRcmdr("Stepwise Variable Selection"))
    modelFrame <- tkframe(top)
    modelType <- eval(parse(text=paste("class(", ActiveModel(), ")[1]", sep="")))
    if(modelType=="lm") {
        modelName <- tclVar(paste("LinearModel.", getRcmdr("modelNumber"), sep=""))
        }
    else if(modelType=="glm") {
        modelName <- tclVar(paste("GLM.", getRcmdr("modelNumber"), sep=""))
        }
    else {
        modelName <- tclVar(paste("Step.", getRcmdr("modelNumber"), sep=""))
        }
    modelField <- tkentry(modelFrame, width="20", textvariable=modelName)
    directionFrame <- tkframe(top)
    radioButtons(name="direction",
      buttons=c("both", "backward", "forward"), labels=gettextRcmdr(c("Both",
      "Backward", "Foward")), title=gettextRcmdr("Search Direction"))
    optionsFrame <- tkframe(top)
    radioButtons(optionsFrame, name="criteria", buttons=c("aic", "bic"), 
      labels=gettextRcmdr(c("AIC", "BIC")), title=gettextRcmdr("Selection Criteria"))
    checkFrame <- tkframe(optionsFrame)
    printSum <- tclVar("1")
    printCB <- tkcheckbutton(checkFrame)
    tkconfigure(printCB, variable=printSum)
    onOK <- function(){
        searchDir <- tclvalue(directionVariable)
        selectCrit <- tclvalue(criteriaVariable)
        sumPrint <- tclvalue(printSum)
        modelValue <- trim.blanks(tclvalue(modelName))
        if (!is.valid.name(modelValue)){
            UpdateModelNumber(-1)
            errorCondition(recall=stepwise, message=
	      sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue),
	      model=TRUE)
            return()
            }
        closeDialog()
        if(selectCrit=="aic") {
            kVal <- 2
            }
        else {
            nob <- eval(parse(text=paste("length(", ActiveModel(), "$residuals)", sep="")))
            kVal <- paste("log(", nob, ")", sep="")
            }
        command <- paste("step(", ActiveModel(), ", direction=", '"', searchDir, '"',
          ", k=", kVal, ")", sep="")
 	    doItAndPrint(paste(modelValue, "<-", command))
        if(sumPrint=="1") {
            doItAndPrint(paste("summary(", modelValue, ")", sep=""))
            if(modelType == "glm") {
                doItAndPrint(paste("1 - (", modelValue, "$deviance/",
                modelValue, "$null.deviance) # McFadden R^2", sep=""))
                }
            }
        activeModel(modelValue)
        activateMenus()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="step")
    tkgrid(modelField, sticky="w")
    tkgrid(tklabel(top, text=gettextRcmdr("Enter the name of the new model:")),
      modelFrame, sticky="w")
    tkgrid(criteriaFrame, sticky="w")
    tkgrid(tklabel(checkFrame, text=gettextRcmdr("Print Summary  ")), printCB,
      sticky="w")
    tkgrid(checkFrame, sticky="w")
    tkgrid(directionFrame, optionsFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=3, columns=2)
    }

## Lift charts and scoring
liftChart <- function() {
	dataSets <- listDataSets()
    parseDataSetGlm <- function(x) {
        y <- eval(parse(text=paste(x, "$call", sep="")))
        out <- unlist(strsplit(as.character(y), "="))[4]
        return(out)
        }
    parseDataSetOther <- function(x) {
        y <- eval(parse(text=paste(x, "$call", sep="")))
        out <- unlist(strsplit(as.character(y), "="))[3]
        return(out)
        }
    .activeDataSet <- ActiveDataSet()
    if(length(listGeneralizedLinearModels()>0)) {
        glmObjects <- listGeneralizedLinearModels()
        testDataSetGlm <- tapply(glmObjects, as.factor(1:length(glmObjects)),
          parseDataSetGlm)
        validGlmObjects <- glmObjects[testDataSetGlm==.activeDataSet]
        }
    else {validGlmObjects <- character(0)}
    if(length(listRpartModels()) > 0 | length(listNnetModels()) > 0) {
        otherObjects <- c(listRpartModels(), listNnetModels())
        testDataSetOther <- tapply(otherObjects,
          as.factor(1:length(otherObjects)), parseDataSetOther)
        validOtherObjects <- otherObjects[testDataSetOther==.activeDataSet]
        }
    else {validOtherObjects <- character(0)}
    validModelObjects <- c(validGlmObjects, validOtherObjects)
    initializeDialog(title="Lift Chart")
    modelBox <- variableListBox(top, validModelObjects, selectmode="multiple",
      title=gettextRcmdr("Select One or More Models"))
    dataSetBox <- variableListBox(top, dataSets, selectmode="single",
      title=gettextRcmdr("Select a Single Data Set"),
      initialSelection=if (is.null(.activeDataSet)) NULL else which(.activeDataSet == dataSets) - 1)
    subsetBoxBCA()
    radioButtons(name="type", buttons=c("cumulative", "incremental"),
      labels=gettextRcmdr(c("Total Cumulative Response",
      "Incremental Response Rate")),
      title=gettextRcmdr("Select the Lift Chart Type"))
    responseFrame <- tkframe(top)
    if(is.null(getRcmdr(".trueResp"))) trueResp <- tclVar("")
    else trueResp <- tclVar(getRcmdr(".trueResp"))
    trueRespField <- tkentry(responseFrame, width="6", textvariable=trueResp)
    if(is.null(getRcmdr(".targetLevel"))) targLevel <- tclVar("")
    else targLevel <- tclVar(getRcmdr(".targetLevel"))
    targLevelField <- tkentry(responseFrame, width="10", textvariable=targLevel)
    sampNameFrame <- tkframe(top)
    if(is.null(getRcmdr(".Sample"))) sampName <- tclVar("")
    else sampName <- tclVar(getRcmdr(".Sample"))
    sampNameField <- tkentry(sampNameFrame, width="15", textvariable=sampName)
    onOK <- function() {
        modelObjcts <- getSelection(modelBox)
        dataSet <- getSelection(dataSetBox)
        subset <- tclvalue(subsetVariable)
        chartType <- tclvalue(typeVariable)
        trueRespRate <- tclvalue(trueResp)
        targLabel <- tclvalue(targLevel)
        sampLabel <- tclvalue(sampName)
        closeDialog()
        if(length(modelObjcts)==0) {
            errorCondition(recall=liftChart, 
              message=gettextRcmdr("No model objects have been selected."))
            return()
            }
        if(length(dataSet)==0) {
            errorCondition(recall=liftChart, 
              message=gettextRcmdr("No comparison data set has been selected."))
            return()
            }
        if(trueRespRate=="") {
            errorCondition(recall=liftChart, 
              message=gettextRcmdr("A true sample response rate must be provided."))
            return()
            }
        if(targLabel=="") {
            errorCondition(recall=liftChart, 
              message=gettextRcmdr("A target level label must be provided."))
            return()
            }
        if(is.equality.prob(trim.blanks(subset))) {
            errorCondition(recall=liftChart, 
            message=gettextRcmdr('Use "==" not "=" for equalties in subsetting.'))
            return()
            }
        if(length(modelObjcts) == 1) {
            modList <- paste('"', modelObjcts, '"', sep="")
            }
        else {
            modList <- paste('"', modelObjcts[1], '"', sep="")
            for(i in 2:length(modelObjcts)) {
                modList <- paste(modList, ', "', modelObjcts[i], '"', sep="")
                }
            }
        modList <- paste("c(", modList, ")", sep="")
        if((trim.blanks(subset) == gettextRcmdr("<all valid cases>")) |
          (trim.blanks(subset) == "")) {
             dataSet <- dataSet
             }
         else {
             dataSet <- paste(dataSet, "[", dataSet, "$",
               trim.blanks(subset),",]", sep="")
             putRcmdr(".subset", trim.blanks(subset))
             }
        putRcmdr(".targetLevel", targLabel)
        putRcmdr(".trueResp", trueRespRate)
        putRcmdr(".Sample", sampLabel)
        command1 <- paste(modList, ", ", dataSet, ", ",  '"', targLabel, '"',
          ", ", trueRespRate, ", ",  '"', chartType, '"',
          ", ", '"', sampLabel, '"', sep="")
        command <- paste("lift.chart(", command1, ")", sep="")
        justDoIt(command)
        logger(command)
        activateMenus()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="lift.chart")
    tkgrid(getFrame(modelBox), getFrame(dataSetBox), sticky="nw")
    tkgrid(typeFrame, subsetFrame, sticky="nw")
    tkgrid(tklabel(top, text=""), sticky="w")
    tkgrid(tklabel(responseFrame, text=gettextRcmdr("True Response Rate: ")),
      trueRespField, sticky="w")
    tkgrid(tklabel(responseFrame, text=gettextRcmdr("Target Level: ")),
      targLevelField, sticky="w")
    tkgrid(tklabel(sampNameFrame, text=gettextRcmdr("Sample Name: ")),
      sampNameField, sticky="w")
    tkgrid(responseFrame, sampNameFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=2, columns=2)
    invisible()
    }

liftChartP <- function() {
    parseDataSetGlm <- function(x) {
        y <- eval(parse(text=paste(x, "$call", sep="")))
        out <- unlist(strsplit(as.character(y), "="))[4]
        return(out)
        }
    parseDataSetOther <- function(x) {
        y <- eval(parse(text=paste(x, "$call", sep="")))
        out <- unlist(strsplit(as.character(y), "="))[3]
        return(out)
        }
    .activeDataSet <- ActiveDataSet()
    if(length(listGeneralizedLinearModels()>0)) {
        glmObjects <- listGeneralizedLinearModels()
        testDataSetGlm <- tapply(glmObjects, as.factor(1:length(glmObjects)),
          parseDataSetGlm)
        validGlmObjects <- glmObjects[testDataSetGlm==.activeDataSet]
        }
    else {validGlmObjects <- character(0)}
    if(length(listRpartModels()) > 0 | length(listNnetModels()) > 0) {
        otherObjects <- c(listRpartModels(), listNnetModels())
        testDataSetOther <- tapply(otherObjects,
          as.factor(1:length(otherObjects)), parseDataSetOther)
        validOtherObjects <- otherObjects[testDataSetOther==.activeDataSet]
        }
    else {validOtherObjects <- character(0)}
    validModelObjects <- c(validGlmObjects, validOtherObjects)
    test <- length(validModelObjects) > 0
    return(test)
    }

rankScoreGUI <- function() {
    parseDataSetGlm <- function(x) {
        y <- eval(parse(text=paste(x, "$call", sep="")))
        out <- unlist(strsplit(as.character(y), "="))[4]
        return(out)
        }
    parseDataSetOther <- function(x) {
        y <- eval(parse(text=paste(x, "$call", sep="")))
        out <- unlist(strsplit(as.character(y), "="))[3]
        return(out)
        }
    .activeDataSet <- ActiveDataSet()
    if(length(listGeneralizedLinearModels()>0)) {
        glmObjects <- listGeneralizedLinearModels()
        testDataSetGlm <- tapply(glmObjects, as.factor(1:length(glmObjects)),
          parseDataSetGlm)
        validGlmObjects <- glmObjects[testDataSetGlm==.activeDataSet]
        }
    else {validGlmObjects <- character(0)}
    if(length(listRpartModels()) > 0 | length(listNnetModels()) > 0) {
        otherObjects <- c(listRpartModels(), listNnetModels())
        testDataSetOther <- tapply(otherObjects,
          as.factor(1:length(otherObjects)), parseDataSetOther)
        validOtherObjects <- otherObjects[testDataSetOther==.activeDataSet]
        }
    else {validOtherObjects <- character(0)}
    validModelObjects <- c(validGlmObjects, validOtherObjects)
    initializeDialog(title="Rank Order Records Based on Expected Response")
    modelBox <- variableListBox(top, validModelObjects, selectmode="single",
      title=gettextRcmdr("Select a Single Model"))
    dataSetBox <- variableListBox(top, listDataSets(), selectmode="single",
      title=gettextRcmdr("Select a Single Data Set"))
    targetFrame <- tkframe(top)
    if(is.null(getRcmdr(".targetLevel"))) targLevel <- tclVar("")
    else targLevel <- tclVar(getRcmdr(".targetLevel"))
    targLevelField <- tkentry(targetFrame, width="10", textvariable=targLevel)
    varNameFrame<-tkframe(top)
    varName <- tclVar("Score")
    varNameField <- tkentry(varNameFrame, width="10", textvariable=varName)
    onOK <- function() {
        modObj <- getSelection(modelBox)
        dataSet <- getSelection(dataSetBox)
        targLabel <- tclvalue(targLevel)
        varLabel <- tclvalue(varName)
        closeDialog()
        if(length(modObj)==0) {
            errorCondition(recall=rankScoreGUI, 
              message=gettextRcmdr("No model has been selected."))
            return()
            }
        if(length(dataSet)==0) {
            errorCondition(recall=rankScoreGUI, 
              message=gettextRcmdr("No data set to score has been selected."))
            return()
            }
        if(targLabel=="") {
            errorCondition(recall=rankScoreGUI, 
              message=gettextRcmdr("A target level label must be provided."))
            return()
            }
        if (!is.valid.name(varLabel)){
            errorCondition(recall=rankScoreGUI, message=
	      paste('"', varLabel, '" ',
                  gettextRcmdr("is not a valid name."), sep=""))
            return()
            }
        c1 <- paste('"', modObj, '", ', dataSet, ', "', targLabel,'"', sep="")
        command <- paste(dataSet, "$", varLabel, " <- rankScore(", c1, ")", sep="")
        justDoIt(command)
        logger(command)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="score")
    tkgrid(getFrame(modelBox), getFrame(dataSetBox), sticky="nw")
    tkgrid(tklabel(targetFrame, text=gettextRcmdr("Target Level:")),
      targLevelField, sticky="w")
    tkgrid(tklabel(varNameFrame, text=gettextRcmdr("Score Variable Name: ")),
      varNameField, sticky="w")
    tkgrid(targetFrame, varNameFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=2, columns=2)
    return(invisible())
    }

rawProbScoreGUI <- function() {
    parseDataSetGlm <- function(x) {
        y <- eval(parse(text=paste(x, "$call", sep="")))
        out <- unlist(strsplit(as.character(y), "="))[4]
        return(out)
        }
    parseDataSetOther <- function(x) {
        y <- eval(parse(text=paste(x, "$call", sep="")))
        out <- unlist(strsplit(as.character(y), "="))[3]
        return(out)
        }
    .activeDataSet <- ActiveDataSet()
    if(length(listGeneralizedLinearModels()>0)) {
        glmObjects <- listGeneralizedLinearModels()
        testDataSetGlm <- tapply(glmObjects, as.factor(1:length(glmObjects)),
          parseDataSetGlm)
        validGlmObjects <- glmObjects[testDataSetGlm==.activeDataSet]
        }
    else {validGlmObjects <- character(0)}
    if(length(listRpartModels()) > 0 | length(listNnetModels()) > 0) {
        otherObjects <- c(listRpartModels(), listNnetModels())
        testDataSetOther <- tapply(otherObjects,
          as.factor(1:length(otherObjects)), parseDataSetOther)
        validOtherObjects <- otherObjects[testDataSetOther==.activeDataSet]
        }
    else {validOtherObjects <- character(0)}
    validModelObjects <- c(validGlmObjects, validOtherObjects)
    initializeDialog(title="Score Records Based on Sample Probabilities")
    modelBox <- variableListBox(top, validModelObjects, selectmode="single",
      title=gettextRcmdr("Select a Single Model"))
    dataSetBox <- variableListBox(top, listDataSets(), selectmode="single",
      title=gettextRcmdr("Select a Single Data Set"))
    targetFrame <- tkframe(top)
    if(is.null(getRcmdr(".targetLevel"))) targLevel <- tclVar("")
    else targLevel <- tclVar(getRcmdr(".targetLevel"))
    targLevelField <- tkentry(targetFrame, width="10", textvariable=targLevel)
    varNameFrame<-tkframe(top)
    varName <- tclVar("Score")
    varNameField <- tkentry(varNameFrame, width="10", textvariable=varName)
    onOK <- function() {
        modObj <- getSelection(modelBox)
        dataSet <- getSelection(dataSetBox)
        targLabel <- tclvalue(targLevel)
        varLabel <- tclvalue(varName)
        closeDialog()
        if(length(modObj)==0) {
            errorCondition(recall=rawProbScoreGUI, 
              message=gettextRcmdr("No model has been selected."))
            return()
            }
        if(length(dataSet)==0) {
            errorCondition(recall=rawProbScoreGUI, 
              message=gettextRcmdr("No data set to score has been selected."))
            return()
            }
        if(targLabel=="") {
            errorCondition(recall=rawProbScoreGUI, 
              message=gettextRcmdr("A target level label must be provided."))
            return()
            }
        if (!is.valid.name(varLabel)){
            errorCondition(recall=rawProbScoreGUI, message=
	      paste('"', varLabel, '" ',
                  gettextRcmdr("is not a valid name."), sep=""))
            return()
            }
        c1 <- paste('"', modObj, '", ', dataSet, ', "', targLabel,'"', sep="")
        command <- paste(dataSet, "$", varLabel, " <- rawProbScore(", c1, ")", sep="")
        justDoIt(command)
        logger(command)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="score")
    tkgrid(getFrame(modelBox), getFrame(dataSetBox), sticky="nw")
    tkgrid(tklabel(targetFrame, text=gettextRcmdr("Target Level:")),
      targLevelField, sticky="w")
    tkgrid(tklabel(varNameFrame, text=gettextRcmdr("Score Variable Name: ")),
      varNameField, sticky="w")
    tkgrid(targetFrame, varNameFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=2, columns=2)
    return(invisible())
    }

adjProbScoreGUI <- function() {
    parseDataSetGlm <- function(x) {
        y <- eval(parse(text=paste(x, "$call", sep="")))
        out <- unlist(strsplit(as.character(y), "="))[4]
        return(out)
        }
    parseDataSetOther <- function(x) {
        y <- eval(parse(text=paste(x, "$call", sep="")))
        out <- unlist(strsplit(as.character(y), "="))[3]
        return(out)
        }
    .activeDataSet <- ActiveDataSet()
    if(length(listGeneralizedLinearModels()>0)) {
        glmObjects <- listGeneralizedLinearModels()
        testDataSetGlm <- tapply(glmObjects, as.factor(1:length(glmObjects)),
          parseDataSetGlm)
        validGlmObjects <- glmObjects[testDataSetGlm==.activeDataSet]
        }
    else {validGlmObjects <- character(0)}
    if(length(listRpartModels()) > 0 | length(listNnetModels()) > 0) {
        otherObjects <- c(listRpartModels(), listNnetModels())
        testDataSetOther <- tapply(otherObjects,
          as.factor(1:length(otherObjects)), parseDataSetOther)
        validOtherObjects <- otherObjects[testDataSetOther==.activeDataSet]
        }
    else {validOtherObjects <- character(0)}
    validModelObjects <- c(validGlmObjects, validOtherObjects)
    initializeDialog(title="Score Records Based on Adjusted Probabilities")
    modelBox <- variableListBox(top, validModelObjects, selectmode="single",
      title=gettextRcmdr("Select a Single Model"))
    dataSetBox <- variableListBox(top, listDataSets(), selectmode="single",
      title=gettextRcmdr("Select a Single Data Set"))
    targetFrame <- tkframe(top)
    trueResp <- tclVar("")
    trueRespField <- tkentry(targetFrame, width="6", textvariable=trueResp)
    if(is.null(getRcmdr(".targetLevel"))) targLevel <- tclVar("")
    else targLevel <- tclVar(getRcmdr(".targetLevel"))
    targLevelField <- tkentry(targetFrame, width="10", textvariable=targLevel)
    varNameFrame<-tkframe(top)
    varName <- tclVar("Score")
    varNameField <- tkentry(varNameFrame, width="10", textvariable=varName)
    onOK <- function() {
        modObj <- getSelection(modelBox)
        dataSet <- getSelection(dataSetBox)
        trueRespRate <- tclvalue(trueResp)
        targLabel <- tclvalue(targLevel)
        varLabel <- tclvalue(varName)
        closeDialog()
        if(length(modObj)==0) {
            errorCondition(recall=adjProbScoreGUI, 
              message=gettextRcmdr("No model has been selected."))
            return()
            }
        if(length(dataSet)==0) {
            errorCondition(recall=adjProbScoreGUI, 
              message=gettextRcmdr("No data set to score has been selected."))
            return()
            }
        if(targLabel=="") {
            errorCondition(recall=adjProbScoreGUI, 
              message=gettextRcmdr("A target level label must be provided."))
            return()
            }
        if (!is.valid.name(varLabel)){
            errorCondition(recall=adjProbScoreGUI, message=
	      paste('"', varLabel, '" ',
                  gettextRcmdr("is not a valid name."), sep=""))
            return()
            }
        c1 <- paste('"', modObj, '", ', dataSet, ', "', targLabel,'", ', trueRespRate, sep="")
        command <- paste(dataSet, "$", varLabel, " <- adjProbScore(", c1, ")", sep="")
        justDoIt(command)
        logger(command)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="score")
    tkgrid(getFrame(modelBox), getFrame(dataSetBox), sticky="nw")
    tkgrid(tklabel(targetFrame, text=gettextRcmdr("True Response Rate: ")),
      trueRespField, sticky="w")
    tkgrid(tklabel(targetFrame, text=gettextRcmdr("Target Level:")),
      targLevelField, sticky="w")
    tkgrid(tklabel(varNameFrame, text=gettextRcmdr("Score Variable Name: ")),
      varNameField, sticky="w")
    tkgrid(targetFrame, varNameFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=2, columns=2)
    return(invisible())
    }

variablesP <- function(n=1) activeDataSetP() && length(listVariables()) >= n

# rglSaveP <- function() {
#    if(any(loadedNamespaces()!="rgl")) chk <- FALSE
#    else if(rgl.cur() == 0) chk <- FALSE
#    else chk <- TRUE
#    return(chk)
#    }

is.equality.prob <- function(subset) {
    problem <- FALSE
    chkString <- unlist(strsplit(subset, "="))
    if(length(chkString) == 2 & is.valid.name(trim.blanks(chkString[1]))) {
        problem <- TRUE
        }
    return(problem)
    }

# The "about" page for the plug-in
helpAboutBCA <- function() {
	if (as.numeric(R.Version()$major) >= 3) print(help("RcmdrPlugin.BCA"))
	else help("RcmdrPlugin.BCA")
}

# The change to the scatterPlot and scatterPlotMatrix functions to set the value
# of loess.threshold parameter to 2 from the default of 5. This allows for the
# plotting of discrete variables.
scatterPlotBCA <- function () {
	defaults <- list(initial.x = NULL, initial.y = NULL, initial.jitterx = 0, initial.jittery = 0, 
			initial.logstringx = 0, initial.logstringy = 0, initial.log = 0, initial.box = 1, 
			initial.line = 1, initial.smooth = 1, initial.spread = 1, initial.span = 50,
			initial.subset = gettextRcmdr ("<all valid cases>"), initial.ylab = gettextRcmdr ("<auto>"), 
			initial.xlab = gettextRcmdr("<auto>"), initial.pch = gettextRcmdr("<auto>"), 
			initial.cexValue = 1, initial.cex.axisValue = 1, initial.cex.labValue = 1, initialGroup=NULL, initial.lines.by.group=1) 
	dialog.values <- getDialog("scatterPlot", defaults)
	initial.group <- dialog.values$initial.group
	.linesByGroup <- if (dialog.values$initial.lines.by.group == 1) TRUE else FALSE
	.groups <- if (is.null(initial.group)) FALSE else initial.group
	initializeDialog(title = gettextRcmdr("Scatterplot"))
	.numeric <- Numeric()
	variablesFrame <- tkframe(top)
	xBox <- variableListBox(variablesFrame, .numeric, title = gettextRcmdr("x-variable (pick one)"), 
			initialSelection = varPosn (dialog.values$initial.x, "numeric"))
	yBox <- variableListBox(variablesFrame, .numeric, title = gettextRcmdr("y-variable (pick one)"), 
			initialSelection = varPosn (dialog.values$initial.y, "numeric"))
	optionsParFrame <- tkframe(top)
	checkBoxes(window = optionsParFrame, frame = "optionsFrame", 
			boxes = c("identify", "jitterX", "jitterY", "logX", "logY", 
					"boxplots", "lsLine", "smoothLine", "spread"), initialValues = c(dialog.values$initial.log, 
					dialog.values$initial.jitterx, dialog.values$initial.jittery, 
					dialog.values$initial.logstringx, dialog.values$initial.logstringy,
					dialog.values$initial.box, dialog.values$initial.line, dialog.values$initial.smooth,
					dialog.values$initial.spread),labels = gettextRcmdr(c("Identify points", 
							"Jitter x-variable", "Jitter y-variable", "Log x-axis", 
							"Log y-axis", "Marginal boxplots", "Least-squares line", 
							"Smooth line", "Show spread")), title = "Options")
	sliderValue <- tclVar(dialog.values$initial.span)
	slider <- tkscale(optionsFrame, from = 0, to = 100, showvalue = TRUE, 
			variable = sliderValue, resolution = 5, orient = "horizontal")
	subsetBox(subset.expression = dialog.values$initial.subset)
	labelsFrame <- tkframe(top)
	xlabVar <- tclVar(dialog.values$initial.xlab)
	ylabVar <- tclVar(dialog.values$initial.ylab)
	xlabFrame <- tkframe(labelsFrame)
	xlabEntry <- ttkentry(xlabFrame, width = "25", textvariable = xlabVar)
	xlabScroll <- ttkscrollbar(xlabFrame, orient = "horizontal", 
			command = function(...) tkxview(xlabEntry, ...))
	tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll, 
						...))
	tkgrid(labelRcmdr(xlabFrame, text = gettextRcmdr("x-axis label"), 
					fg = "blue"), sticky = "w")
	tkgrid(xlabEntry, sticky = "w")
	tkgrid(xlabScroll, sticky = "ew")
	ylabFrame <- tkframe(labelsFrame)
	ylabEntry <- ttkentry(ylabFrame, width = "25", textvariable = ylabVar)
	ylabScroll <- ttkscrollbar(ylabFrame, orient = "horizontal", 
			command = function(...) tkxview(ylabEntry, ...))
	tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll, 
						...))
	tkgrid(labelRcmdr(ylabFrame, text = gettextRcmdr("y-axis label"), 
					fg = "blue"), sticky = "w")
	tkgrid(ylabEntry, sticky = "w")
	tkgrid(ylabScroll, sticky = "ew")
	tkgrid(xlabFrame, labelRcmdr(labelsFrame, text = "     "), 
			ylabFrame, sticky = "w")
	parFrame <- tkframe(optionsParFrame)
	pchVar <- tclVar(dialog.values$initial.pch)
	pchEntry <- ttkentry(parFrame, width = 25, textvariable = pchVar)
	cexValue <- tclVar(dialog.values$initial.cexValue)
	cex.axisValue <- tclVar(dialog.values$initial.cex.axisValue)
	cex.labValue <- tclVar(dialog.values$initial.cex.labValue)
	cexSlider <- tkscale(parFrame, from = 0.5, to = 2.5, showvalue = TRUE, 
			variable = cexValue, resolution = 0.1, orient = "horizontal")
	cex.axisSlider <- tkscale(parFrame, from = 0.5, to = 2.5, 
			showvalue = TRUE, variable = cex.axisValue, resolution = 0.1, 
			orient = "horizontal")
	cex.labSlider <- tkscale(parFrame, from = 0.5, to = 2.5, 
			showvalue = TRUE, variable = cex.labValue, resolution = 0.1, 
			orient = "horizontal")
	onOK <- function() {
		x <- getSelection(xBox)
		y <- getSelection(yBox)
		jitter <- if ("1" == tclvalue(jitterXVariable) && "1" == 
						tclvalue(jitterYVariable)) 
					", jitter=list(x=1, y=1)"
				else if ("1" == tclvalue(jitterXVariable)) 
					", jitter=list(x=1)"
				else if ("1" == tclvalue(jitterYVariable)) 
					", jitter=list(y=1)"
				else ""
		logstring <- ""
		if ("1" == tclvalue(logXVariable)) 
			logstring <- paste(logstring, "x", sep = "")
		if ("1" == tclvalue(logYVariable)) 
			logstring <- paste(logstring, "y", sep = "")
		log <- tclvalue(identifyVariable)
		box <- tclvalue(boxplotsVariable)
		line <- tclvalue(lsLineVariable)
		smooth <-  tclvalue(smoothLineVariable)
		spread <- tclvalue(spreadVariable)
		span <- as.numeric(tclvalue(sliderValue))
		initial.subset <- subset <- tclvalue(subsetVariable)
		subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) 
					""
				else paste(", subset=", subset, sep = "")
		cex.axis <- as.numeric(tclvalue(cex.axisValue))
		cex <- as.numeric(tclvalue(cexValue))
		cex.lab <- as.numeric(tclvalue(cex.labValue))
		xlab <- trim.blanks(tclvalue(xlabVar))
		xlab <- if (xlab == gettextRcmdr("<auto>")) 
					""
				else paste(", xlab=\"", xlab, "\"", sep = "")
		ylab <- trim.blanks(tclvalue(ylabVar))
		ylab <- if (ylab == gettextRcmdr("<auto>")) 
					""
				else paste(", ylab=\"", ylab, "\"", sep = "")
		pch <- gsub(" ", ",", tclvalue(pchVar))
		putDialog ("scatterPlot", list (initial.x = x, initial.y = y, initial.jitterx = tclvalue(jitterXVariable),
						initial.jittery = tclvalue(jitterYVariable), initial.logstringx = tclvalue(logXVariable),
						initial.logstringy = tclvalue(logYVariable), initial.log = log, initial.box = box, 
						initial.line = line, initial.smooth = smooth, initial.spread = spread,
						initial.span = span, initial.subset = initial.subset, initial.xlab = tclvalue(xlabVar),
						initial.ylab = tclvalue(ylabVar), initial.cexValue = tclvalue(cexValue), 
						initial.cex.axisValue = tclvalue(cex.axisValue), initial.cex.labValue = tclvalue(cex.labValue), 
						initial.pch = pch, initial.group=if (.groups == FALSE) NULL else .groups,
						initial.lines.by.group=if (.linesByGroup) 1 else 0))
		closeDialog()
		if ("" == pch) {
			errorCondition(recall = scatterPlot, message = gettextRcmdr("No plotting characters."))
			return()
		}
		pch <- if (trim.blanks(pch) == gettextRcmdr("<auto>")) 
					""
				else paste(", pch=c(", pch, ")", sep = "")
		if (length(x) == 0 || length(y) == 0) {
			errorCondition(recall = scatterPlot, message = gettextRcmdr("You must select two variables"))
			return()
		}
		if (x == y) {
			errorCondition(recall = scatterPlot, message = gettextRcmdr("x and y variables must be different"))
			return()
		}
		.activeDataSet <- ActiveDataSet()
		log <- if (logstring != "") 
					paste(", log=\"", logstring, "\"", sep = "")
				else ""
		if ("1" == tclvalue(identifyVariable)) {
			RcmdrTkmessageBox(title = "Identify Points", message = paste(gettextRcmdr("Use left mouse button to identify points,\n"), 
							gettextRcmdr(if (MacOSXP()) 
												"esc key to exit."
											else "right button to exit."), sep = ""), icon = "info", 
					type = "ok")
			idtext <- ", id.method=\"identify\""
		}
		else idtext <- ""
		box <- if ("1" == tclvalue(boxplotsVariable)) 
					"'xy'"
				else "FALSE"
		line <- if ("1" == tclvalue(lsLineVariable)) 
					"lm"
				else "FALSE"
		smooth <- as.character("1" == tclvalue(smoothLineVariable))
		spread <- as.character("1" == tclvalue(spreadVariable))
		cex <- if (cex == 1) 
					""
				else paste(", cex=", cex, sep = "")
		cex.axis <- if (cex.axis == 1) 
					""
				else paste(", cex.axis=", cex.axis, sep = "")
		cex.lab <- if (cex.lab == 1) 
					""
				else paste(", cex.lab=", cex.lab, sep = "")
		if (.groups == FALSE) {
			doItAndPrint(paste("scatterplot(", y, "~", x, log, 
							", reg.line=", line, ", smooth=", smooth,
							", loess.threshold=2, spread=", 
							spread, idtext, ", boxplots=", box, ", span=", 
							span/100, jitter, xlab, ylab, cex, cex.axis, 
							cex.lab, pch, ", data=", .activeDataSet, subset, 
							")", sep = ""))
		}
		else {
			doItAndPrint(paste("scatterplot(", y, "~", x, " | ", 
							.groups, ", reg.line=", line, ", smooth=", smooth, 
							", loess.threshold=2, spread=", spread, idtext,
							", boxplots=", box,	", span=", span/100, jitter,
							xlab, ylab, cex, cex.axis, cex.lab, pch,
							", by.groups=", .linesByGroup, 
							", data=", .activeDataSet, subset, ")", sep = ""))
		}
		activateMenus()
		tkfocus(CommanderWindow())
	}
	groupsBox(scatterPlot, plotLinesByGroup = TRUE, initialGroup=initial.group, initialLinesByGroup=dialog.values$initial.lines.by.group,
			initialLabel=if (is.null(initial.group)) gettextRcmdr("Plot by groups") else paste(gettextRcmdr("Plot by:"), initial.group))
	OKCancelHelp(helpSubject = "scatterplot", reset = "scatterPlot")
	tkgrid(getFrame(xBox), getFrame(yBox), sticky = "nw")
	tkgrid(variablesFrame, sticky = "w")
	tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Span for smooth")), 
			slider, sticky = "w")
	tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Plotting Parameters"), 
					fg = "blue"), sticky = "w")
	tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Plotting characters")), 
			pchEntry, stick = "w")
	tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Point size")), 
			cexSlider, sticky = "w")
	tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Axis text size")), 
			cex.axisSlider, sticky = "w")
	tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Axis-labels text size")), 
			cex.labSlider, sticky = "w")
	tkgrid(optionsFrame, parFrame, sticky = "nw")
	tkgrid(optionsParFrame, sticky = "w")
	tkgrid(labelsFrame, sticky = "w")
	tkgrid(subsetFrame, sticky = "w")
	tkgrid(groupsFrame, sticky = "w")
	tkgrid(labelRcmdr(top, text = " "))
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	dialogSuffix(rows = 8, columns = 2)
}

scatterPlotMatrixBCA <- function () {
	defaults <- list(initial.variables = NULL, initial.line = 1, initial.smooth = 1, initial.spread = 0, 
			initial.span = 50, initial.diag = "density", initial.subset = gettextRcmdr ("<all valid cases>"),
			initialGroup=NULL, initial.lines.by.group=1) 
	dialog.values <- getDialog("scatterPlotMatrix", defaults)
	initial.group <- dialog.values$initial.group
	.linesByGroup <- if (dialog.values$initial.lines.by.group == 1) TRUE else FALSE
	.groups <- if (is.null(initial.group)) FALSE else initial.group
	initializeDialog(title = gettextRcmdr("Scatterplot Matrix"))
	variablesBox <- variableListBox(top, Numeric(), title = gettextRcmdr("Select variables (three or more)"), 
			selectmode = "multiple", initialSelection = varPosn (dialog.values$initial.variables, "numeric"))
	checkBoxes(frame = "optionsFrame", boxes = c("lsLine", "smoothLine", 
					"spread"), initialValues = c(dialog.values$initial.line, dialog.values$initial.smooth,
					dialog.values$initial.spread), labels = gettextRcmdr(c("Least-squares lines", 
							"Smooth lines", "Show spread")))
	sliderValue <- tclVar(dialog.values$initial.span)
	slider <- tkscale(optionsFrame, from = 0, to = 100, showvalue = TRUE, 
			variable = sliderValue, resolution = 5, orient = "horizontal")
	radioButtons(name = "diagonal", buttons = c("density", "histogram", 
					"boxplot", "oned", "qqplot", "none"), labels = gettextRcmdr(c("Density plots", 
							"Histograms", "Boxplots", "One-dimensional scatterplots", 
							"Normal QQ plots", "Nothing (empty)")), title = gettextRcmdr("On Diagonal"), 
			initialValue = dialog.values$initial.diag)
	subsetBox(subset.expression = dialog.values$initial.subset)
	onOK <- function() {
		variables <- getSelection(variablesBox)
		closeDialog()
		if (length(variables) < 3) {
			errorCondition(recall = scatterPlotMatrix, message = gettextRcmdr("Fewer than 3 variable selected."))
			return()
		}
		line <- if ("1" == tclvalue(lsLineVariable)) 
					"lm"
				else "FALSE"
		smooth <- as.character("1" == tclvalue(smoothLineVariable))
		spread <- as.character("1" == tclvalue(spreadVariable))
		span <- as.numeric(tclvalue(sliderValue))
		diag <- as.character(tclvalue(diagonalVariable))
		initial.subset <- subset <- tclvalue(subsetVariable)
		subset <- if (trim.blanks(subset) == gettextRcmdr("<all valid cases>")) ""
				else paste(", subset=", subset, sep="")
		putDialog("scatterPlotMatrix", list(initial.variables = variables, initial.line = tclvalue (lsLineVariable), 
						initial.smooth = tclvalue(smoothLineVariable),initial.spread = tclvalue (spreadVariable), 
						initial.span = span, initial.diag = diag, initial.subset = initial.subset, 
						initial.group=if (.groups == FALSE) NULL else .groups,
						initial.lines.by.group=if (.linesByGroup) 1 else 0))
		.activeDataSet <- ActiveDataSet()
		if (.groups == FALSE) {
			command <- paste("scatterplotMatrix(~", paste(variables, 
							collapse = "+"), ", reg.line=", line, ", smooth=", 
					smooth, ", loess.threshold=2, spread=", spread, ", span=", span/100, 
					", diagonal = '", diag, "', data=", .activeDataSet, 
					subset, ")", sep = "")
			logger(command)
			justDoIt(command)
		}
		else {
			command <- paste("scatterplotMatrix(~", paste(variables, 
							collapse = "+"), " | ", .groups, ", reg.line=", 
					line, ", smooth=", smooth, ", loess.threshold=2, spread=", spread, 
					", span=", span/100, ", diagonal= '", diag, "', by.groups=", 
					.linesByGroup, ", data=", .activeDataSet, subset, 
					")", sep = "")
			logger(command)
			justDoIt(command)
		}
		activateMenus()
		tkfocus(CommanderWindow())
	}
	groupsBox(scatterPlot, plotLinesByGroup = TRUE, initialGroup=initial.group, initialLinesByGroup=dialog.values$initial.lines.by.group,
			initialLabel=if (is.null(initial.group)) gettextRcmdr("Plot by groups") else paste(gettextRcmdr("Plot by:"), initial.group))
	OKCancelHelp(helpSubject = "scatterplotMatrix", reset = "scatterPlotMatrix")
	tkgrid(getFrame(variablesBox), sticky = "nw")
	tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Span for smooth")), 
			slider, sticky = "w")
	tkgrid(optionsFrame, sticky = "w")
	tkgrid(diagonalFrame, sticky = "w")
	tkgrid(subsetFrame, sticky = "w")
	tkgrid(groupsFrame, sticky = "w")
	tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
	dialogSuffix(rows = 6, columns = 2)
}

## Added as a workaround for R CMD check --as-cran to allow for an *.RD file
## example that would allow the package to not issue a warning for "No examples,
## no tests, no vignette"
# The neural network functions
nnSub <- function(data, subset) {
    rowInd <- 1:nrow(data)
    if(length(unlist(strsplit(subset, "=="))) == 2) {
        compType <- "Equal"
        subSplit <- unlist(strsplit(subset, "=="))
        subVar <- trim.blanks(subSplit[1])
        subVal <- eval(parse(text=trim.blanks(subSplit[2])))
        }
    else if(length(unlist(strsplit(subset, "!="))) == 2) {
        compType <- "NotEqual"
        subSplit <- unlist(strsplit(subset, "!="))
        subVar <- trim.blanks(subSplit[1])
        subVal <- eval(parse(text=trim.blanks(subSplit[2])))
        }
    else if(length(unlist(strsplit(subset, ">="))) == 2) {
        compType <- "GrtrEqual"
        subSplit <- unlist(strsplit(subset, ">="))
        subVar <- trim.blanks(subSplit[1])
        subVal <- eval(parse(text=trim.blanks(subSplit[2])))
        }
    else if(length(unlist(strsplit(subset, "<="))) == 2) {
        compType <- "LessEqual"
        subSplit <- unlist(strsplit(subset, "<="))
        subVar <- trim.blanks(subSplit[1])
        subVal <- eval(parse(text=trim.blanks(subSplit[2])))
        }
    else if(length(unlist(strsplit(subset, ">"))) == 2) {
        compType <- "Grtr"
        subSplit <- unlist(strsplit(subset, ">"))
        subVar <- trim.blanks(subSplit[1])
        subVal <- eval(parse(text=trim.blanks(subSplit[2])))
        }
    else if(length(unlist(strsplit(subset, "<"))) == 2) {
        compType <- "Less"
        subSplit <- unlist(strsplit(subset, "<"))
        subVar <- trim.blanks(subSplit[1])
        subVal <- eval(parse(text=trim.blanks(subSplit[2])))
        }
    else stop("No appropriate comparison operator was found.")
    subSet <- data[[subVar]]
    if(compType=="Equal") rowInd <- rowInd[subSet==subVal]
    else if(compType=="NotEqual") rowInd <- rowInd[subSet!=subVal]
    else if(compType=="GrtrEqual") rowInd <- rowInd[subSet>=subVal]
    else if(compType=="LessEqual") rowInd <- rowInd[subSet<=subVal]
    else if(compType=="Grtr") rowInd <- rowInd[subSet>subVal]
    else rowInd <- rowInd[subSet<subVal]
    rowInd <- rowInd[!is.na(rowInd)]
    return(rowInd)
    }

Nnet <- function (formula, data, decay, size, subset="") {
    set.seed(1)
    if(subset=="") {
        nnetObj <- nnet(formula=formula, data=data, decay=decay, size=size,
          maxit=400)
        for(i in 2:10) {
            set.seed(i)
            newNnet <- nnet(formula=formula, data=data, decay=decay, size=size,
              maxit=400)
            if(newNnet$value < nnetObj$value) nnetObj <- newNnet
            }
        }
    else {
        selectRow <- nnSub(data, subset)
        nData <- data[selectRow,]
        nnetObj <- nnet(formula=formula, data=nData, decay=decay, size=size,
          maxit=400)
        for(i in 2:10) {
            set.seed(i)
            newNnet <- nnet(formula=formula, data=nData, decay=decay, size=size,
              maxit=400)
            if(newNnet$value < nnetObj$value) nnetObj <- newNnet
            }
        }
    nnetObj$call <- call("Nnet", formula=formula, data=substitute(data), decay=decay,
      size=size, subset=subset)
    return(nnetObj)
    }
