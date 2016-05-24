# Rcmdr plugin seeg
# Miguel F. Acevedo

ecdfPlot <- function(){
# dialog to pick a variable from a dataset and call ecdf function of seeg
    top <- NULL; buttonsFrame <- NULL 
    initializeDialog(title= gettextRcmdr("Empirical Cumulative (ECDF)"))
    xBox <- variableListBox(top, Numeric(), title="Variable (pick one)")
    onOK <- function(){
        x <- getSelection(xBox)
        closeDialog()
        if (length(x) == 0){
            errorCondition(recall=ecdfPlot, message="You must select a variable")
            return()
            }
        command <- paste("plot(ecdf(", ActiveDataSet(), "$", x, '))', sep="")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="ecdf")
    tkgrid(getFrame(xBox), sticky="nw")    
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=2)
    }

tsPlot <- function(){
# dialog to pick a variable from a dataset and call ts.plot function
    top <- NULL; buttonsFrame <- NULL 
    initializeDialog(title="Time series Plot")
    xBox <- variableListBox(top, Numeric(), title="Variable (pick one)")
    onOK <- function(){
        x <- getSelection(xBox)
        closeDialog()
        if (length(x) == 0){
            errorCondition(recall=tsPlot, message="You must select a variable")
            return()
            }
        command <- paste("ts.plot(", ActiveDataSet(), "$", x, ')', sep="")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="ts.plot")
    tkgrid(getFrame(xBox), sticky="nw")    
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=2)
    }

acfPlot <- function(){
# dialog to pick a variable from a dataset and call acf function
    top <- NULL; buttonsFrame <- NULL 
    initializeDialog(title="Autocorrelation function Plot")
    xBox <- variableListBox(top, Numeric(), title="Variable (pick one)")
    onOK <- function(){
        x <- getSelection(xBox)
        closeDialog()
        if (length(x) == 0){
            errorCondition(recall=acfPlot, message="You must select a variable")
            return()
            }
        command <- paste("acf(", ActiveDataSet(), "$", x, ')', sep="")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="acf")
    tkgrid(getFrame(xBox), sticky="nw")    
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=2)
    }

edaPlot <- function(){
# dialog to pick a variable from a dataset and call eda6 function of seeg
    top <- NULL; buttonsFrame <- NULL 
    initializeDialog(title="Exploratory Analysis (EDA)")
    xBox <- variableListBox(top, Numeric(), title="Variable (pick one)")
    labelFrame <- tkframe(top)
    labelName <- tclVar("X")
    entrylabel <- tkentry(labelFrame, width="20", textvariable=labelName)
    pdfVariable <- tclVar("0")
    pdfCheckBox <- tkcheckbutton(labelFrame, variable=pdfVariable)
    onOK <- function(){
        x <- getSelection(xBox)
        closeDialog()
        if (length(x) == 0){
            errorCondition(recall=edaPlot, message="You must select a variable")
            return()
            }
	pdfvalue <- tclvalue(pdfVariable)
        labelNameValue <- as.character(trim.blanks(tclvalue(labelName)))
        command <- paste("eda6(", ActiveDataSet(), "$", x,',"', labelNameValue,'")', sep="")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="eda6")
    tkgrid(getFrame(xBox), sticky="nw")    
    tkgrid(tklabel(labelFrame, text="Variable name\nto label graphs:"), entrylabel, sticky="w")
    tkgrid(labelFrame,sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=2)

    }

scangeoEASppp <- function() {
# dialog to set a name for dataset to be read from a geoEAS file using scan.geoeas.ppp of seeg
    top <- NULL; buttonsFrame <- NULL 
    initializeDialog(title=gettextRcmdr("Read as geoEAS file"))
    optionsFrame <- tkframe(top)
    dsname <- tclVar(gettextRcmdr("Dataset"))
    entryDsname <- tkentry(optionsFrame, width="20", textvariable=dsname)
    onOK <- function(){
        dsnameValue <- trim.blanks(tclvalue(dsname))
        closeDialog()
        if (dsnameValue == ""){
            errorCondition(recall=scangeoEASppp,
                message=gettextRcmdr("You must enter a name for the data set."))
                return()
                }
        if (!is.valid.name(dsnameValue)){
            errorCondition(recall=scangeoEASppp,
                message=paste('"', dsnameValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
            }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
                scangeoEASppp()
                return()
                }
            }
        file <- tclvalue(tkgetOpenFile(filetypes=
            gettextRcmdr('{"Text Files" {".txt" ".TXT" ".dat" ".DAT" }} {"All Files" {"*"}}')))
        if (file == "") {
            if (getRcmdr("grab.focus")) tkgrab.release(top)
            tkdestroy(top)
            return()
            }
        command <- paste(dsnameValue, '<- scan.geoeas.ppp("',file,'")', sep="")
        doItAndPrint(command)
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="scan.geoeas.ppp")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Enter name for dataset:")), entryDsname, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1)
}

scangeoEASras <- function() {
# dialog to set a name for dataset to be read from a geoEAS file using scan.map.ras
    top <- NULL; buttonsFrame <- NULL 
    initializeDialog(title=gettextRcmdr("Read as geoEAS raster file"))
    optionsFrame <- tkframe(top)
    dsname <- tclVar(gettextRcmdr("Dataset"))
    entryDsname <- tkentry(optionsFrame, width="20", textvariable=dsname)
    onOK <- function(){
        dsnameValue <- trim.blanks(tclvalue(dsname))
        closeDialog()
        if (dsnameValue == ""){
            errorCondition(recall=scangeoEASras,
                message=gettextRcmdr("You must enter a name for the data set."))
                return()
                }
        if (!is.valid.name(dsnameValue)){
            errorCondition(recall=scangeoEASras,
                message=paste('"', dsnameValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
            }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
                scangeoEASras()
                return()
                }
            }
        file <- tclvalue(tkgetOpenFile(filetypes=
            gettextRcmdr('{"Text Files" {".txt" ".TXT" ".dat" ".DAT" }} {"All Files" {"*"}}')))
        if (file == "") {
            if (getRcmdr("grab.focus")) tkgrab.release(top)
            tkdestroy(top)
            return()
            }
        command <- paste(dsnameValue, '<- scan.map.ras("',file,'")', sep="")
        doItAndPrint(command)
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="scan.map.ras")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Enter name for dataset:")), entryDsname, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1)
}

convertToppp <- function(){
# dialog to select a dataset and convert to ppp using makeppp of seeg 
    top <- NULL; buttonsFrame <- NULL 
    initializeDialog(title=gettextRcmdr("Make a ppp (Planar Point Pattern)"))
    optionsFrame <- tkframe(top)
    dsname <- tclVar(gettextRcmdr("pppSet"))
    entryDsname <- tkentry(optionsFrame, width="20", textvariable=dsname)
    onOK <- function(){
        dsnameValue <- trim.blanks(tclvalue(dsname))
        closeDialog()
        if (dsnameValue == ""){
            errorCondition(recall=convertToppp,
                message=gettextRcmdr("You must enter a name for the ppp."))
                return()
                }
        command <- paste(dsnameValue, " <- makeppp(", ActiveDataSet(), ')', sep="")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="makeppp")
    tkgrid(tklabel(optionsFrame, text=gettextRcmdr("Enter name for ppp (it will not become active data set)")), entryDsname, sticky="w")
    tkgrid.configure(entryDsname, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=2)

}

quadChisqppp <- function(){
# dialog to select a dataset and perform quadrat analysis using Chisquare  
    top <- NULL; buttonsFrame <- NULL 
    initializeDialog(title="Quadrat Chisq Analysis")
    muVar <- tclVar("5")
    muEntry <- tkentry(top, width="6", textvariable=muVar)
    onOK <- function(){
        muValue <- as.numeric(tclvalue(muVar))
        closeDialog()
        command <- paste("quad.chisq.ppp(", ActiveDataSet(), ",", muValue, ')', sep="")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="quad.chisq.ppp")
    tkgrid(tklabel(top, text="target intensity"), muEntry, sticky="e")
    tkgrid.configure(muEntry, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=2)

}

quadPoisppp <- function() {
# dialog to select a dataset and perform quadrat analysis using Poisson  
    top <- NULL; buttonsFrame <- NULL 
    initializeDialog(title="Quadrat Poisson Analysis")
    muVar <- tclVar("0.1")
    muEntry <- tkentry(top, width="6", textvariable=muVar)
    onOK <- function(){
        muValue <- as.numeric(tclvalue(muVar))
        closeDialog()
        command <- paste("quad.poisson.ppp(", ActiveDataSet(), ",", muValue, ')', sep="")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="quad.poisson.ppp")
    tkgrid(tklabel(top, text="target intensity"), muEntry, sticky="e")
    tkgrid.configure(muEntry, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=2)

}

nnGKppp <- function() {
# invokes function nnGK.ppp of seeg to perform Nearest Neighbor analysis 
        command <- paste("nnGK.ppp(", ActiveDataSet(), ')', sep="")
        doItAndPrint(command)
}

nnGKenvppp <- function() {
# dialog to select number of Monte Carlo realizations and calculate envelopes for G and K 
    top <- NULL; buttonsFrame <- NULL 
    initializeDialog(title="Monte Carlo Nearest Neighbor")
    nsimVar <- tclVar("100")
    nsimEntry <- tkentry(top, width="6", textvariable=nsimVar)
    onOK <- function(){
        nsimValue <- as.numeric(tclvalue(nsimVar))
        closeDialog()
        command <- paste("nnGKenv.ppp(", ActiveDataSet(), ",", nsimValue, ')', sep="")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="nnGKenv.ppp")
    tkgrid(tklabel(top, text="Number Realizations"), nsimEntry, sticky="e")
    tkgrid.configure(nsimEntry, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=2)
}

varCal <- function() {
# dialog to set arguments for semivariogram calculations
    top <- NULL; buttonsFrame <- NULL 
    initializeDialog(title="Variogram Calculation")
    nlagsVar <- tclVar("10")
    nlagsEntry <- tkentry(top, width="6", textvariable=nlagsVar)
    typeVar <- tclVar("isotropic")
    isoButton <- tkradiobutton(top, variable=typeVar, value="isotropic")
    anisoButton <- tkradiobutton(top, variable=typeVar, value="anisotropic")
    thetaVar <- tclVar("0")
    thetaEntry <- tkentry(top, width="6", textvariable=thetaVar)
    dthetaVar <- tclVar("7.5")
    dthetaEntry <- tkentry(top, width="6", textvariable=dthetaVar)
    maxdistVar <- tclVar("10")
    maxdistEntry <- tkentry(top, width="6", textvariable=maxdistVar)
    dsnameVar <- tclVar("Dataset")
    dsnameEntry <- tkentry(top, width="6", textvariable=dsnameVar)
    onOK <- function(){
        nlagsValue <- as.numeric(tclvalue(nlagsVar))
        typeValue <- tclvalue(typeVar)
        thetaValue <- as.numeric(tclvalue(thetaVar))
        dthetaValue <- as.numeric(tclvalue(dthetaVar))
        maxdistValue <- as.numeric(tclvalue(maxdistVar))
        dsnameValue <- trim.blanks(tclvalue(dsnameVar))
        closeDialog()
        if (length(ActiveDataSet())==0){
            errorCondition(recall=varCal,
                message=gettextRcmdr("Select Active Data Set"))
                return()
                }
        if (dsnameValue == ActiveDataSet()){
            errorCondition(recall=varCal,
                message=gettextRcmdr("Variogram name should be different from the Active Data Set name"))
                return()
                }
        command <- paste(dsnameValue, " <- vario(", ActiveDataSet(), ",", nlagsValue, ',"',
                  typeValue, '",', thetaValue, ',', dthetaValue,',', maxdistValue,')', sep="")
        doItAndPrint(command)
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="vario")
    tkgrid(tklabel(top, text="Number lags"), nlagsEntry, sticky="e")
    tkgrid.configure(nlagsEntry, sticky="w")
    tkgrid(tklabel(top, text="Max Distance (determine based on domain)"), maxdistEntry,
           sticky="e")
    tkgrid.configure(maxdistEntry, sticky="w")
    tkgrid(tklabel(top, text="Dataset name to store the variogram"), dsnameEntry, sticky="e")
    tkgrid.configure(dsnameEntry, sticky="w")
    tkgrid(tklabel(top, text=gettextRcmdr("OmniDirectional")), isoButton, sticky="e")
    tkgrid(tklabel(top, text=gettextRcmdr("Directional")), anisoButton,
           tklabel(top, text="Theta"), thetaEntry,
           tklabel(top, text="Tol Theta"), dthetaEntry,
           sticky="e")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=2)
}

SemiVarCov <- function() {
# Dialog to set model of Semi-Variance and Covariance 
    top <- NULL; buttonsFrame <- NULL 
    initializeDialog(title="Model Semi-Variance and Covariance")
    nlagsVar <- tclVar("10")
    nlagsEntry <- tkentry(top, width="6", textvariable=nlagsVar)
    n0Var <- tclVar("0")
    n0Entry <- tkentry(top, width="6", textvariable=n0Var)
    c0Var <- tclVar("1")
    c0Entry <- tkentry(top, width="6", textvariable=c0Var)
    rangeVar <- tclVar("1")
    rangeEntry <- tkentry(top, width="6", textvariable=rangeVar)
    onOK <- function(){
        nlagsValue <- as.numeric(tclvalue(nlagsVar))
        n0Value <- as.numeric(tclvalue(n0Var))
        c0Value <- as.numeric(tclvalue(c0Var))
        rangeValue <- as.numeric(tclvalue(rangeVar))
        closeDialog()
        if (length(ActiveDataSet())==0){
            errorCondition(recall=SemiVarCov,
                message=gettextRcmdr("Select a variogram as Active Data Set"))
                return()
                }
        command <- paste("model.semivar.cov(", ActiveDataSet(), ",",
                         nlagsValue, ',', n0Value, ',', c0Value, ',', rangeValue,')', sep="")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="model.semivar.cov")
    tkgrid(tklabel(top, text="Number lags"), nlagsEntry, sticky="e")
    tkgrid.configure(nlagsEntry, sticky="w")
    tkgrid(tklabel(top, text="Nugget"), n0Entry, sticky="e")
    tkgrid.configure(n0Entry, sticky="w")
    tkgrid(tklabel(top, text="Sill"), c0Entry, sticky="e")
    tkgrid.configure(c0Entry, sticky="w")
    tkgrid(tklabel(top, text="Range"), rangeEntry, sticky="e")
    tkgrid.configure(rangeEntry, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=2)
}

MakeSemiVarCov <- function() {
# Dialog to make a model of Semi-Variance and Covariance 
    top <- NULL; buttonsFrame <- NULL; .x <- NULL 
    initializeDialog(title="Make Spherical Model of Semi-Variance and Covariance")
    n0Var <- tclVar("0")
    n0Entry <- tkentry(top, width="6", textvariable=n0Var)
    c0Var <- tclVar("1")
    c0Entry <- tkentry(top, width="6", textvariable=c0Var)
    rangeVar <- tclVar("1")
    rangeEntry <- tkentry(top, width="6", textvariable=rangeVar)
    dsnameVar <- tclVar("VarModelName")
    dsnameEntry <- tkentry(top, width="6", textvariable=dsnameVar)
    onOK <- function(){
        n0Value <- as.numeric(tclvalue(n0Var))
        c0Value <- as.numeric(tclvalue(c0Var))
        rangeValue <- as.numeric(tclvalue(rangeVar))
        dsnameValue <- trim.blanks(tclvalue(dsnameVar))
        closeDialog()
        if (dsnameValue == ActiveDataSet()){
            errorCondition(recall=varCal,
                message=gettextRcmdr("Variogram model name should be different from the Active Data Set name"))
                return()
                }
        command <- paste(dsnameValue, " <- make.variogram(",n0Value, ',',
                  c0Value, ',', rangeValue,')', sep="")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="make.variogram")
    tkgrid(tklabel(top, text="Nugget"), n0Entry, sticky="e")
    tkgrid.configure(n0Entry, sticky="w")
    tkgrid(tklabel(top, text="Sill"), c0Entry, sticky="e")
    tkgrid.configure(c0Entry, sticky="w")
    tkgrid(tklabel(top, text="Range"), rangeEntry, sticky="e")
    tkgrid.configure(rangeEntry, sticky="w")
    tkgrid(tklabel(top, text="Name to store the variogram model"), dsnameEntry, sticky="e")
    tkgrid.configure(dsnameEntry, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=2)
}

OKrigingCmd <- function(){
# dialog to set arguments for Ordinary Kriging
    top <- NULL; buttonsFrame <- NULL 
    initializeDialog(title="Ordinary Kriging")
    dsnameSemiVar <- tclVar("Semivariance")
    dsnameSemiEntry <- tkentry(top, width="15", textvariable=dsnameSemiVar)
    nstepVar <- tclVar("1")
    nstepEntry <- tkentry(top, width="6", textvariable=nstepVar)
    maxdistVar <- tclVar("10")
    maxdistEntry <- tkentry(top, width="6", textvariable=maxdistVar)
    dsnameVar <- tclVar("Dataset")
    dsnameEntry <- tkentry(top, width="15", textvariable=dsnameVar)
    onOK <- function(){
        dsnameSemiValue <- trim.blanks(tclvalue(dsnameSemiVar))
        nstepValue <- as.numeric(tclvalue(nstepVar))
        maxdistValue <- as.numeric(tclvalue(maxdistVar))
        dsnameValue <- trim.blanks(tclvalue(dsnameVar))
        closeDialog()
        command <- paste(dsnameValue, " <- Okriging(", ActiveDataSet(), ",",
                          dsnameSemiValue, ",", nstepValue, ",", maxdistValue,')', sep="")
        doItAndPrint(command)
        activeDataSet(dsnameValue)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="Okriging")
    tkgrid(tklabel(top, text="Dataset name of Semivariance model"), dsnameSemiEntry, sticky="e")
    tkgrid.configure(dsnameSemiEntry, sticky="w")
    tkgrid(tklabel(top, text="Step"), nstepEntry, sticky="e")
    tkgrid.configure(nstepEntry, sticky="w")
    tkgrid(tklabel(top, text="Max Distance (determine based on domain)"), maxdistEntry,
           sticky="e")
    tkgrid.configure(maxdistEntry, sticky="w")
    tkgrid(tklabel(top, text="Dataset name to store the kriged results"), dsnameEntry, sticky="e")
    tkgrid.configure(dsnameEntry, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=2)
}

PlotKriged <- function(){
# dialog to set arguments to plot a dataset together with its kriged map
    top <- NULL; buttonsFrame <- NULL 
    initializeDialog(title="Plot Kriged")
    dsnameKriged <- tclVar(gettextRcmdr("KrigedDataSet"))
    dsnameKrigedEntry <- tkentry(top, width="15", textvariable=dsnameKriged)
    filename <- tclVar(gettextRcmdr('"out.pdf"'))
    filenameEntry <- tkentry(top, width="20", textvariable=filename)
    onOK <- function(){
        dsnameKrigedValue <- trim.blanks(tclvalue(dsnameKriged))
        filenameValue <- trim.blanks(tclvalue(filename))
        closeDialog()
        if (dsnameKrigedValue == ""){
            errorCondition(recall=PlotKriged,
                message=gettextRcmdr("You must enter a name for the data set."))
                return()
                }
        if (!is.valid.name(dsnameKrigedValue)){
            errorCondition(recall=PlotKriged,
                message=paste('"', dsnameKrigedValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
            }
        command <- paste("plotkriged(", ActiveDataSet(), ",",
                         dsnameKrigedValue,",", filenameValue,')', sep="")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="plotkriged")
    tkgrid(tklabel(top, text="Dataset name of kriged results"), dsnameKrigedEntry, sticky="e")
    tkgrid.configure(dsnameKrigedEntry, sticky="w")
    tkgrid(tklabel(top, text="Name of PDF out file"), filenameEntry, sticky="e")
    tkgrid.configure(filenameEntry, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=2)
}

