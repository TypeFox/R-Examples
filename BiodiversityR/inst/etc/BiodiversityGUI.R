# BiodiversityR
# Biodiversity Analysis Functions for R
# Developed by Roeland Kindt
#
# This software accompanies Kindt R and Coe R. 2005. Tree Diversity Analysis. A manual 
# and software for some common statistical methods for biodiversity and ecological 
# analysis. World Agroforestry Centre (ICRAF), Nairobi. vi+196 pp. This is also the suggested citation
# for this software.
#
# Many of the user interface functions were based on John Fox's R Commander (Rcmdr)
# Other functions are often based on the vegan package (Oksanen, Kindt, Legendre & O'Hara)
#
# To use most of the BiodiversityR functions, a community and an environmental dataset need to be identified first.
# Both datasets have the same number of rows (= number of sample units)
# Columns in the community dataset usually represent species
# Columns in the environmental represent environmental characteristics of sample units
#
# Roeland Kindt takes no liability for any direct, special, indirect or consequential damages resulting
# from loss of use, data or profits arising in connection with the use or performance of this software.
#
# The software can be quoted or reproduced without charge, provided the source is
# acknowledged. You must adhere to conditions of copyright of R software documented in
# rw2011\COPYING. Use the citation() or loaded.citations() function for acknowledgments in publications for 
# any package that you made use of.

putRcmdr(".communityDataSet", NULL)
putRcmdr("operatorFont", tkfont.create(family="courier", size=getRcmdr("log.font.size")))


#changed based on R Commander 1.9-6 (data-Menu.R)

selectCommunityDataSet <- function(){
	dataSets <- listDataSets()
	.communityDataSet <- CommunityDataSet()
	if ((length(dataSets) == 1) && !is.null(.communityDataSet)) {
		Message(message=gettextRcmdr("There is only one dataset in memory."),
				type="warning")
		tkfocus(CommanderWindow())
		return()
	}
	if (length(dataSets) == 0){
		Message(message=gettextRcmdr("There are no data sets from which to choose."),
				type="error")
		tkfocus(CommanderWindow())
		return()
	}
	initializeDialog(title=gettextRcmdr("Select Community Data Set"))
	dataSetsBox <- variableListBox(top, dataSets, title=gettextRcmdr("Data Sets (pick one)"),
			initialSelection=if (is.null(.communityDataSet)) NULL else which(.communityDataSet == dataSets) - 1)
	onOK <- function(){
		communityDataSet(getSelection(dataSetsBox))
		closeDialog()
		tkfocus(CommanderWindow())
	}
	OKCancelHelp()
	tkgrid(getFrame(dataSetsBox), sticky="nw")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=2, columns=1)
}


#changed based on R Commander 1.9-6 (utilities.R)


communityDataSet <- function(dsname, flushModel=TRUE, flushDialogMemory=TRUE){
	.communityDataSet <- CommunityDataSet()
	if (missing(dsname)) {
		if (is.null(.communityDataSet)){
			Message(message=gettextRcmdr("There is no community data set."), type="error")
			return(FALSE)
		}
		else return(.communityDataSet)
	}
	if (!is.data.frame(ds <- get(dsname, envir=.GlobalEnv))){
		if (!exists.method("as.data.frame", ds, default=FALSE)){
			Message(message=paste(dsname, gettextRcmdr(" is not a data frame and cannot be attached."),
							sep=""), type="error")
			tkfocus(CommanderWindow())
			return()
		}
		command <- paste(dsname, " <- as.data.frame(", dsname, ")", sep="")
		justDoIt(command)
		logger(command)
		Message(message=paste(dsname, gettextRcmdr(" has been coerced to a data frame"), sep=""),
				type="warning")
	}
	varnames <- names(get(dsname, envir=.GlobalEnv))
	newnames <- make.names(varnames)
	badnames <- varnames != newnames
	if (any(badnames)){
		command <- paste("names(", dsname, ") <- make.names(names(",
				dsname, "))", sep="")
		doItAndPrint(command)
	}
	if (!is.null(.communityDataSet) && getRcmdr("attach.data.set")
			&& (length(grep(.communityDataSet, search())) !=0)) {
		detach(pos = match(.communityDataSet, search()))
		logger(paste("detach(", .communityDataSet, ")", sep=""))
	}
	if (flushModel) {
		putRcmdr(".activeModel", NULL)
		RcmdrTclSet("modelName", gettextRcmdr("<No active model>"))
		if (!is.SciViews()) tkconfigure(getRcmdr("modelLabel"), foreground="red") else refreshStatus()
	}
	if (flushDialogMemory) putRcmdr("dialog.values", list())
	# -PhG tkconfigure(.modelLabel, foreground="red")
	CommunityDataSet(dsname)
	Message(sprintf(gettextRcmdr("The dataset %s has %d rows and %d columns."), dsname,
					nrow(get(dsname, envir=.GlobalEnv)), ncol(get(dsname, envir=.GlobalEnv))), type="note")
	if (any(badnames)) Message(message=paste(dsname, gettextRcmdr(" contains non-standard variable names:\n"),
						paste(varnames[badnames], collapse=", "),
						gettextRcmdr("\nThese have been changed to:\n"), paste(newnames[badnames], collapse=", "),
						sep=""), type="warning")
	CVariables(listCVariables())
#	Numeric(listNumeric())
#	Factors(listFactors())
#	TwoLevelFactors(listTwoLevelFactors())

# changed 2014 all to #
#	RcmdrTclSet("dataSetName", paste(" ", dsname, " "))
	# -PhG tkconfigure(.dataSetLabel, foreground="blue")
#	if (!is.SciViews()) tkconfigure(getRcmdr("dataSetLabel"), foreground="blue") else refreshStatus() # +PhG
# 	if (getRcmdr("attach.data.set")){
# 		attach(get(dsname, envir=.GlobalEnv), name=dsname)
# 		logger(paste("attach(", dsname, ")", sep=""))
# 	}
#	if (is.SciViews()) refreshStatus() else if (flushModel) tkconfigure(getRcmdr("modelLabel"), foreground="red") # +PhG (& J.Fox, 25Dec04)
	activateMenus()
	dsname
}




# changed based on R Commander 1.9-6 (utilities.R)

checkCommunityDataSet <- function(){
	if (communityDataSet() == FALSE) {
		tkfocus(CommanderWindow())
		FALSE
	}
	else TRUE
}


# changed based on R Commander 1.9-6 (utilities.R)

CommunityDataSet <- function(name){
	if (missing(name)) {
		temp <- getRcmdr(".communityDataSet")
		if (is.null(temp))
			return(NULL)
		else
		if (!exists(temp) || !is.data.frame(get(temp,envir=.GlobalEnv))) {
			Message(sprintf(gettextRcmdr("the dataset %s is no longer available"),
							temp), type="error")
			putRcmdr(".communityDataSet", NULL)
			RcmdrTclSet("dataSetName", gettextRcmdr("<No active dataset>"))
			putRcmdr(".activeModel", NULL)
			RcmdrTclSet("modelName", gettextRcmdr("<No active model>"))
			if (!is.SciViews()) {
				tkconfigure(getRcmdr("dataSetLabel"), foreground="red") 
				tkconfigure(getRcmdr("modelLabel"), foreground="red") 
			} 
			else refreshStatus()
			activateMenus()
			if (getRcmdr("suppress.menus") && RExcelSupported()) return(NULL)
		}
		return(temp)
	}
	else putRcmdr(".communityDataSet", name)
}


CVariables <- function(cnames){
    if (missing(cnames)) getRcmdr("cvariables")
    else putRcmdr("cvariables", cnames)
    }

listCVariables <- function(dataSet=CommunityDataSet()) {
    cvars <- eval(parse(text=paste("names(", dataSet,")")), envir=.GlobalEnv)
    if (getRcmdr("sort.names")) sort(cvars) else cvars
    }

communityDataSetP <- function() !is.null(CommunityDataSet())

data(dune)
data(dune.env)
dune2 <- dune
dune.env2 <- dune.env
seq <- c(2,13,4,16,6,1,8,5,17,15,10,11,9,18,3,20,14,19,12,7)
dune2[seq,] <- dune[1:20,]
dune.env2[seq,] <- dune.env[1:20,]
rownames(dune2)[seq] <- rownames(dune)[1:20]
rownames(dune.env2)[seq] <- rownames(dune)[1:20]
dune3 <- dune2
seq <- order(colnames(dune2))
dune3[,1:30] <- dune2[,seq]
colnames(dune3)[1:30] <- colnames(dune2)[seq]
dune <- dune3
dune.env <- dune.env2
rownames(dune) <- rownames(dune.env) <- c("X01","X02","X03","X04","X05","X06","X07","X08","X09","X10",
    "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20")
remove(dune2,dune3,dune.env2)


makecommunityGUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Make community matrix")
    .activeDataSet <- ActiveDataSet()
    .fvariables <- Factors()
    fvariables <- paste(.fvariables, ifelse(is.element(.fvariables, Factors()), "[factor]", ""))
    .nvariables <- Numeric()
    nvariables <- paste(.nvariables)
    modelName <- tclVar("Community.1")
    modelFrame <- tkframe(top, relief="groove", borderwidth=2)
    model <- tkentry(modelFrame, width=40, textvariable=modelName)
    siteFrame <- tkframe(top, relief="groove", borderwidth=2)
    siteBox <- tklistbox(siteFrame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    siteScroll <- tkscrollbar(siteFrame, repeatinterval=5, command=function(...) tkyview(siteBox, ...))
    tkconfigure(siteBox, yscrollcommand=function(...) tkset(siteScroll, ...))
    for (x in fvariables) tkinsert(siteBox, "end", x)
    specFrame <- tkframe(top, relief="groove", borderwidth=2)
    specBox <- tklistbox(specFrame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    specScroll <- tkscrollbar(specFrame, repeatinterval=5, command=function(...) tkyview(specBox, ...))
    tkconfigure(specBox, yscrollcommand=function(...) tkset(specScroll, ...))
    for (x in fvariables) tkinsert(specBox, "end", x)
    valueFrame <- tkframe(top, relief="groove", borderwidth=2)
    valueBox <- tklistbox(valueFrame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    valueScroll <- tkscrollbar(valueFrame, repeatinterval=5, command=function(...) tkyview(valueBox, ...))
    tkconfigure(valueBox, yscrollcommand=function(...) tkset(valueScroll, ...))
    for (x in nvariables) tkinsert(valueBox, "end", x)
    subsetFrame <- tkframe(top, relief="groove", borderwidth=2)
    subset1Frame <- tkframe(subsetFrame)
    subset2Frame <- tkframe(subsetFrame)
    subsetBox <- tklistbox(subset1Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    subsetScroll <- tkscrollbar(subset1Frame, repeatinterval=5, command=function(...) tkyview(subsetBox, ...))
    tkconfigure(subsetBox, yscrollcommand=function(...) tkset(subsetScroll, ...))
    variables <- c("all",fvariables)
    for (x in variables) tkinsert(subsetBox, "end", x)
    subset <- tclVar("")
    subsetEntry <- tkentry(subset2Frame, width=10, textvariable=subset)
    onOK <- function(){
        modelValue <- tclvalue(modelName)
        site <- .fvariables[as.numeric(tkcurselection(siteBox))+1]
        spec <- .fvariables[as.numeric(tkcurselection(specBox))+1]
        value <- .nvariables[as.numeric(tkcurselection(valueBox))+1]
        var <- variables[as.numeric(tkcurselection(subsetBox))+1]
        sub <- tclvalue(subset)
        if (var == "all") {
            command <- paste("makecommunitydataset(", .activeDataSet, ", row='", site, "', column='", spec, "', value='", value, "')", sep="")
        }else{
            var <- .fvariables[as.numeric(tkcurselection(subsetBox))]
            command <- paste("makecommunitydataset(", .activeDataSet, ", row='", site, "', column='", spec, "', value='", value, "', factor='", var, "', level='", sub, "')", sep="")
        }
        logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        communityDataSet(modelValue)       
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('makecommunitydataset', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(modelFrame, text="Save result as:    ", width=15), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(tklabel(siteFrame, text="Site variable (rows)"), sticky="w")
    tkgrid(siteBox, siteScroll, sticky="w")
    tkgrid(siteFrame, sticky="w")
    tkgrid(tklabel(specFrame, text="Species variable (columns)"), sticky="w")
    tkgrid(specBox, specScroll, sticky="w")
    tkgrid(specFrame, sticky="w")
    tkgrid(tklabel(valueFrame, text="Abundance variable"), sticky="w")
    tkgrid(valueBox, valueScroll, sticky="w")
    tkgrid(valueFrame, sticky="w")
    tkgrid(tklabel(subsetFrame, text="Subset options"), sticky="w")
    tkgrid(subsetBox, subsetScroll, sticky="w")
    tkgrid(tklabel(subset2Frame, text="subset:                ", width=15), subsetEntry, sticky="w")
    tkgrid(subset1Frame, sticky="w")
    tkgrid(subset2Frame, sticky="w")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(subsetScroll, sticky="ns")
    tkgrid.configure(siteScroll, sticky="ns")
    tkgrid.configure(specScroll, sticky="ns")
    tkgrid.configure(valueScroll, sticky="ns")
    tkselection.set(subsetBox, 0)
    tkselection.set(siteBox, 0)
    tkselection.set(specBox, 0)
    tkselection.set(valueBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(siteBox)
    tkwait.window(top)
}


importfromExcelGUI <- function() {
    initializeDialog(title="Read Community and Environmental data From Excel")
    optionsFrame <- tkframe(top, relief="groove", borderwidth=2)
    comdsname <- tclVar("CommunityData")
    entrycomDsname <- tkentry(optionsFrame, width="20", textvariable=comdsname)
    envdsname <- tclVar("EnvironmentalData")
    entryenvDsname <- tkentry(optionsFrame, width="20", textvariable=envdsname)
    dsites <- tclVar("sites")
    entrysites <- tkentry(optionsFrame, width="20", textvariable=dsites)
    stackedFrame <- tkframe(top, relief="groove", borderwidth=2)
    stackedVariable <- tclVar("0")
    stackedCheckBox <- tkcheckbutton(stackedFrame, variable=stackedVariable)
    scolumn <- tclVar("species")
    entrycol <- tkentry(stackedFrame, width="20", textvariable=scolumn)
    sval <- tclVar("abundance")
    entryval <- tkentry(stackedFrame, width="20", textvariable=sval)
    sfactor <- tclVar("all")
    entryfactor <- tkentry(stackedFrame, width="20", textvariable=sfactor)
    slevel <- tclVar("all")
    entrylevel <- tkentry(stackedFrame, width="20", textvariable=slevel)
    onOK <- function(){
        closeDialog()
        comdsnameValue <- tclvalue(comdsname)
        envdsnameValue <- tclvalue(envdsname)
        sitesValue <- tclvalue(dsites)
        colValue <- tclvalue(scolumn)
        valValue <- tclvalue(sval)
        factorValue <- tclvalue(sfactor)
        levelValue <- tclvalue(slevel)
        file <- tclvalue(tkgetOpenFile(filetypes='{"Excel Files" {".xls" ".XLS"}} {"All Files" {"*"}}'))
        if (file == "") {
            if (getRcmdr("grab.focus")) tkgrab.release(top)
            tkdestroy(top)
            return()
            }
        justDoIt(paste("library(RODBC)", sep=""))
        logger(paste("library(RODBC)", sep=""))
        stacked <- tclvalue(stackedVariable) == "1"
        if (stacked==F) {
            command <- paste("import.from.Excel('", file, "', data.type='community', sheet='community', sitenames='", sitesValue, "', cepnames=F, write.csv=F, csv.file='community.csv')", sep="")
        }else{
            if (factorValue=="all") { 
                command <- paste("import.from.Excel('", file, "', data.type='stacked', sheet='stacked', sitenames='", sitesValue, "', column='", colValue, "', value='", valValue, "', cepnames=F, write.csv=F, csv.file='community.csv')", sep="")
            }else{
                command <- paste("import.from.Excel('", file, "', data.type='stacked', sheet='stacked', sitenames='", sitesValue, "', column='", colValue, "', value='", valValue, "', factor='", factorValue, "', level='", levelValue, "', cepnames=F, write.csv=F, csv.file='community.csv')", sep="")
            }
        }
        logger(paste(comdsnameValue, " <- ", command, sep=""))
        assign(comdsnameValue, justDoIt(command), envir=.GlobalEnv)
        communityDataSet(comdsnameValue)
        command <- paste("import.from.Excel('", file, "', data.type='environmental', sheet='environmental', sitenames='", sitesValue, "', write.csv=F, csv.file='environmental.csv')", sep="")
        logger(paste(envdsnameValue, " <- ", command, sep=""))
        assign(envdsnameValue, justDoIt(command), envir=.GlobalEnv)
        activeDataSet(envdsnameValue)
        tkfocus(CommanderWindow())
        }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('import.from.Excel', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(optionsFrame, text="Names for new datasets"), sticky="w")
    tkgrid(tklabel(optionsFrame, text="Enter name for community data set:"), entrycomDsname, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Enter name for environmental data set:"), entryenvDsname, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Enter name for variable with sites:"), entrysites, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Options for stacked data entry"), sticky="w")
    tkgrid(tklabel(stackedFrame, text="Import community dataset from stacked format:"), stackedCheckBox, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Enter variable for species:"), entrycol, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Enter variable for abundance:"), entryval, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Enter factor for subset:"), entryfactor, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Enter level for subset:"), entrylevel, sticky="w")
    tkgrid(stackedFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1)
}

importfromExcel2007GUI <- function() {
    initializeDialog(title="Read Community and Environmental data From Excel 2007")
    optionsFrame <- tkframe(top, relief="groove", borderwidth=2)
    comdsname <- tclVar("CommunityData")
    entrycomDsname <- tkentry(optionsFrame, width="20", textvariable=comdsname)
    envdsname <- tclVar("EnvironmentalData")
    entryenvDsname <- tkentry(optionsFrame, width="20", textvariable=envdsname)
    dsites <- tclVar("sites")
    entrysites <- tkentry(optionsFrame, width="20", textvariable=dsites)
    stackedFrame <- tkframe(top, relief="groove", borderwidth=2)
    stackedVariable <- tclVar("0")
    stackedCheckBox <- tkcheckbutton(stackedFrame, variable=stackedVariable)
    scolumn <- tclVar("species")
    entrycol <- tkentry(stackedFrame, width="20", textvariable=scolumn)
    sval <- tclVar("abundance")
    entryval <- tkentry(stackedFrame, width="20", textvariable=sval)
    sfactor <- tclVar("all")
    entryfactor <- tkentry(stackedFrame, width="20", textvariable=sfactor)
    slevel <- tclVar("all")
    entrylevel <- tkentry(stackedFrame, width="20", textvariable=slevel)
    onOK <- function(){
        closeDialog()
        comdsnameValue <- tclvalue(comdsname)
        envdsnameValue <- tclvalue(envdsname)
        sitesValue <- tclvalue(dsites)
        colValue <- tclvalue(scolumn)
        valValue <- tclvalue(sval)
        factorValue <- tclvalue(sfactor)
        levelValue <- tclvalue(slevel)
        file <- tclvalue(tkgetOpenFile(filetypes='{"Excel Files" {".xlsx" ".XLSX"}} {"All Files" {"*"}}'))
        if (file == "") {
            if (getRcmdr("grab.focus")) tkgrab.release(top)
            tkdestroy(top)
            return()
            }
        justDoIt(paste("library(RODBC)", sep=""))
        logger(paste("library(RODBC)", sep=""))
        stacked <- tclvalue(stackedVariable) == "1"
        if (stacked==F) {
            command <- paste("import.from.Excel2007('", file, "', data.type='community', sheet='community', sitenames='", sitesValue, "', cepnames=F, write.csv=F, csv.file='community.csv')", sep="")
        }else{
            if (factorValue=="all") { 
                command <- paste("import.from.Excel2007('", file, "', data.type='stacked', sheet='stacked', sitenames='", sitesValue, "', column='", colValue, "', value='", valValue, "', cepnames=F, write.csv=F, csv.file='community.csv')", sep="")
            }else{
                command <- paste("import.from.Excel2007('", file, "', data.type='stacked', sheet='stacked', sitenames='", sitesValue, "', column='", colValue, "', value='", valValue, "', factor='", factorValue, "', level='", levelValue, "', cepnames=F, write.csv=F, csv.file='community.csv')", sep="")
            }
        }
        logger(paste(comdsnameValue, " <- ", command, sep=""))
        assign(comdsnameValue, justDoIt(command), envir=.GlobalEnv)
        communityDataSet(comdsnameValue)
        command <- paste("import.from.Excel2007('", file, "', data.type='environmental', sheet='environmental', sitenames='", sitesValue, "', write.csv=F, csv.file='environmental.csv')", sep="")
        logger(paste(envdsnameValue, " <- ", command, sep=""))
        assign(envdsnameValue, justDoIt(command), envir=.GlobalEnv)
        activeDataSet(envdsnameValue)
        tkfocus(CommanderWindow())
        }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('import.from.Excel', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(optionsFrame, text="Names for new datasets"), sticky="w")
    tkgrid(tklabel(optionsFrame, text="Enter name for community data set:"), entrycomDsname, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Enter name for environmental data set:"), entryenvDsname, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Enter name for variable with sites:"), entrysites, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Options for stacked data entry"), sticky="w")
    tkgrid(tklabel(stackedFrame, text="Import community dataset from stacked format:"), stackedCheckBox, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Enter variable for species:"), entrycol, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Enter variable for abundance:"), entryval, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Enter factor for subset:"), entryfactor, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Enter level for subset:"), entrylevel, sticky="w")
    tkgrid(stackedFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1)
}


importfromAccessGUI <- function() {
    initializeDialog(title="Read Community and Environmental data From Access")
    optionsFrame <- tkframe(top, relief="groove", borderwidth=2)
    comdsname <- tclVar("CommunityData")
    entrycomDsname <- tkentry(optionsFrame, width="20", textvariable=comdsname)
    envdsname <- tclVar("EnvironmentalDataset")
    entryenvDsname <- tkentry(optionsFrame, width="20", textvariable=envdsname)
    dsites <- tclVar("sites")
    entrysites <- tkentry(optionsFrame, width="20", textvariable=dsites)
    stackedFrame <- tkframe(top, relief="groove", borderwidth=2)
    stackedVariable <- tclVar("0")
    stackedCheckBox <- tkcheckbutton(stackedFrame, variable=stackedVariable)
    scolumn <- tclVar("species")
    entrycol <- tkentry(stackedFrame, width="20", textvariable=scolumn)
    sval <- tclVar("abundance")
    entryval <- tkentry(stackedFrame, width="20", textvariable=sval)
    sfactor <- tclVar("all")
    entryfactor <- tkentry(stackedFrame, width="20", textvariable=sfactor)
    slevel <- tclVar("all")
    entrylevel <- tkentry(stackedFrame, width="20", textvariable=slevel)
    onOK <- function(){
        closeDialog()
        comdsnameValue <- tclvalue(comdsname)
        envdsnameValue <- tclvalue(envdsname)
        sitesValue <- tclvalue(dsites)
        colValue <- tclvalue(scolumn)
        valValue <- tclvalue(sval)
        factorValue <- tclvalue(sfactor)
        levelValue <- tclvalue(slevel)
        file <- tclvalue(tkgetOpenFile(filetypes='{"Access Files" {".mdb" ".MDB"}} {"All Files" {"*"}}'))
        if (file == "") {
            if (getRcmdr("grab.focus")) tkgrab.release(top)
            tkdestroy(top)
            return()
            }
        justDoIt(paste("library(RODBC)", sep=""))
        logger(paste("library(RODBC)", sep=""))
        stacked <- tclvalue(stackedVariable) == "1"
        if (stacked==F) {
            command <- paste("import.from.Access('", file, "', data.type='community', table='community',sitenames='", sitesValue, "')", sep="")
        }else{
            if (factorValue=="all") { 
                command <- paste("import.from.Access('", file, "', data.type='stacked', table='stacked', sitenames='", sitesValue, "',column='", colValue, "',value='", valValue, "')", sep="")
            }else{
                command <- paste("import.from.Access('", file, "', data.type='stacked', table='stacked', sitenames='", sitesValue, "',column='", colValue, "',value='", valValue, "',factor='", factorValue, "',level='", levelValue, "')", sep="")
            }
        }
        logger(paste(comdsnameValue, " <- ", command, sep=""))
        assign(comdsnameValue, justDoIt(command), envir=.GlobalEnv)
        communityDataSet(comdsnameValue)
        command <- paste("import.from.Access('", file, "', data.type='environmental', table='environmental',sitenames='", sitesValue, "')", sep="")
        logger(paste(envdsnameValue, " <- ", command, sep=""))
        assign(envdsnameValue, justDoIt(command), envir=.GlobalEnv)
        activeDataSet(envdsnameValue)
        tkfocus(CommanderWindow())
        }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('import.from.Access', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(optionsFrame, text="Names for new datasets"), sticky="w")
    tkgrid(tklabel(optionsFrame, text="Enter name for community data set:"), entrycomDsname, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Enter name for environmental data set:"), entryenvDsname, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Enter name for variable with sites:"), entrysites, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Options for stacked data entry"), sticky="w")
    tkgrid(tklabel(stackedFrame, text="Import community dataset from stacked format:"), stackedCheckBox, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Enter variable for species:"), entrycol, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Enter variable for abundance:"), entryval, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Enter factor for subset:"), entryfactor, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Enter level for subset:"), entrylevel, sticky="w")
    tkgrid(stackedFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1)
}

importfromAccess2007GUI <- function() {
    initializeDialog(title="Read Community and Environmental data From Access 2007")
    optionsFrame <- tkframe(top, relief="groove", borderwidth=2)
    comdsname <- tclVar("CommunityData")
    entrycomDsname <- tkentry(optionsFrame, width="20", textvariable=comdsname)
    envdsname <- tclVar("EnvironmentalData")
    entryenvDsname <- tkentry(optionsFrame, width="20", textvariable=envdsname)
    dsites <- tclVar("sites")
    entrysites <- tkentry(optionsFrame, width="20", textvariable=dsites)
    stackedFrame <- tkframe(top, relief="groove", borderwidth=2)
    stackedVariable <- tclVar("0")
    stackedCheckBox <- tkcheckbutton(stackedFrame, variable=stackedVariable)
    scolumn <- tclVar("species")
    entrycol <- tkentry(stackedFrame, width="20", textvariable=scolumn)
    sval <- tclVar("abundance")
    entryval <- tkentry(stackedFrame, width="20", textvariable=sval)
    sfactor <- tclVar("all")
    entryfactor <- tkentry(stackedFrame, width="20", textvariable=sfactor)
    slevel <- tclVar("all")
    entrylevel <- tkentry(stackedFrame, width="20", textvariable=slevel)
    onOK <- function(){
        closeDialog()
        comdsnameValue <- tclvalue(comdsname)
        envdsnameValue <- tclvalue(envdsname)
        sitesValue <- tclvalue(dsites)
        colValue <- tclvalue(scolumn)
        valValue <- tclvalue(sval)
        factorValue <- tclvalue(sfactor)
        levelValue <- tclvalue(slevel)
        file <- tclvalue(tkgetOpenFile(filetypes='{"Access Files" {".mdbx" ".MDBX"}} {"All Files" {"*"}}'))
        if (file == "") {
            if (getRcmdr("grab.focus")) tkgrab.release(top)
            tkdestroy(top)
            return()
            }
        justDoIt(paste("library(RODBC)", sep=""))
        logger(paste("library(RODBC)", sep=""))
        stacked <- tclvalue(stackedVariable) == "1"
        if (stacked==F) {
            command <- paste("import.from.Access2007('", file, "', data.type='community', table='community',sitenames='", sitesValue, "')", sep="")
        }else{
            if (factorValue=="all") { 
                command <- paste("import.from.Access2007('", file, "', data.type='stacked', table='stacked', sitenames='", sitesValue, "',column='", colValue, "',value='", valValue, "')", sep="")
            }else{
                command <- paste("import.from.Access2007('", file, "', data.type='stacked', table='stacked', sitenames='", sitesValue, "',column='", colValue, "',value='", valValue, "',factor='", factorValue, "',level='", levelValue, "')", sep="")
            }
        }
        logger(paste(comdsnameValue, " <- ", command, sep=""))
        assign(comdsnameValue, justDoIt(command), envir=.GlobalEnv)
        communityDataSet(comdsnameValue)
        command <- paste("import.from.Access2007('", file, "', data.type='environmental', table='environmental', sitenames='", sitesValue, "')", sep="")
        logger(paste(envdsnameValue, " <- ", command, sep=""))
        assign(envdsnameValue, justDoIt(command), envir=.GlobalEnv)
        activeDataSet(envdsnameValue)
        tkfocus(CommanderWindow())
        }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('import.from.Access', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(optionsFrame, text="Names for new datasets"), sticky="w")
    tkgrid(tklabel(optionsFrame, text="Enter name for community data set:"), entrycomDsname, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Enter name for environmental data set:"), entryenvDsname, sticky="w")
    tkgrid(tklabel(optionsFrame, text="Enter name for variable with sites:"), entrysites, sticky="w")
    tkgrid(optionsFrame, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Options for stacked data entry"), sticky="w")
    tkgrid(tklabel(stackedFrame, text="Import community dataset from stacked format:"), stackedCheckBox, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Enter variable for species:"), entrycol, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Enter variable for abundance:"), entryval, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Enter factor for subset:"), entryfactor, sticky="w")
    tkgrid(tklabel(stackedFrame, text="Enter level for subset:"), entrylevel, sticky="w")
    tkgrid(stackedFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=4, columns=1)
}



samesitesGUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Same rows for community and environmental")
    .activeDataSet <- ActiveDataSet()
    .communityDataSet <- CommunityDataSet()
    saveFrame <- tkframe(top, relief="groove", borderwidth=2)
    saveVariable <- tclVar("1")
    saveCheckBox <- tkcheckbutton(saveFrame, variable=saveVariable)
    onOK <- function(){
        doItAndPrint(paste("check.datasets(", .communityDataSet, ",", .activeDataSet, ")", sep="")) 
        sav <- tclvalue(saveVariable) == "1"
        if (sav==T) {
            DataSet <- eval(parse(text=paste(.communityDataSet, sep="")), envir=.GlobalEnv)
            newname <- paste(.communityDataSet, ".orig", sep="")
            logger(paste(newname, " <- ", .communityDataSet, sep=""))
            assign(newname,DataSet, envir=.GlobalEnv)
        }
        logger(paste(.communityDataSet, " <- ", "same.sites(", .communityDataSet, "," , .activeDataSet, ")", sep=""))
        assign(.communityDataSet, justDoIt(paste("same.sites(", .communityDataSet, ",", .activeDataSet, ")", sep="")), envir=.GlobalEnv)
        communityDataSet(.communityDataSet)
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('same.sites', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(saveCheckBox, tklabel(saveFrame, text="save original community matrix"), sticky="w")
    tkgrid(saveFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkwait.window(top)
}

removezeroes <- function(){
    .communityDataSet <- CommunityDataSet()
    command <- paste("removezerospecies(", .communityDataSet, ")", sep="")
    logger(paste(.communityDataSet, " <- ", command, sep=""))
    assign(.communityDataSet, justDoIt(command), envir=.GlobalEnv)
    communityDataSet(.communityDataSet)
    invisible()
}

replaceNAs <- function(){
    .communityDataSet <- CommunityDataSet()
    command <- paste("replaceNAcomm(", .communityDataSet, ")", sep="")
    logger(paste(.communityDataSet, " <- ", command, sep=""))
    assign(.communityDataSet, justDoIt(command), envir=.GlobalEnv)
    communityDataSet(.communityDataSet)
    invisible()
}

vegemite.table <- function(){
    .communityDataSet <- CommunityDataSet()
    command <- paste("vegemite(", .communityDataSet, ", use=cca(", .communityDataSet, "), scale='Braun.Blanquet')", sep="")
    doItAndPrint(paste(command))
}

tabasco.table <- function(){
    .communityDataSet <- CommunityDataSet()
    command <- paste("tabasco(", .communityDataSet, ", use=cca(", .communityDataSet, "))", sep="")
    doItAndPrint(paste(command))
}

beals.smoothing <- function(){
    .communityDataSet <- CommunityDataSet()
    command <- paste("beals(", .communityDataSet, ", include=F)", sep="")
    doItAndPrint(paste(command))
}

ind.power <- function(){
    .communityDataSet <- CommunityDataSet()
    command <- paste("indpower(", .communityDataSet, ", type=0)", sep="")
    doItAndPrint(paste(command))
}

nested.checks <- function(){
    .communityDataSet <- CommunityDataSet()
    doItAndPrint(paste("nestedchecker(", .communityDataSet, ")", sep=""))
    doItAndPrint(paste("nestedtemp(", .communityDataSet, ")", sep=""))
    doItAndPrint(paste("nestednodf(", .communityDataSet, ")", sep=""))
    doItAndPrint(paste("nestedbetasor(", .communityDataSet, ")", sep=""))
    doItAndPrint(paste("nestedbetajac(", .communityDataSet, ")", sep=""))
}


viewcommunity <- function(){
    command <- justDoIt(paste("invisible(edit(", communityDataSet(), "))", sep=""))
}

editcommunity <- function(){
    .communityDataSet <- CommunityDataSet()
    justDoIt(paste("fix(", .communityDataSet, ")", sep=""))
    communityDataSet(.communityDataSet)
    invisible()
}

viewenvironmental <- function(){
    justDoIt(paste("invisible(edit(", ActiveDataSet(), "))", sep=""))
}

editenvironmental <- function(){
    .activeDataSet <- ActiveDataSet()
    justDoIt(paste("fix(", .activeDataSet, ")", sep=""))
    activeDataSet(.activeDataSet)
    invisible()
}

checkdatasets <- function(){
    .activeDataSet <- ActiveDataSet()
    .communityDataSet <- CommunityDataSet()
    doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
}


removeNAGUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "remove NA cases")
    .activeDataSet <- ActiveDataSet()
    .communityDataSet <- CommunityDataSet()
    .variables <- Variables()
    variables <- paste(.variables, ifelse(is.element(.variables, Factors()), "[factor]", ""))
    varFrame <- tkframe(top, relief="groove", borderwidth=2)
    subsetBox <- tklistbox(varFrame, width=27, height=7,
        selectmode="single", background="white", exportselection="FALSE") 
    subsetScroll <- tkscrollbar(varFrame, repeatinterval=5, command=function(...) tkyview(subsetBox, ...))
    tkconfigure(subsetBox, yscrollcommand=function(...) tkset(subsetScroll, ...))
    for (x in variables) tkinsert(subsetBox, "end", x)
    onOK <- function(){
        doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep="")) 
        var <- .variables[as.numeric(tkcurselection(subsetBox))+1]
        command <- paste("removeNAcomm(", .communityDataSet, ",", .activeDataSet, ",'", var, "')", sep="")
        logger(paste(.communityDataSet, " <- ", command, sep=""))
        assign(.communityDataSet, justDoIt(command), envir=.GlobalEnv)
        command <- paste("removezerospecies(", .communityDataSet, ")", sep="")
        logger(paste(.communityDataSet, " <- ", command, sep=""))
        assign(.communityDataSet, justDoIt(command), envir=.GlobalEnv)
        command <- paste("removeNAenv(", .activeDataSet, ",'", var, "')", sep="")
        logger(paste(.activeDataSet, " <- ", command, sep=""))
        assign(.activeDataSet, justDoIt(command), envir=.GlobalEnv)
        activeDataSet(.activeDataSet)
        communityDataSet(.communityDataSet)
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(varFrame, text="Select variable"), sticky="w")
    tkgrid(subsetBox, subsetScroll, sticky="w")
    tkgrid(varFrame, sticky="w")
    tkgrid(OKbutton, tklabel(buttonsFrame, text="    "), cancelButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(subsetScroll, sticky="ns")
    tkselection.set(subsetBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(subsetBox)
    tkwait.window(top)
}


disttransGUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Community matrix transformation")
    .activeDataSet <- ActiveDataSet()
    .communityDataSet <- CommunityDataSet()
    methodFrame <- tkframe(top, relief="groove", borderwidth=2)
    methodBox <- tklistbox(methodFrame, width=50, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    methodScroll <- tkscrollbar(methodFrame, repeatinterval=5, command=function(...) tkyview(methodBox, ...))
    tkconfigure(methodBox, yscrollcommand=function(...) tkset(methodScroll, ...))
    methods <- c("hellinger","chord","profiles","chi.square","log","square","pa",
        "Braun.Blanquet", "Domin", "Hult", "Hill", "fix", "coverscale.log")
    for (x in methods) tkinsert(methodBox, "end", x)
    saveFrame <- tkframe(top, relief="groove", borderwidth=2)
    saveVariable <- tclVar("1")
    saveCheckBox <- tkcheckbutton(saveFrame, variable=saveVariable)
    onOK <- function(){
        method <- methods[as.numeric(tkcurselection(methodBox))+1]
        sav <- tclvalue(saveVariable) == "1"
        if (sav==T) {
            DataSet <- eval(parse(text=paste(.communityDataSet, sep="")), envir=.GlobalEnv)
            newname <- paste(.communityDataSet, ".orig", sep="")
            logger(paste(newname, " <- ", .communityDataSet, sep=""))
            assign(newname,DataSet, envir=.GlobalEnv)
        }
        logger(paste(.communityDataSet, " <- ", "disttransform(", .communityDataSet, ", method='", method, "')", sep=""))
        assign(.communityDataSet, justDoIt(paste("disttransform(", .communityDataSet, ", method='", method, "')", sep="")), envir=.GlobalEnv)
        communityDataSet(.communityDataSet)
        }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(methodFrame, text="Method"), sticky="w")
    tkgrid(methodBox, methodScroll,sticky="w")
    tkgrid(methodFrame, sticky="w")
    tkgrid(saveCheckBox, tklabel(saveFrame, text="save original community matrix"), sticky="w")
    tkgrid(saveFrame, sticky="w")
    tkgrid(OKbutton, cancelButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(methodScroll, sticky="ns")
    tkselection.set(methodBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(methodBox)
    tkwait.window(top)
}


# BiodiversityR 2.3 reload .activeDataSet

envirosummaryGUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Summary of environmental variables")
    .activeDataSet <- ActiveDataSet()
    .variables <- Variables()
    variables <- paste(.variables, ifelse(is.element(.variables, Factors()), "[factor]", ""))
    varFrame <- tkframe(top, relief="groove", borderwidth=2)
    subsetBox <- tklistbox(varFrame, width=27, height=7,
        selectmode="single", background="white", exportselection="FALSE") 
    subsetScroll <- tkscrollbar(varFrame, repeatinterval=5, command=function(...) tkyview(subsetBox, ...))
    tkconfigure(subsetBox, yscrollcommand=function(...) tkset(subsetScroll, ...))
    variables <- c("all",variables)
    for (x in variables) tkinsert(subsetBox, "end", x)
    onOK <- function(){
        var <- variables[as.numeric(tkcurselection(subsetBox))+1]
        if (var == "all") {
            doItAndPrint(paste("summary(", .activeDataSet, ")", sep=""))
        }else{
            var <- .variables[as.numeric(tkcurselection(subsetBox))]
            doItAndPrint(paste("summary(", .activeDataSet, "$", var, ")", sep=""))
        }
    }
    onPlot <- function(){
        var <- variables[as.numeric(tkcurselection(subsetBox))+1]
        if (var == "all") {
            doItAndPrint(paste("pairs(", .activeDataSet, ")", sep=""))
        }else{
            var <- .variables[as.numeric(tkcurselection(subsetBox))]
            varfactor <- eval(parse(text=paste("is.factor(", .activeDataSet, "$", var, ")", sep="")), envir=.GlobalEnv)
            if (varfactor==T) {
                doItAndPrint(paste("plot(", .activeDataSet, "$", var,",xlab='", var, "',ylab='n')",sep=""))
            }else{
                doItAndPrint(paste("boxplot(", .activeDataSet, "$", var,",xlab='", var, "')", sep=""))
                doItAndPrint(paste("points(mean(", .activeDataSet, "$", var,"), pch=19, cex=1.5)",sep=""))
            }
        }
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    plotButton <- tkbutton(buttonsFrame, text="Plot", width="12", command=onPlot)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(varFrame, text="Select variable"), sticky="w")
    tkgrid(subsetBox, subsetScroll, sticky="w")
    tkgrid(varFrame, sticky="w")
    tkgrid(OKbutton, plotButton, tklabel(buttonsFrame, text="    "), cancelButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(subsetScroll, sticky="ns")
    tkselection.set(subsetBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(subsetBox)
    tkwait.window(top)
}


boxcoxGUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Box-Cox transformation")
    .activeDataSet <- ActiveDataSet()
    .variables <- Variables()
    variables <- paste(.variables, ifelse(is.element(.variables, Factors()), "[factor]", ""))
    varFrame <- tkframe(top, relief="groove", borderwidth=2)
    subsetBox <- tklistbox(varFrame, width=27, height=7,
        selectmode="single", background="white", exportselection="FALSE") 
    subsetScroll <- tkscrollbar(varFrame, repeatinterval=5, command=function(...) tkyview(subsetBox, ...))
    tkconfigure(subsetBox, yscrollcommand=function(...) tkset(subsetScroll, ...))
    for (x in variables) tkinsert(subsetBox, "end", x)
    onOK <- function(){
        var <- .variables[as.numeric(tkcurselection(subsetBox))+1]
        doItAndPrint(paste("par(mfrow=c(1,2))", sep=""))
        doItAndPrint(paste("qqPlot(", .activeDataSet, "$", var, ")", sep=""))
        doItAndPrint(paste("shapiro.test(", .activeDataSet, "$", var, ")", sep=""))
        doItAndPrint(paste("ks.test(", .activeDataSet, "$", var, ",pnorm)", sep=""))
        doItAndPrint(paste("summary(powerTransform(na.omit(", .activeDataSet, ")$", var, "))", sep=""))
        justDoIt(paste(.activeDataSet, "$", var, ".boxcox <- ", .activeDataSet, "$", var, "^ powerTransform(na.omit(", .activeDataSet, ")$", var, ")$lambda", sep=""))
        logger(paste(.activeDataSet, "$", var, ".boxcox <- ", .activeDataSet, "$", var, "^ powerTransform(na.omit(", .activeDataSet, ")$", var, ")$lambda", sep=""))
        activeDataSet(.activeDataSet)
        doItAndPrint(paste("qqPlot(", .activeDataSet, "$", var, ".boxcox)", sep=""))
        doItAndPrint(paste("shapiro.test(", .activeDataSet, "$", var, ".boxcox)", sep=""))
        doItAndPrint(paste("ks.test(", .activeDataSet, "$", var, ".boxcox, pnorm)", sep=""))
        doItAndPrint(paste("par(mfrow=c(1,1))", sep=""))
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(varFrame, text="Select variable"), sticky="w")
    tkgrid(subsetBox, subsetScroll, sticky="w")
    tkgrid(varFrame, sticky="w")
    tkgrid(OKbutton, tklabel(buttonsFrame, text="    "), cancelButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(subsetScroll, sticky="ns")
    tkselection.set(subsetBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(subsetBox)
    tkwait.window(top)
}


accumGUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Species accumulation curves")
    .activeDataSet <- ActiveDataSet()
    .communityDataSet <- CommunityDataSet()
    .variables <- Factors()
    variables <- paste(.variables, ifelse(is.element(.variables, Factors()), "[factor]", ""))
    .svariables <- Numeric()
    svariables <- paste(.svariables)
    .cvariables <- CVariables()
    cvariables <- paste(.cvariables)
    modelName <- tclVar("Accum.1")
    modelFrame <- tkframe(top, relief="groove", borderwidth=2)
    model <- tkentry(modelFrame, width=40, textvariable=modelName)
    choicesFrame <- tkframe(top, relief="groove", borderwidth=2)
    methodFrame <- tkframe(choicesFrame)
    method1Frame <- tkframe(methodFrame)
    method2Frame <- tkframe(methodFrame)
    methodBox <- tklistbox(method1Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    methodScroll <- tkscrollbar(method1Frame, repeatinterval=5, command=function(...) tkyview(methodBox, ...))
    tkconfigure(methodBox, yscrollcommand=function(...) tkset(methodScroll, ...))
    methods <- c("exact","exact (unconditioned)","random","rarefaction", "coleman","collector","poolaccum","estaccumR","drarefy")
    permVariable <- tclVar("100")
    permutation <- tkentry(method2Frame, width=10, textvariable=permVariable)
    for (x in methods) tkinsert(methodBox, "end", x)
    optionFrame <- tkframe(choicesFrame)
    addVariable <- tclVar("0")
    addCheckBox <- tkcheckbutton(optionFrame, variable=addVariable)
    xlist <- tclVar("")
    xEntry <- tkentry(optionFrame, width=10, textvariable=xlist)
    ylist <- tclVar("")
    yEntry <- tkentry(optionFrame, width=10, textvariable=ylist)
    symbol <- tclVar("1")
    symbolEntry <- tkentry(optionFrame, width=10, textvariable=symbol)
    cia <- tclVar("2")
    ciEntry <- tkentry(optionFrame, width=10, textvariable=cia)
    cexa <- tclVar("1")
    cexEntry <- tkentry(optionFrame, width=10, textvariable=cexa)
    colour <- tclVar("1")
    colourEntry <- tkentry(optionFrame, width=10, textvariable=colour)
    option2Frame <- tkframe(choicesFrame)   
    scaleBox <- tklistbox(option2Frame, width=27, height=6,
        selectmode="single", background="white", exportselection="FALSE") 
    scaleScroll <- tkscrollbar(option2Frame, repeatinterval=5, command=function(...) tkyview(scaleBox, ...))
    tkconfigure(scaleBox, yscrollcommand=function(...) tkset(scaleScroll, ...))
    svariables <- c("sites",svariables)
    for (x in svariables) tkinsert(scaleBox, "end", x)
    subsetFrame <- tkframe(choicesFrame)
    subset1Frame <- tkframe(subsetFrame)
    subset2Frame <- tkframe(subsetFrame)
    subsetBox <- tklistbox(subset1Frame, width=27, height=8,
        selectmode="single", background="white", exportselection="FALSE") 
    subsetScroll <- tkscrollbar(subset1Frame, repeatinterval=5, command=function(...) tkyview(subsetBox, ...))
    tkconfigure(subsetBox, yscrollcommand=function(...) tkset(subsetScroll, ...))
    variables <- c("all",variables)
    for (x in variables) tkinsert(subsetBox, "end", x)
    subset <- tclVar(".")
    subsetEntry <- tkentry(subset2Frame, width=10, textvariable=subset)
    onOK <- function(){
        doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
        modelValue <- tclvalue(modelName)
        method <- methods[as.numeric(tkcurselection(methodBox))+1]
        var <- variables[as.numeric(tkcurselection(subsetBox))+1]
        sub <- tclvalue(subset)
        xlim <- tclvalue(xlist)
        if (xlim != "") {xlim <- paste(", xlim=c(", xlim, ")", sep="")}
        ylim <- tclvalue(ylist)
        if (ylim != "") {ylim <- paste(", ylim=c(", ylim, ")", sep="")}
        perm <- as.numeric(tclvalue(permVariable))
        ci <- tclvalue(cia)
        cex <- tclvalue(cexa)
        var2 <- svariables[as.numeric(tkcurselection(scaleBox))+1]
        if (var2 == "sites") {
            scale <- paste(", scale=''", sep="")
            xlab <- paste(", xlab='sites'", sep="")
        }else{
            var2 <- .svariables[as.numeric(tkcurselection(scaleBox))]
            scale <- paste(", scale='", var2, "'", sep="")
            xlab <- paste(", xlab='", var2, "'", sep="")
        }
        if (method %in% c("exact","exact (unconditioned)","random","rarefaction", "coleman","collector")) {
            if (var == "all") {
                if (method == "exact (unconditioned)") {
                    command <- paste("accumresult(", .communityDataSet, ", y=", .activeDataSet, ", method='exact', conditioned =F, gamma = 'boot', permutations=", perm, scale, ")", sep="")
                }else{
                    command <- paste("accumresult(", .communityDataSet, ", y=", .activeDataSet, ", method='", method, "', conditioned =T, gamma = 'boot', permutations=", perm, scale, ")", sep="")
                }
            }else{
                var <- .variables[as.numeric(tkcurselection(subsetBox))]
                if (sub == ".") {
                    if (method == "exact (unconditioned)") {
                        command <- paste("accumcomp(", .communityDataSet, ", y=", .activeDataSet, ", factor='", var, "', method='exact', conditioned =F, gamma = 'boot', permutations=", perm, ", legend=F, rainbow=T, ci=", ci, ", ci.type='bar', cex=", cex, xlab, xlim, ylim, scale, ")", sep="")
                    }else{
                        command <- paste("accumcomp(", .communityDataSet, ", y=", .activeDataSet, ", factor='", var, "', method='", method, "', conditioned =T, gamma = 'boot', permutations=", perm, ", legend=F, rainbow=T, ci=", ci, ", ci.type='bar', cex=", cex, xlab, xlim, ylim, scale, ")", sep="")
                    }
                }else{
                    if (method == "exact (unconditioned)") {
                        command <- paste("accumresult(", .communityDataSet, ", y=", .activeDataSet, ", factor='", var, "', level='", sub, "', method='exact', conditioned =F, gamma = 'boot' , permutations=", perm , scale, ")", sep="")
                    }else{
                        command <- paste("accumresult(", .communityDataSet, ", y=", .activeDataSet, ", factor='", var, "', level='", sub, "', method='", method, "', conditioned =T, gamma = 'boot' , permutations=", perm , scale, ")", sep="")
                    }
                }
            }
        }
        if (method == "poolaccum") {
                command <- paste("poolaccum(", .communityDataSet, ", permutations=", perm, ", minsize=1)", sep="")
        }
        if (method == "estaccumR") {
                command <- paste("estaccumR(", .communityDataSet, ", permutations=", perm, ")", sep="")
        }
        if (method == "drarefy") {
                command <- paste("drarefy(", .communityDataSet, ", sample=round(nrow(", .communityDataSet, ")/4))", sep="")
        }
        logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        doItAndPrint(paste(modelValue))
    }
    onPlot <- function(){
        modelValue <- tclvalue(modelName)
        method <- methods[as.numeric(tkcurselection(methodBox))+1]
        addit <- tclvalue(addVariable) == "1" 
        xlim <- tclvalue(xlist)
        if (xlim != "") {xlim <- paste(", xlim=c(", xlim, ")", sep="")}
        ylim <- tclvalue(ylist)
        if (ylim != "") {ylim <- paste(", ylim=c(", ylim, ")", sep="")}
        pch <- tclvalue(symbol)
        ci <- tclvalue(cia)
        cex <- tclvalue(cexa)
        col <- tclvalue(colour)
        sub <- tclvalue(subset)
        if (sub == ".") {sub <- ""}
        var2 <- svariables[as.numeric(tkcurselection(scaleBox))+1]
        if (var2 == "sites") {
            xlab <- paste(", xlab='sites'", sep="")
        }else{
            var2 <- .variables[as.numeric(tkcurselection(scaleBox))]
            xlab <- paste(", xlab='", var2, "'", sep="")
        }
        if (method %in% c("exact","exact (unconditioned)","random","rarefaction", "coleman","collector")) {
            doItAndPrint(paste("accumplot(", modelValue, ", addit=", addit, ", ci=", ci, ", ci.type='bar', col='", col, "', cex=", cex, xlab, ", ylab='species richness'", xlim, ylim, ", pch=", pch, ", labels='", sub ,"')", sep=""))
        }
        if (method %in% c("poolaccum", "estaccumR")) {
            doItAndPrint(paste("plot(", modelValue, ")", sep=""))
        }
        activateMenus()
        tkfocus(CommanderWindow())
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('accumresult', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    plotButton <- tkbutton(buttonsFrame, text="Plot", width="12", command=onPlot)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(modelFrame, text="Save result as:    ", width=15), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(tklabel(method1Frame, text="Accumulation method"), sticky="w")
    tkgrid(methodBox, methodScroll,sticky="w")
    tkgrid(tklabel(method2Frame, text="permutations", width=10), permutation, sticky="w")
    tkgrid(method1Frame, sticky="w")
    tkgrid(method2Frame, sticky="w")
    tkgrid(tklabel(option2Frame, text="scale of x axis"), sticky="w")
    tkgrid(scaleBox, scaleScroll, sticky="w")
    tkgrid(tklabel(subsetFrame, text="Subset options"), sticky="w")
    tkgrid(subsetBox, subsetScroll, sticky="w")
    tkgrid(tklabel(subset2Frame, text="subset:                ", width=15), subsetEntry, sticky="w")
    tkgrid(subset1Frame, sticky="w")
    tkgrid(subset2Frame, sticky="w")
    tkgrid(tklabel(optionFrame, text="Plot options"), sticky="w")
    tkgrid(addCheckBox, tklabel(optionFrame, text="add plot"), sticky="w")
    tkgrid(tklabel(optionFrame, text="x limits:               ", width=15), xEntry, sticky="w")
    tkgrid(tklabel(optionFrame, text="y limits:               ", width=15), yEntry, sticky="w")
    tkgrid(tklabel(optionFrame, text="ci:                    ", width=15), ciEntry, sticky="w")
    tkgrid(tklabel(optionFrame, text="symbol:                 ", width=15), symbolEntry, sticky="w")
    tkgrid(tklabel(optionFrame, text="cex:                    ", width=15), cexEntry, sticky="w")
    tkgrid(tklabel(optionFrame, text="colour:                 ", width=15), colourEntry, sticky="w")
    tkgrid(methodFrame, tklabel(choicesFrame, text="", width=1), option2Frame, sticky="w")
    tkgrid(subsetFrame, tklabel(choicesFrame, text="", width=1), optionFrame, sticky="w")
    tkgrid(choicesFrame, sticky="w")
    tkgrid(OKbutton, plotButton, tklabel(buttonsFrame, text="    "), cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(methodScroll, sticky="ns")
    tkgrid.configure(subsetScroll, sticky="ns")
    tkgrid.configure(scaleScroll, sticky="ns")
    tkselection.set(methodBox, 0)
    tkselection.set(subsetBox, 0)
    tkselection.set(scaleBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(methodBox)
    tkwait.window(top)
}


diversityGUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Diversity calculation")
    .activeDataSet <- ActiveDataSet()
    .communityDataSet <- CommunityDataSet()
    .variables <- Factors()
    variables <- paste(.variables, ifelse(is.element(.variables, Factors()), "[factor]", ""))
    .cvariables <- CVariables()
    cvariables <- paste(.cvariables)
    modelName <- tclVar("Diversity.1")
    modelFrame <- tkframe(top, relief="groove", borderwidth=2)
    model <- tkentry(modelFrame, width=40, textvariable=modelName)
    choicesFrame <- tkframe(top, relief="groove", borderwidth=2)
    indexFrame <- tkframe(choicesFrame)
    indexBox <- tklistbox(indexFrame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    indexScroll <- tkscrollbar(indexFrame, repeatinterval=5, command=function(...) tkyview(indexBox, ...))
    tkconfigure(indexBox, yscrollcommand=function(...) tkset(indexScroll, ...))
    indices <- c("richness","abundance","Shannon","Simpson","inverseSimpson","Logalpha","Berger","Jevenness","Eevenness",
        "jack1","jack2","chao","boot","richness (contribdiv)", "simpson (contribdiv)", "eventstar")
    for (x in indices) tkinsert(indexBox, "end", x)
    methodFrame <- tkframe(choicesFrame)
    methodBox <- tklistbox(methodFrame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    methodScroll <- tkscrollbar(methodFrame, repeatinterval=5, command=function(...) tkyview(methodBox, ...))
    tkconfigure(methodBox, yscrollcommand=function(...) tkset(methodScroll, ...))
    methods <- c("all (pooled)","separate per site","mean","sd","jackknife")
    for (x in methods) tkinsert(methodBox, "end", x)
    optionFrame <- tkframe(choicesFrame)
    dataVariable <- tclVar("0")
    dataCheckBox <- tkcheckbutton(optionFrame, variable=dataVariable)
    sortVariable <- tclVar("0")
    sortCheckBox <- tkcheckbutton(optionFrame, variable=sortVariable)
    labelVariable <- tclVar("0")
    labelCheckBox <- tkcheckbutton(optionFrame, variable=labelVariable)
    addVariable <- tclVar("0")
    addCheckBox <- tkcheckbutton(optionFrame, variable=addVariable)
    ylist <- tclVar("0,5")
    yEntry <- tkentry(optionFrame, width=10, textvariable=ylist)
    symbol <- tclVar("1")
    symbolEntry <- tkentry(optionFrame, width=10, textvariable=symbol)
    subsetFrame <- tkframe(choicesFrame)
    subset1Frame <- tkframe(subsetFrame)
    subset2Frame <- tkframe(subsetFrame)
    subsetBox <- tklistbox(subset1Frame, width=27, height=7,
        selectmode="single", background="white", exportselection="FALSE") 
    subsetScroll <- tkscrollbar(subset1Frame, repeatinterval=5, command=function(...) tkyview(subsetBox, ...))
    tkconfigure(subsetBox, yscrollcommand=function(...) tkset(subsetScroll, ...))
    variables <- c("all",variables)
    for (x in variables) tkinsert(subsetBox, "end", x)
    subset <- tclVar(".")
    subsetEntry <- tkentry(subset2Frame, width=10, textvariable=subset)
    onOK <- function(){
        doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
        modelValue <- tclvalue(modelName)
        index <- indices[as.numeric(tkcurselection(indexBox))+1]
        method <- methods[as.numeric(tkcurselection(methodBox))+1]
        if (method == "all (pooled)") {method <- "all"}
        data <- tclvalue(dataVariable) == "1"
        sortit <- tclvalue(sortVariable) == "1"
        if (data==T) {sortit <- F}
        var <- variables[as.numeric(tkcurselection(subsetBox))+1]
        sub <- tclvalue(subset)
        if (index %in% c("richness","abundance","Shannon","Simpson","inverseSimpson","Logalpha","Berger","Jevenness","Eevenness",
            "jack1","jack2","chao","boot")) {
            if (var == "all") {
                command <- paste("diversityresult(", .communityDataSet, ", index='", index,
                    "' ,method='", method, "', sortit=", sortit, ", digits=3)", sep="")
            }else{
                var <- .variables[as.numeric(tkcurselection(subsetBox))]
                if (sub == ".") {
                    command <- paste("diversitycomp(", .communityDataSet, ", y=", .activeDataSet, ", factor1='", var, "', factor2='', index='", index,
                        "' ,method='", method, "', sortit=", sortit, ", digits=3)", sep="")
                }else{
                    command <- paste("diversityresult(", .communityDataSet, ", y=", .activeDataSet, ", factor='", var, "', level='", sub, "', index='", index,
                        "' ,method='", method, "', sortit=", sortit, ", digits=3)", sep="")
                }
            }
        }
        if (index == "richness (contribdiv)") {
                command <- paste("contribdiv(", .communityDataSet, ", index='richness')", sep="")
        }
        if (index == "simpson (contribdiv)") {
                command <- paste("contribdiv(", .communityDataSet, ", index='simpson')", sep="")
        }
        if (index == "eventstar") {
                command <- paste("eventstar(", .communityDataSet, ", qmax=5)", sep="")
        }
        logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        doItAndPrint(paste(modelValue))
        if (data==T && var=="all" && method=="separate per site" && index %in% c("richness","abundance","Shannon","Simpson","inverseSimpson","Logalpha","Berger","Jevenness","Eevenness","jack1","jack2","chao","boot")) {
            justDoIt(paste(.activeDataSet, "$", index, " <- diversityresult(", .communityDataSet, ", index='", index,"' ,method='", method,"')[,1]", sep=""))
            logger(paste(.activeDataSet, "$", index, " <- diversityresult(", .communityDataSet, ", index='", index,"' ,method='", method,"')[,1]", sep=""))
            activeDataSet(.activeDataSet)
        }
    }
    onPlot <- function() {
        modelValue <- tclvalue(modelName)
        index <- indices[as.numeric(tkcurselection(indexBox))+1]
        method <- methods[as.numeric(tkcurselection(methodBox))+1]
        if (method == "all (pooled)") {method <- "all"}
        var <- variables[as.numeric(tkcurselection(subsetBox))+1]
        labelit <- tclvalue(labelVariable) == "1"
        addit <- tclvalue(addVariable) == "1" 
        ylim <- tclvalue(ylist)
        if (ylim != "") {ylim <- paste(", ylim=c(", ylim, ")", sep="")}
        pch <- tclvalue(symbol)
        sub <- tclvalue(subset)
        if (index %in% c("richness","abundance","Shannon","Simpson","inverseSimpson","Logalpha","Berger","Jevenness","Eevenness")) {
            if (var!="all" && sub=="." && method!="separate per site") {
                if (addit==F) {
                    justDoIt(paste("plot(rep(-90, nrow(", modelValue, ")) ~ as.factor(rownames(", modelValue, ")), xlab='", method, "', ylab=colnames(", modelValue, "), type='n'", ylim, ")", sep=""))
                    logger(paste("plot(rep(-90, nrow(", modelValue, ")) ~ as.factor(rownames(", modelValue, ")), xlab='", method, "', ylab=colnames(", modelValue, "), type='n'", ylim, ")", sep=""))
                }
                doItAndPrint(paste("points(", modelValue, "[,2] ~ c(1:nrow(", modelValue, ")), pch=", pch, ")", sep=""))                
                if (labelit==T) {doItAndPrint(paste("text(c(1:nrow(", modelValue, "))," , modelValue, "[,2], labels=rownames(", modelValue, "), pos=3)", sep="")) }
                if (labelit==T) {doItAndPrint(paste("text(c(1:nrow(", modelValue, "))," , modelValue, "[,2], labels=", modelValue, "[,1], pos=1)", sep="")) }
            }
        }
        if (index %in% c("richness (contribdiv)", "simpson (contribdiv)")) {
            doItAndPrint(paste("plot(", modelValue, ")", sep=""))
        } 
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
    }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('diversityresult', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    plotButton <- tkbutton(buttonsFrame, text="Plot", width="12", command=onPlot)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(modelFrame, text="Save result as:    ", width=15), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(tklabel(indexFrame, text="Diversity index"), sticky="w")
    tkgrid(indexBox, indexScroll,sticky="w")
    tkgrid(tklabel(methodFrame, text="Calculation method"), sticky="w")
    tkgrid(methodBox, methodScroll,sticky="w")
    tkgrid(tklabel(subsetFrame, text="Subset options"), sticky="w")
    tkgrid(subsetBox, subsetScroll, sticky="w")
    tkgrid(tklabel(subset2Frame, text="subset:                ", width=15), subsetEntry, sticky="w")
    tkgrid(subset1Frame, sticky="w")
    tkgrid(subset2Frame, sticky="w")
    tkgrid(tklabel(optionFrame, text="Output options"), sticky="w")
    tkgrid(dataCheckBox, tklabel(optionFrame, text="save results"), sticky="w")
    tkgrid(sortCheckBox, tklabel(optionFrame, text="sort results"), sticky="w")
    tkgrid(labelCheckBox, tklabel(optionFrame, text="label results"), sticky="w")
    tkgrid(addCheckBox, tklabel(optionFrame, text="add plot"), sticky="w")
    tkgrid(tklabel(optionFrame, text="y limits:               ", width=20), yEntry, sticky="w")
    tkgrid(tklabel(optionFrame, text="symbol:                 ", width=20), symbolEntry, sticky="w")
    tkgrid(indexFrame, tklabel(choicesFrame, text="", width=1), methodFrame, sticky="w")
    tkgrid(subsetFrame, tklabel(choicesFrame, text="", width=1), optionFrame, sticky="w")
    tkgrid(choicesFrame, sticky="w")
    tkgrid(OKbutton, plotButton, tklabel(buttonsFrame, text="    "), cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(indexScroll, sticky="ns")
    tkgrid.configure(methodScroll, sticky="ns")
    tkgrid.configure(subsetScroll, sticky="ns")
    tkselection.set(methodBox, 0)
    tkselection.set(indexBox, 0)
    tkselection.set(subsetBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(methodBox)
    tkwait.window(top)
}

diversityvars <- function(){
    .communityDataSet <- CommunityDataSet()
    .activeDataSet <- ActiveDataSet()
    justDoIt(paste(.activeDataSet, " <- diversityvariables(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
    logger(paste(.activeDataSet, " <- diversityvariables(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
    activeDataSet(.activeDataSet)
}


rankabunGUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Rank abundance curves")
    .activeDataSet <- ActiveDataSet()
    .communityDataSet <- CommunityDataSet()
    .variables <- Factors()
    variables <- paste(.variables, ifelse(is.element(.variables, Factors()), "[factor]", ""))
    .cvariables <- CVariables()
    cvariables <- paste(.cvariables)
    modelName <- tclVar("RankAbun.1")
    modelFrame <- tkframe(top, relief="groove", borderwidth=2)
    model <- tkentry(modelFrame, width=40, textvariable=modelName)
    choicesFrame <- tkframe(top, relief="groove", borderwidth=2)
    optionFrame <- tkframe(choicesFrame)
    option1Frame <- tkframe(optionFrame)
    option2Frame <- tkframe(optionFrame)
    scaleBox <- tklistbox(option1Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    scaleScroll <- tkscrollbar(option1Frame, repeatinterval=5, command=function(...) tkyview(scaleBox, ...))
    tkconfigure(scaleBox, yscrollcommand=function(...) tkset(scaleScroll, ...))
    scales <- c("abundance","proportion","logabun","accumfreq")
    for (x in scales) tkinsert(scaleBox, "end", x)
    radVariable <- tclVar("0")
    radCheckBox <- tkcheckbutton(option2Frame, variable=radVariable)
    addVariable <- tclVar("0")
    addCheckBox <- tkcheckbutton(option2Frame, variable=addVariable)
    xlist <- tclVar("")
    xEntry <- tkentry(option2Frame, width=10, textvariable=xlist)
    ylist <- tclVar("")
    yEntry <- tkentry(option2Frame, width=10, textvariable=ylist)
    subsetFrame <- tkframe(choicesFrame)
    subset1Frame <- tkframe(subsetFrame)
    subset2Frame <- tkframe(subsetFrame)
    subsetBox <- tklistbox(subset1Frame, width=27, height=9,
        selectmode="single", background="white", exportselection="FALSE") 
    subsetScroll <- tkscrollbar(subset1Frame, repeatinterval=5, command=function(...) tkyview(subsetBox, ...))
    tkconfigure(subsetBox, yscrollcommand=function(...) tkset(subsetScroll, ...))
    variables <- c("all",variables)
    for (x in variables) tkinsert(subsetBox, "end", x)
    subset <- tclVar(".")
    subsetEntry <- tkentry(subset2Frame, width=10, textvariable=subset)
    onOK <- function(){
        doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
        modelValue <- tclvalue(modelName)
        var <- variables[as.numeric(tkcurselection(subsetBox))+1]
        sub <- tclvalue(subset)
        scale <- scales[as.numeric(tkcurselection(scaleBox))+1]
        xlim <- tclvalue(xlist)
        if (xlim != "") {xlim <- paste(", xlim=c(", xlim, ")", sep="")}
        ylim <- tclvalue(ylist)
        if (ylim != "") {ylim <- paste(", ylim=c(", ylim, ")", sep="")}
        if (var == "all") {
            command <- paste("rankabundance(", .communityDataSet, ")", sep="")
        }else{
            var <- .variables[as.numeric(tkcurselection(subsetBox))]
            if (sub == ".") {
                command <- paste("rankabuncomp(", .communityDataSet, ", y=", .activeDataSet, ", factor='", var, "', scale='", scale, "'", xlim, ylim, ", legend=F, rainbow=T)", sep="")
            }else{
                command <- paste("rankabundance(", .communityDataSet, ", y=", .activeDataSet, ", factor='", var, "', level='", sub, "')", sep="")
            }
        }
        logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        doItAndPrint(paste(modelValue))
    }
    onPlot <- function(){
        modelValue <- tclvalue(modelName)
        var <- variables[as.numeric(tkcurselection(subsetBox))+1]
        radfit <- tclvalue(radVariable) == "1"
        addit <- tclvalue(addVariable) == "1"
        scale <- scales[as.numeric(tkcurselection(scaleBox))+1]
        xlim <- tclvalue(xlist)
        if (xlim != "") {xlim <- paste(", xlim=c(", xlim, ")", sep="")}
        ylim <- tclvalue(ylist)
        if (ylim != "") {ylim <- paste(", ylim=c(", ylim, ")", sep="")}
        sub <- tclvalue(subset)
        if (radfit==T) {
            if (var == "all") {
                doItAndPrint(paste("radfitresult(", .communityDataSet, ")", sep=""))
            }else{
                var <- .variables[as.numeric(tkcurselection(subsetBox))]
                doItAndPrint(paste("radfitresult(", .communityDataSet, ", y=", .activeDataSet, ", factor='", var, "', level='", sub, "')", sep=""))
            }
        }else{
            if (sub == ".") {
                doItAndPrint(paste("rankabunplot(", modelValue, ",scale='", scale, "', addit=", addit, xlim, ylim, ", specnames=c(1,2,3))", sep=""))
            }else{
                doItAndPrint(paste("rankabunplot(", modelValue, ",scale='", scale, "', addit=", addit, ", labels='", sub, "'", xlim, ylim, ", specnames=c(1,2,3))", sep=""))
            }
        }
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
    }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('rankabundance', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    plotButton <- tkbutton(buttonsFrame, text="Plot", width="12", command=onPlot)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(modelFrame, text="Save result as:    ", width=15), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(tklabel(subsetFrame, text="Subset options"), sticky="w")
    tkgrid(subsetBox, subsetScroll, sticky="w")
    tkgrid(tklabel(subset2Frame, text="subset:                ", width=15), subsetEntry, sticky="w")
    tkgrid(subset1Frame, sticky="w")
    tkgrid(subset2Frame, sticky="w")
    tkgrid(tklabel(option1Frame, text="Plot options"), sticky="w")
    tkgrid(scaleBox, scaleScroll,sticky="w")
    tkgrid(radCheckBox, tklabel(option2Frame, text="fit RAD"), sticky="w")
    tkgrid(addCheckBox, tklabel(option2Frame, text="add plot"), sticky="w")
    tkgrid(tklabel(option2Frame, text="x limits:               ", width=20), xEntry, sticky="w")
    tkgrid(tklabel(option2Frame, text="y limits:               ", width=20), yEntry, sticky="w")
    tkgrid(option1Frame, sticky="w")
    tkgrid(option2Frame, sticky="w")
    tkgrid(subsetFrame, tklabel(choicesFrame, text="", width=1), optionFrame, sticky="w")
    tkgrid(choicesFrame, sticky="w")
    tkgrid(OKbutton, plotButton, tklabel(buttonsFrame, text="    "), cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(scaleScroll, sticky="ns")
    tkgrid.configure(subsetScroll, sticky="ns")
    tkselection.set(scaleBox, 0)
    tkselection.set(subsetBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkwait.window(top)
    }


renyiGUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Renyi diversity profile")
    .activeDataSet <- ActiveDataSet()
    .communityDataSet <- CommunityDataSet()
    .variables <- Factors()
    variables <- paste(.variables, ifelse(is.element(.variables, Factors()), "[factor]", ""))
    .cvariables <- CVariables()
    cvariables <- paste(.cvariables)
    modelName <- tclVar("Renyi.1")
    modelFrame <- tkframe(top, relief="groove", borderwidth=2)
    model <- tkentry(modelFrame, width=40, textvariable=modelName)
    choicesFrame <- tkframe(top, relief="groove", borderwidth=2)
    methodFrame <- tkframe(choicesFrame)
    methodBox <- tklistbox(methodFrame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    methodScroll <- tkscrollbar(methodFrame, repeatinterval=3, command=function(...) tkyview(methodBox, ...))
    tkconfigure(methodBox, yscrollcommand=function(...) tkset(methodScroll, ...))
    methods <- c("all","separate per site","accumulation")
    for (x in methods) tkinsert(methodBox, "end", x)
    scalelist <- tclVar("0,0.25,0.5,1,2,4,8,Inf")
    scaleFrame <- tkframe(choicesFrame)
    scaleEntry <- tkentry(scaleFrame, width=40, textvariable=scalelist)
    permVariable <- tclVar("100")
    permutation <- tkentry(scaleFrame, width=10, textvariable=permVariable)
    optionFrame <- tkframe(choicesFrame)
    evenVariable <- tclVar("0")
    evenCheckBox <- tkcheckbutton(optionFrame, variable=evenVariable)
    addVariable <- tclVar("0")
    addCheckBox <- tkcheckbutton(optionFrame, variable=addVariable)
    ylist <- tclVar("")
    yEntry <- tkentry(optionFrame, width=40, textvariable=ylist)
    symbol <- tclVar("1")
    symbolEntry <- tkentry(optionFrame, width=40, textvariable=symbol)
    colour <- tclVar("1")
    colourEntry <- tkentry(optionFrame, width=40, textvariable=colour)
    cexa <- tclVar("1")
    cexEntry <- tkentry(optionFrame, width=40, textvariable=cexa)
    subsetFrame <- tkframe(choicesFrame)
    subset1Frame <- tkframe(subsetFrame)
    subset2Frame <- tkframe(subsetFrame)
    subsetBox <- tklistbox(subset1Frame, width=27, height=7,
        selectmode="single", background="white", exportselection="FALSE") 
    subsetScroll <- tkscrollbar(subset1Frame, repeatinterval=5, command=function(...) tkyview(subsetBox, ...))
    tkconfigure(subsetBox, yscrollcommand=function(...) tkset(subsetScroll, ...))
    variables <- c("all",variables)
    for (x in variables) tkinsert(subsetBox, "end", x)
    subset <- tclVar(".")
    subsetEntry <- tkentry(subset2Frame, width=10, textvariable=subset)
    onOK <- function(){
        doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
        modelValue <- tclvalue(modelName)
        method <- methods[as.numeric(tkcurselection(methodBox))+1]
        scales <- tclvalue(scalelist)
        evenness <- tclvalue(evenVariable) == "1"
        ylim <- tclvalue(ylist)
        if (ylim != "") {ylim <- paste(", ylim=c(", ylim, ")", sep="")}
        var <- variables[as.numeric(tkcurselection(subsetBox))+1]
        sub <- tclvalue(subset)
        perm <- as.numeric(tclvalue(permVariable))
        if (var == "all") {
            if (method=="accumulation") {
                command <- paste("renyiaccum(", .communityDataSet, ", scales=c(", scales, "), permutations=", perm, ")", sep="")
            }else{
                command <- paste("renyiresult(", .communityDataSet, ", scales=c(", scales, "), method='", method, "')", sep="")
            }
        }else{
            var <- .variables[as.numeric(tkcurselection(subsetBox))]
            if (sub == ".") {
                command <- paste("renyicomp(", .communityDataSet, ", evenness=", evenness, ", y=", .activeDataSet, ", factor='", var, "', scales=c(", scales, "), permutations=", perm, ylim, ", legend=F)", sep="")
            }else{
                if (method=="accumulation") {
                    command <- paste("renyiaccumresult(", .communityDataSet, ", y=", .activeDataSet, ", factor='", var, "', level='", sub, "', scales=c(", scales, "), permutations=", perm, ")", sep="")
                }else{
                    command <- paste("renyiresult(", .communityDataSet, ", y=", .activeDataSet, ", factor='", var, "', level='", sub, "', scales=c(", scales, "), method='", method, "')", sep="")
                }
            }
        }
        logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        doItAndPrint(paste(modelValue))
    }
    onPlot <- function(){
        modelValue <- tclvalue(modelName)
        method <- methods[as.numeric(tkcurselection(methodBox))+1]
        evenness <- tclvalue(evenVariable) == "1"
        addit <- tclvalue(addVariable) == "1"
        ylim <- tclvalue(ylist)
        if (ylim != "") {ylim <- paste(", ylim=c(", ylim, ")", sep="")}
        pch <- tclvalue(symbol)
        col <- tclvalue(colour)
        cex <- tclvalue(cexa)
        var <- variables[as.numeric(tkcurselection(subsetBox))+1]
        sub <- tclvalue(subset)
        if (method=="accumulation") {
            justDoIt(paste("persp(", modelValue, ")", sep=""))
            logger(paste("persp(", modelValue, ")", sep=""))
        }else{
            if (var == "all" || sub != ".") {
                justDoIt(paste("renyiplot(", modelValue, ", evenness=", evenness, ", addit=", addit, ", rainbow=T, legend=F, pch=", pch, ",col='", col, "', cex=", cex, ylim, ")", sep=""))
                logger(paste("renyiplot(", modelValue, ", evenness=", evenness, ", addit=", addit, ", rainbow=T, legend=F, pch=", pch, ",col='", col, "', cex=", cex, ylim, ")", sep=""))
            }
        }
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
    }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('reyniresult', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    plotButton <- tkbutton(buttonsFrame, text="Plot", width="12", command=onPlot)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(modelFrame, text="Save result as:    ", width=15), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(tklabel(methodFrame, text="Calculation method"), sticky="w")
    tkgrid(methodBox, methodScroll,sticky="w")
    tkgrid(tklabel(scaleFrame, text=" "), sticky="w")
    tkgrid(tklabel(scaleFrame, text="scale parameters:               ", width=20), scaleEntry, sticky="w")
    tkgrid(tklabel(scaleFrame, text="permutations", width=10), permutation, sticky="w")
    tkgrid(tklabel(subsetFrame, text="Subset options"), sticky="w")
    tkgrid(subsetBox, subsetScroll, sticky="w")
    tkgrid(tklabel(subset2Frame, text="subset:                ", width=15), subsetEntry, sticky="w")
    tkgrid(subset1Frame, sticky="w")
    tkgrid(subset2Frame, sticky="w")
    tkgrid(tklabel(optionFrame, text="Plot options"), sticky="w")
    tkgrid(evenCheckBox, tklabel(optionFrame, text="evenness profile"), sticky="w")
    tkgrid(addCheckBox, tklabel(optionFrame, text="add plot"), sticky="w")
    tkgrid(tklabel(optionFrame, text="y limits:               ", width=20), yEntry, sticky="w")
    tkgrid(tklabel(optionFrame, text="symbol:                 ", width=20), symbolEntry, sticky="w")
    tkgrid(tklabel(optionFrame, text="colour:                 ", width=20), colourEntry, sticky="w")
    tkgrid(tklabel(optionFrame, text="cex:                    ", width=20), cexEntry, sticky="w")
    tkgrid(methodFrame, tklabel(choicesFrame, text="", width=1), scaleFrame, sticky="w")
    tkgrid(subsetFrame, tklabel(choicesFrame, text="", width=1), optionFrame, sticky="w")
    tkgrid(choicesFrame, sticky="w")
    tkgrid(OKbutton, plotButton, tklabel(buttonsFrame, text="    "), cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(methodScroll, sticky="ns")
    tkgrid.configure(subsetScroll, sticky="ns")
    tkselection.set(methodBox, 0)
    tkselection.set(subsetBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(methodBox)
    tkwait.window(top)
}


countGUI <- function(){
    contrasts <- c("contr.treatment", "contr.poly")
    checkAddOperator <- function(rhs){
        rhs.chars <- rev(strsplit(rhs, "")[[1]])
        if (length(rhs.chars) < 1) return(FALSE)
        check.char <- if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1)) 
                rhs.chars[1] else rhs.chars[2]
        !is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%"))
        }
    top <- tktoplevel()
    tkwm.title(top, "Analysis of species abundance")
    .activeDataSet <- ActiveDataSet()
    .communityDataSet <- CommunityDataSet()
    .variables <- Variables()
    variables <- paste(.variables, ifelse(is.element(.variables, Factors()), "[factor]", ""))
    allvars <- ""
    if (length(.variables) > 1) {
        for (i in 1:(length(.variables)-1)) {
            allvars <- paste(allvars, .variables[i], "+")
        }
        allvars <- paste(allvars, .variables[length(.variables)])
    }else{
        allvars <- paste(allvars, .variables[1])
    }
    .cvariables <- CVariables()
    cvariables <- paste(.cvariables)
    xFrame <- tkframe(top, relief="groove", borderwidth=2)
    x1Frame <- tkframe(xFrame)
    x4Frame <- tkframe(xFrame)
    x2Frame <- tkframe(x4Frame)
    x3Frame <- tkframe(x4Frame)
    xBox <- tklistbox(x2Frame, width=28, height=5,
        selectmode="single", background="white", exportselection="FALSE")
    xScroll <- tkscrollbar(x2Frame, repeatinterval=5, command=function(...) tkyview(xBox, ...))
    tkconfigure(xBox, yscrollcommand=function(...) tkset(xScroll, ...))
    for (x in variables) tkinsert(xBox, "end", x)
    resFrame <- tkframe(top, relief="groove", borderwidth=2)
    yFrame <- tkframe(resFrame)
    yBox <- tklistbox(yFrame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE")
    yScroll <- tkscrollbar(yFrame, repeatinterval=5, command=function(...) tkyview(yBox, ...), width=18)
    tkconfigure(yBox, yscrollcommand=function(...) tkset(yScroll, ...))
    for (x in cvariables) tkinsert(yBox, "end", x)
    lhsVariable <- tclVar("")
    lhsFrame <- tkframe(resFrame)
    lhsEntry <- tkentry(lhsFrame, width=28, textvariable=lhsVariable)
    rhsVariable <- tclVar("")
    rhsEntry <- tkentry(x1Frame, width=60, textvariable=rhsVariable)
    modelName <- tclVar("Count.model1")
    modelFrame <- tkframe(top, relief="groove", borderwidth=2)
    model <- tkentry(modelFrame, width=40, textvariable=modelName)
    subsetVariable <- tclVar("")
    subsetFrame <- tkframe(top, relief="groove", borderwidth=2)
    subsetEntry <- tkentry(subsetFrame, width=40, textvariable=subsetVariable)
    plotFrame <- tkframe(top, relief="groove", borderwidth=2)
    plot1Frame <- tkframe(plotFrame)
    plot2Frame <- tkframe(plotFrame)
    typeBox <- tklistbox(plot1Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    typeScroll <- tkscrollbar(plot1Frame, repeatinterval=5, command=function(...) tkyview(typeBox, ...))
    tkconfigure(typeBox, yscrollcommand=function(...) tkset(typeScroll, ...))
    types <- c("diagnostic plots","levene test (factor)","term plot","effect plot", "qq plot", "result plot (new)", 
        "result plot (add)","result plot (interpolate)", "cr plot","av plot","influence plot","multcomp (factor)","rpart")
    for (x in types) tkinsert(typeBox, "end", x)
    axisBox <- tklistbox(plot2Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    axisScroll <- tkscrollbar(plot2Frame, repeatinterval=5, command=function(...) tkyview(axisBox, ...))
    tkconfigure(axisBox, yscrollcommand=function(...) tkset(axisScroll, ...))
    for (x in variables) tkinsert(axisBox, "end", x)
    optionFrame <- tkframe(top, relief="groove", borderwidth=2)
    option1Frame <- tkframe(optionFrame)
    option2Frame <- tkframe(optionFrame)
    optionBox <- tklistbox(option1Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    optionScroll <- tkscrollbar(option1Frame, repeatinterval=5, command=function(...) tkyview(optionBox, ...))
    tkconfigure(optionBox, yscrollcommand=function(...) tkset(optionScroll, ...))
    options <- c("linear model","Poisson model","quasi-Poisson model","negative binomial model","gam model",
        "gam negbinom model","glmmPQL","rpart")
    for (x in options) tkinsert(optionBox, "end", x)    
    standardVariable <- tclVar("0")
    standardCheckBox <- tkcheckbutton(option2Frame, variable=standardVariable)
    summaryVariable <- tclVar("1")
    summaryCheckBox <- tkcheckbutton(option2Frame, variable=summaryVariable)
    anovaVariable <- tclVar("0")
    anovaCheckBox <- tkcheckbutton(option2Frame, variable=anovaVariable)
    dataVariable <- tclVar("0")
    dataCheckBox <- tkcheckbutton(option2Frame, variable=dataVariable)
    onDoubleClick <- function(){
        var <- as.character(tkget(xBox, "active"))[1]
            tkfocus(rhsEntry)
            rhs <- tclvalue(rhsVariable)
            rhs.chars <- rev(strsplit(rhs, "")[[1]])
            check.char <- if (length(rhs.chars) > 0){
                if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1)) 
                    rhs.chars[1] else rhs.chars[2]
                }
                else ""
            tclvalue(rhsVariable) <- if (rhs == "" || 
                is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%")))
                    paste(rhs, var, sep="")
                else paste(rhs, "+", var)
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onDoubleClick2 <- function(){
        var <- as.character(tkget(yBox, "active"))[1]
        lhs <- tclvalue(lhsVariable)
        tclvalue(lhsVariable) <- var
        }
    onPlus <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "+ ")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onTimes <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "*", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onColon <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, ":", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onSlash <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "/",  sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onIn <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "%in% ")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onMinus <- function(){
        rhs <- tclvalue(rhsVariable)
        tclvalue(rhsVariable) <- paste(rhs, "+I(")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onMinus2 <- function(){
        rhs <- tclvalue(rhsVariable)
        tclvalue(rhsVariable) <- paste(rhs, "s(")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onPower <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "^", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onLeftParen <- function(){
        tkfocus(rhsEntry)
        rhs <- tclvalue(rhsVariable)
        tclvalue(rhsVariable) <- paste(rhs, "(", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onRightParen <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, ")", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onOK <- function(){
        doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
        check.empty <- gsub(" ", "", tclvalue(lhsVariable))
        if ("" == check.empty) {
            tkmessageBox(message="Left-hand side of model empty.", 
                icon="error", type="ok")
            }
        check.empty <- gsub(" ", "", tclvalue(rhsVariable))
        if ("" == check.empty) {
            tkmessageBox(message="Right-hand side of model empty.", 
                icon="error", type="ok")
            }
        modelValue <- tclvalue(modelName)
        if (!is.valid.name(modelValue)){
            tkmessageBox(message=paste('"', modelValue, '" is not a valid name.', 
                sep=""), icon="error", type="ok")
            }
        right <- tclvalue(rhsVariable)
        if (right == ".") right <- allvars
        formula <- paste(tclvalue(lhsVariable), right, sep=" ~ ")
        subsetval <- tclvalue(subsetVariable)
        if (subsetval != "") {
            DataSet1 <- eval(parse(text=paste(.activeDataSet, sep="")), envir=.GlobalEnv)
            DataSet2 <- eval(parse(text=paste(.communityDataSet, sep="")), envir=.GlobalEnv)
            list <- (rownames(DataSet1) != subsetval)
            DataSet1 <- DataSet1[list,]
            DataSet2 <- DataSet2[list,]
            name1 <- paste(.activeDataSet,".m.", subsetval,sep="")
            name2 <- paste(.communityDataSet,".m.", subsetval,sep="")
            assign(name1,DataSet1, envir=.GlobalEnv)
            assign(name2,DataSet2, envir=.GlobalEnv)
            activeDataSet(name1)
            communityDataSet(name2)
        }
        stan <- tclvalue(standardVariable) == "1"
        if (stan==T) {
            DataSet1 <- eval(parse(text=paste(.activeDataSet, sep="")), envir=.GlobalEnv)
            standard <- paste(.activeDataSet, ".standard",sep="")
            for (j in 1:ncol(DataSet1)) {
                if (is.factor(DataSet1[,j]) == F) {DataSet1[,j] <- scale(DataSet1[,j])}
            }
            assign(standard,DataSet1, envir=.GlobalEnv)
            activeDataSet(standard)
        }      
        justDoIt(paste(.activeDataSet, "$", tclvalue(lhsVariable), "<- ", .communityDataSet, "$",tclvalue(lhsVariable), sep=""))
        logger(paste(.activeDataSet, "$", tclvalue(lhsVariable), "<- ", .communityDataSet, "$",tclvalue(lhsVariable), sep=""))
        option <- options[as.numeric(tkcurselection(optionBox))+1]
        if (option=="negative binomial model") {
            justDoIt(paste("library(MASS)"))
            logger(paste("library(MASS)"))
        }
        if (option=="gam model" || option=="gam negbinom model") {
            justDoIt(paste("library(mgcv)"))
            logger(paste("library(mgcv)"))
        }
        if (option=="rpart") {
            justDoIt(paste("library(rpart)"))
            logger(paste("library(rpart)"))
        }
        if (option == "linear model"){
            command <- paste("lm(", formula, ", data=",.activeDataSet, ", na.action=na.exclude)", sep="")
        }
        if (option == "Poisson model"){
            command <- paste("glm(", formula, ", family=poisson(link=log), data=",.activeDataSet, ", na.action=na.exclude)", sep="")
        }
        if (option == "quasi-Poisson model"){
            command <- paste("glm(", formula, ", family=quasipoisson(link=log), data=",.activeDataSet, ", na.action=na.exclude)", sep="")
        }
        if (option == "negative binomial model"){
            command <- paste("glm.nb(", formula, ", init.theta=1, data=",.activeDataSet, ", na.action=na.exclude)", sep="")
        }        
        if (option == "gam model"){
            command <- paste("gam(", formula, ", family=poisson(link=log), data=",.activeDataSet, ", na.action=na.exclude)", sep="")
        }
        if (option == "gam negbinom model"){
            command <- paste("gam(", formula, ", family=negbin(1), data=",.activeDataSet, ", na.action=na.exclude)", sep="")
        }
        if (option == "glmmPQL"){
            command <- paste("glmmPQL(", formula, ", family=quasipoisson(link=log), data=",.activeDataSet, ", na.action=na.exclude)", sep="")
        }
        if (option == "rpart"){
            command <- paste("rpart(", formula, ", data=",.activeDataSet, ", na.action=na.rpart, method='anova')", sep="")
        }
        logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        sum <- tclvalue(summaryVariable) == "1"
        if (sum==T) {
            doItAndPrint(paste("summary(", modelValue, ")", sep=""))
            if (option != "linear model" && option!="gam model"  && option!="gam negbinom model"  && option!="glmmPQL" && option!="rpart") {
                doItAndPrint(paste("deviancepercentage(", modelValue, ", na.omit(", .activeDataSet, "), digits=2)", sep=""))
            }
        }
        anov <- tclvalue(anovaVariable) == "1"
        if (anov==T && option!="glmmPQL" && option!="rpart") {
            doItAndPrint(paste("anova(", modelValue, ",test='F')", sep=""))
                doItAndPrint(paste("car::vif(lm(", formula, ", data=na.omit(",.activeDataSet, ")))", sep=""))
                if (option=="linear model") {
                    doItAndPrint(paste("drop1(", modelValue, ",test='F')", sep=""))
                    doItAndPrint(paste("car::Anova(", modelValue, ",type='II')", sep=""))
                }
                if (option=="Poisson model"  || option=="quasi-Poisson model"  || option=="negative binomial model") {
                    doItAndPrint(paste("drop1(", modelValue, ",test='F')", sep=""))
                    doItAndPrint(paste("car::Anova(", modelValue, ",type='II', test='F', error.estimate='deviance')", sep=""))
                    doItAndPrint(paste("car::Anova(", modelValue, ",type='II', test='Wald')", sep=""))
                }
        }        
        data <- tclvalue(dataVariable) =="1"
        if (data==T) {
            if (option=="rpart") {
                 justDoIt(paste(.activeDataSet, "$", modelValue, ".fit <- predict(", modelValue, ", type='vector', na.action=na.fail)[,2]", sep=""))
                 logger(paste(.activeDataSet, "$", modelValue, ".fit <- predict(", modelValue, ", type='vector', na.action=na.fail)[,2]", sep=""))
            }else{
                 justDoIt(paste(.activeDataSet, "$", modelValue, ".fit <- predict(", modelValue, ", type='response')", sep=""))
                 logger(paste(.activeDataSet, "$", modelValue, ".fit <- predict(", modelValue, ", type='response')", sep=""))
            }
            activeDataSet(.activeDataSet)
        }
    }
    onPlot <- function(){
        modelValue <- tclvalue(modelName)
        y <- tclvalue(lhsVariable)
        right <- tclvalue(rhsVariable)
        if (right == ".") right <- allvars
        formula <- paste(tclvalue(lhsVariable), right, sep=" ~ ")
        axisvar <- .variables[as.numeric(tkcurselection(axisBox))+1]
        varfactor <- eval(parse(text=paste("is.factor(",.activeDataSet, "$", axisvar, ")", sep="")), envir=.GlobalEnv)
        option <- options[as.numeric(tkcurselection(optionBox))+1]
        plottype <- types[as.numeric(tkcurselection(typeBox))+1]
        if (plottype == "diagnostic plots"){
            if (option=="gam model"  || option=="gam negbinom model") {
                doItAndPrint(paste("gam.check(", modelValue, ")", sep=""))
            }
            if (option=="rpart"){
                doItAndPrint(paste("plot(predict(", modelValue, ",type='vector'),residuals(", modelValue, "), xlab='predictions', ylab='residuals')", sep=""))
                doItAndPrint(paste("abline(h=0,lty=3)"))
            }
            if (option=="linear model" || option=="Poisson model"  || option=="quasi-Poisson model"  || option=="negative binomial model") {
                doItAndPrint(paste("par(mfrow=c(2,2))"))
                doItAndPrint(paste("plot(", modelValue, ")", sep=""))
                doItAndPrint(paste("par(mfrow=c(1,1))"))
            }
        }
        if (plottype == "levene test (factor)" && option !="rpart" && varfactor==T) {
            doItAndPrint(paste("leveneTest(residuals(", modelValue, "), ", .activeDataSet ,"$", axisvar, ")", sep=""))
            justDoIt(paste("plot(residuals(", modelValue, ") ~ ", .activeDataSet ,"$", axisvar, ")", sep=""))
            logger(paste("plot(residuals(", modelValue, ") ~ ", .activeDataSet ,"$", axisvar, ")", sep=""))
            doItAndPrint(paste("points(", .activeDataSet ,"$", axisvar, ",residuals(", modelValue, "))", sep=""))
        }
        if (plottype == "term plot" && option !="rpart"){
            if (option == "gam model" || option == "gam negbinom model"){
                doItAndPrint(paste("plot(", modelValue, ", se=T, rug=T)", sep=""))       
            }else{
                doItAndPrint(paste("termplot(", modelValue, ", se=T, partial.resid=T, rug=T, terms='", axisvar, "')", sep=""))
            }    
        }
        if (plottype == "qq plot" && option !="rpart") {
            doItAndPrint(paste("qqPlot(residuals(", modelValue, "))", sep=""))
            doItAndPrint(paste("shapiro.test(residuals(", modelValue, "))", sep=""))
            doItAndPrint(paste("ks.test(residuals(", modelValue, "), pnorm)", sep=""))
        }
        if (plottype == "effect plot" && option !="rpart") {
            justDoIt(paste("library(effects)", sep=""))
            logger(paste("library(effects)", sep=""))
            doItAndPrint(paste("as.data.frame(effect('", axisvar, "', ", modelValue, "))", sep=""))
            doItAndPrint(paste("plot(effect('", axisvar, "', ", modelValue, ", xlevels=500))", sep=""))
        }
        if (plottype == "result plot (new)" || plottype =="result plot (add)" || plottype == "result plot (interpolate)"){
            if (plottype == "result plot (new)"){
                justDoIt(paste("plot(", .activeDataSet, "$", y, "~ ", .activeDataSet, "$", axisvar, ", xlab='", axisvar, "', ylab='", y, "')", sep=""))
                logger(paste("plot(", .activeDataSet, "$", y, "~ ", .activeDataSet, "$", axisvar, ", xlab='", axisvar, "', ylab='", y, "')", sep=""))
            }
            if (plottype=="result plot (interpolate)" && varfactor==F) {            
                varmin <- eval(parse(text=paste("min(",.activeDataSet, "$", axisvar, ")", sep="")), envir=.GlobalEnv)
                varmax <- eval(parse(text=paste("max(",.activeDataSet, "$", axisvar, ")", sep="")), envir=.GlobalEnv)
                prdata <- paste(.activeDataSet, ".pred", sep="")
                prdatacont <- data.frame(seq(varmin,varmax,length=1000))
                colnames(prdatacont) <- axisvar
                assign(prdata, prdatacont, envir=.GlobalEnv)
            }
            if (option=="rpart") {
                if (plottype=="result plot (interpolate)" && varfactor==F) {
                    doItAndPrint(paste("points(predict(", modelValue, ", newdata=", prdata, ", type='vector') ~ ", prdata, "$", axisvar, ", type='l', lwd=2, col='red')", sep=""))
                }else{
                    doItAndPrint(paste("points(predict(", modelValue, ", newdata=", .activeDataSet, ", type='vector') ~ ", .activeDataSet, "$", axisvar, ", col='red', cex=1.5)", sep=""))
                }
            }

            if (option=="linear model" && plottype!="result plot (interpolate)") {
                prmodel <- paste(modelValue, ".pred", sep="")         
                logger(paste(prmodel, " <- data.frame(predict(", modelValue, ", newdata=", .activeDataSet, ", interval='confidence'))", sep=""))
                assign(prmodel, justDoIt(paste("data.frame(predict(", modelValue, ", newdata=", .activeDataSet, ", interval='confidence'))", sep="")), envir=.GlobalEnv)
                doItAndPrint(paste("points(", prmodel, "$fit ~ ", .activeDataSet, "$", axisvar, ", col='red', cex=1.5)", sep=""))
                doItAndPrint(paste("segments(as.numeric(", .activeDataSet, "$", axisvar, "),", prmodel, "$upr, as.numeric(", .activeDataSet, "$", axisvar, "),", prmodel, "$lwr, lty=2, col='red')", sep=""))
            }
            if (option=="linear model" && plottype=="result plot (interpolate)" && varfactor==F) {
                prmodel <- paste(modelValue, ".pred", sep="")         
                logger(paste(prmodel, " <- data.frame(predict(", modelValue, ", newdata=", prdata, ", interval='confidence'))", sep=""))
                assign(prmodel, justDoIt(paste("data.frame(predict(", modelValue, ", newdata=", prdata, ", interval='confidence'))", sep="")), envir=.GlobalEnv)
                doItAndPrint(paste("points(", prmodel, "$fit ~ ", prdata, "$", axisvar, ", type='l', lwd=2, col='red')", sep=""))
                doItAndPrint(paste("points(", prmodel, "$upr ~ ", prdata, "$", axisvar, ", type='l', lty=2, col='red')", sep=""))
                doItAndPrint(paste("points(", prmodel, "$lwr ~ ", prdata, "$", axisvar, ", type='l', lty=2, col='red')", sep=""))
            }

            if (option!="rpart" && option!="linear model" && plottype!="result plot (interpolate)") {
                prmodel <- paste(modelValue, ".pred", sep="")         
                logger(paste(prmodel, " <- predict(", modelValue, ", newdata=", .activeDataSet, ", type='response', se.fit=T)", sep=""))
                assign(prmodel, justDoIt(paste("predict(", modelValue, ", newdata=", .activeDataSet, ", type='response', se.fit=T)", sep="")), envir=.GlobalEnv)
                doItAndPrint(paste("points(", prmodel, "$fit ~ ", .activeDataSet, "$", axisvar, ", col='red', cex=1.5)", sep=""))
                doItAndPrint(paste("segments(as.numeric(", .activeDataSet, "$", axisvar, "),", prmodel, "$fit + 2*", prmodel, "$se.fit, as.numeric(", .activeDataSet, "$", axisvar, "),", prmodel, "$fit - 2*", prmodel, "$se.fit, lty=2, col='red')", sep=""))
            }
            if (option!="rpart" && option!="linear model" && plottype=="result plot (interpolate)" && varfactor==F) {
                prmodel <- paste(modelValue, ".pred", sep="")         
                logger(paste(prmodel, " <- predict(", modelValue, ", newdata=", prdata, ", type='response', se.fit=T)", sep=""))
                assign(prmodel, justDoIt(paste("predict(", modelValue, ", newdata=", prdata, ", type='response', se.fit=T)", sep="")), envir=.GlobalEnv)
                doItAndPrint(paste("points(", prmodel, "$fit ~ ", prdata, "$", axisvar, ", type='l', lwd=2, col='red')", sep=""))
                doItAndPrint(paste("points((", prmodel, "$fit + 2*", prmodel, "$se.fit) ~ ", prdata, "$", axisvar, ", type='l', lty=2, col='red')", sep=""))
                doItAndPrint(paste("points((", prmodel, "$fit - 2*", prmodel, "$se.fit) ~ ", prdata, "$", axisvar, ", type='l', lty=2, col='red')", sep=""))
            }
        }
        if (plottype == "cr plot" && option !="rpart") {
            doItAndPrint(paste("crPlots(", modelValue, ",'", axisvar, "')", sep=""))
        }
        if (plottype == "av plot" && option !="rpart") {
            doItAndPrint(paste("avPlots(", modelValue, ", ask=F, identify.points=F)", sep=""))
        }
        if (plottype == "influence plot" && option !="rpart") {
            doItAndPrint(paste("influencePlot(", modelValue, ", labels=F)", sep=""))
            doItAndPrint(paste("influence.measures(", modelValue, ")", sep=""))
        }
        if (plottype == "multcomp (factor)" && option !="rpart" && varfactor==T) {
            justDoIt(paste("library(multcomp)", sep=""))
            logger(paste("library(multcomp)", sep=""))
            doItAndPrint(paste("plot(print(confint(glht(", modelValue, ", linfct = mcp(", axisvar, "= 'Tukey')))))", sep=""))
        }
        if (plottype == "rpart" && option =="rpart") {
            justDoIt(paste("par(xpd=NA)"))
            logger(paste("par(xpd=NA)"))
            justDoIt(paste("plot(", modelValue, ", compress=T, uniform=F, branch=0.7)", sep=""))
            logger(paste("plot(", modelValue, ", compress=T, uniform=F, branch=0.7)", sep=""))
            doItAndPrint(paste("text(", modelValue, ", use.n=T, all=T, col='blue', cex=1, pretty=0, fancy=T, fwidth=0.99, fheight=0.99)", sep=""))
            justDoIt(paste("par(xpd=F)"))
            logger(paste("par(xpd=F)"))
        }

    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
    }
    .operatorFont <- getRcmdr("operatorFont")
    plusButton <- tkbutton(x3Frame, text="+", width="3", command=onPlus, 
        font=.operatorFont)
    timesButton <- tkbutton(x3Frame, text="*", width="3", command=onTimes, 
        font=.operatorFont)
    colonButton <- tkbutton(x3Frame, text=":", width="3", command=onColon, 
        font=.operatorFont)
    slashButton <- tkbutton(x3Frame, text="/", width="3", command=onSlash, 
        font=.operatorFont)
    inButton <- tkbutton(xFrame, text="%in%", width="3", command=onIn,
        font=.operatorFont)
    minusButton <- tkbutton(x3Frame, text="I(", width="3", command=onMinus, 
        font=.operatorFont)
    minus2Button <- tkbutton(x3Frame, text="s(", width="3", command=onMinus2, 
        font=.operatorFont)
    powerButton <- tkbutton(x3Frame, text="^", width="3", command=onPower, 
        font=.operatorFont)
    leftParenButton <- tkbutton(x3Frame, text="(", width="3", command=onLeftParen, 
        font=.operatorFont)
    rightParenButton <- tkbutton(x3Frame, text=")", width="3", command=onRightParen, 
        font=.operatorFont)
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    plotButton <- tkbutton(buttonsFrame, text="Plot", width="12", command=onPlot)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        help(generalizedLinearModel)
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(modelFrame, text="Save model as:               ", width=20), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(tklabel(option1Frame, text="Model options"), sticky="w")
    tkgrid(optionBox, optionScroll,sticky="w")
    tkgrid(standardCheckBox, tklabel(option2Frame, text="standardise"), sticky="w")
    tkgrid(summaryCheckBox, tklabel(option2Frame, text="print summary"), sticky="w")
    tkgrid(anovaCheckBox, tklabel(option2Frame, text="print anova"), sticky="w")
    tkgrid(dataCheckBox, tklabel(option2Frame, text="add predictions to dataframe"), sticky="w")
    tkgrid(option1Frame, tklabel(optionFrame, text="", width=1), option2Frame, sticky="w")
    tkgrid(optionFrame, sticky="w")
    tkgrid(tklabel(lhsFrame, text="Response"), sticky="w")
    tkgrid(lhsEntry, sticky="nw")
    tkgrid(yBox, yScroll, sticky="nw")
    tkgrid(lhsFrame,tklabel(resFrame, text="", width=1), yFrame)
    tkgrid(resFrame, sticky="w")
    tkgrid(rhsEntry, sticky="w")
    tkgrid(xBox, xScroll,sticky="w") 
    tkgrid(plusButton, timesButton, colonButton, slashButton, inButton, sticky="w")
    tkgrid(minusButton,powerButton, leftParenButton, rightParenButton, minus2Button, sticky="w")
    tkgrid(tklabel(xFrame, text="Explanatory"), sticky="w")
    tkgrid(x1Frame, sticky="w")
    tkgrid(x2Frame, tklabel(xFrame, text="", width=1), x3Frame, sticky="w")
    tkgrid(x4Frame, sticky="w")
    tkgrid(xFrame, sticky="w")
    tkgrid(tklabel(subsetFrame, text="Remove site with name: ", width=20), subsetEntry, sticky="w")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(tklabel(plot1Frame, text="Plot options"), sticky="w")
    tkgrid(typeBox, typeScroll, sticky="nw")
    tkgrid(tklabel(plot2Frame, text="Plot variable"), sticky="w")
    tkgrid(axisBox, axisScroll, sticky="nw")
    tkgrid(plot1Frame, tklabel(plotFrame, text="", width=1), plot2Frame, sticky="w")
    tkgrid(plotFrame, sticky="w")
    tkgrid(OKbutton, plotButton, cancelButton, tklabel(buttonsFrame, text="            "), 
        helpButton, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(xScroll, sticky="ns")
    tkgrid.configure(yScroll, sticky="ns")
    tkgrid.configure(typeScroll, sticky="ns")
    tkgrid.configure(axisScroll, sticky="ns")
    tkgrid.configure(optionScroll, sticky="ns")
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkselection.set(typeBox, 0)
    tkselection.set(optionBox, 0)
    tkselection.set(axisBox, 0)
    tkbind(top, "<Return>", onOK)
    tkbind(xBox, "<Double-ButtonPress-1>", onDoubleClick)
    tkbind(yBox, "<Double-ButtonPress-1>", onDoubleClick2)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(lhsEntry)
    tkwait.window(top)
}    


presabsGUI <- function(){
    contrasts <- c("contr.treatment", "contr.poly")
    checkAddOperator <- function(rhs){
        rhs.chars <- rev(strsplit(rhs, "")[[1]])
        if (length(rhs.chars) < 1) return(FALSE)
        check.char <- if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1)) 
                rhs.chars[1] else rhs.chars[2]
        !is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%"))
        }
    top <- tktoplevel()
    tkwm.title(top, "Analysis of presence/absence")
    .activeDataSet <- ActiveDataSet()
    .communityDataSet <- CommunityDataSet()
    .variables <- Variables()
    variables <- paste(.variables, ifelse(is.element(.variables, Factors()), "[factor]", ""))
    allvars <- ""
    if (length(.variables) > 1) {
        for (i in 1:(length(.variables)-1)) {
            allvars <- paste(allvars, .variables[i], "+")
        }
        allvars <- paste(allvars, .variables[length(.variables)])
    }else{
        allvars <- paste(allvars, .variables[1])
    }
    .cvariables <- CVariables()
    cvariables <- paste(.cvariables)
    xFrame <- tkframe(top, relief="groove", borderwidth=2)
    x1Frame <- tkframe(xFrame)
    x4Frame <- tkframe(xFrame)
    x2Frame <- tkframe(x4Frame)
    x3Frame <- tkframe(x4Frame)
    xBox <- tklistbox(x2Frame, width=28, height=5,
        selectmode="single", background="white", exportselection="FALSE")
    xScroll <- tkscrollbar(x2Frame, repeatinterval=5, command=function(...) tkyview(xBox, ...))
    tkconfigure(xBox, yscrollcommand=function(...) tkset(xScroll, ...))
    for (x in variables) tkinsert(xBox, "end", x)
    resFrame <- tkframe(top, relief="groove", borderwidth=2)
    yFrame <- tkframe(resFrame)
    yBox <- tklistbox(yFrame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE")
    yScroll <- tkscrollbar(yFrame, repeatinterval=5, command=function(...) tkyview(yBox, ...), width=18)
    tkconfigure(yBox, yscrollcommand=function(...) tkset(yScroll, ...))
    for (x in cvariables) tkinsert(yBox, "end", x)
    lhsVariable <- tclVar("")
    lhsFrame <- tkframe(resFrame)
    lhsEntry <- tkentry(lhsFrame, width=28, textvariable=lhsVariable)
    rhsVariable <- tclVar("")
    rhsEntry <- tkentry(x1Frame, width=60, textvariable=rhsVariable)
    modelName <- tclVar("Presabs.model1")
    modelFrame <- tkframe(top, relief="groove", borderwidth=2)
    model <- tkentry(modelFrame, width=40, textvariable=modelName)
    subsetVariable <- tclVar("")
    subsetFrame <- tkframe(top, relief="groove", borderwidth=2)
    subsetEntry <- tkentry(subsetFrame, width=40, textvariable=subsetVariable)
    plotFrame <- tkframe(top, relief="groove", borderwidth=2)
    plot1Frame <- tkframe(plotFrame)
    plot2Frame <- tkframe(plotFrame)
    typeBox <- tklistbox(plot1Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    typeScroll <- tkscrollbar(plot1Frame, repeatinterval=5, command=function(...) tkyview(typeBox, ...))
    tkconfigure(typeBox, yscrollcommand=function(...) tkset(typeScroll, ...))
    types <- c("tabular","diagnostic plots","levene test (factor)","term plot","effect plot","qq plot","result plot (new)", 
        "result plot (add)","result plot (interpolate)", "cr plot","av plot","influence plot","multcomp (factor)","rpart")
    for (x in types) tkinsert(typeBox, "end", x)
    axisBox <- tklistbox(plot2Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    axisScroll <- tkscrollbar(plot2Frame, repeatinterval=5, command=function(...) tkyview(axisBox, ...))
    tkconfigure(axisBox, yscrollcommand=function(...) tkset(axisScroll, ...))
    for (x in variables) tkinsert(axisBox, "end", x)
    optionFrame <- tkframe(top, relief="groove", borderwidth=2)
    option1Frame <- tkframe(optionFrame)
    option2Frame <- tkframe(optionFrame)
    optionBox <- tklistbox(option1Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    optionScroll <- tkscrollbar(option1Frame, repeatinterval=5, command=function(...) tkyview(optionBox, ...))
    tkconfigure(optionBox, yscrollcommand=function(...) tkset(optionScroll, ...))
    options <- c("crosstab","binomial model","quasi-binomial model","gam model", "gam quasi-binomial model","rpart","nnetrandom")
    for (x in options) tkinsert(optionBox, "end", x)    
    standardVariable <- tclVar("0")
    standardCheckBox <- tkcheckbutton(option2Frame, variable=standardVariable)
    summaryVariable <- tclVar("1")
    summaryCheckBox <- tkcheckbutton(option2Frame, variable=summaryVariable)
    anovaVariable <- tclVar("0")
    anovaCheckBox <- tkcheckbutton(option2Frame, variable=anovaVariable)
    dataVariable <- tclVar("0")
    dataCheckBox <- tkcheckbutton(option2Frame, variable=dataVariable)
    onDoubleClick <- function(){
        var <- as.character(tkget(xBox, "active"))[1]
            tkfocus(rhsEntry)
            rhs <- tclvalue(rhsVariable)
            rhs.chars <- rev(strsplit(rhs, "")[[1]])
            check.char <- if (length(rhs.chars) > 0){
                if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1)) 
                    rhs.chars[1] else rhs.chars[2]
                }
                else ""
            tclvalue(rhsVariable) <- if (rhs == "" || 
                is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%")))
                    paste(rhs, var, sep="")
                else paste(rhs, "+", var)
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onDoubleClick2 <- function(){
        var <- as.character(tkget(yBox, "active"))[1]
        lhs <- tclvalue(lhsVariable)
        tclvalue(lhsVariable) <- var
        }
    onPlus <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "+ ")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onTimes <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "*", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onColon <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, ":", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onSlash <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "/",  sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onIn <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "%in% ")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onMinus <- function(){
        rhs <- tclvalue(rhsVariable)
        tclvalue(rhsVariable) <- paste(rhs, "+I(")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onMinus2 <- function(){
        rhs <- tclvalue(rhsVariable)
        tclvalue(rhsVariable) <- paste(rhs, "s(")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onPower <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "^", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onLeftParen <- function(){
        tkfocus(rhsEntry)
        rhs <- tclvalue(rhsVariable)
        tclvalue(rhsVariable) <- paste(rhs, "(", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onRightParen <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, ")", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onOK <- function(){
        doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
        check.empty <- gsub(" ", "", tclvalue(lhsVariable))
        if ("" == check.empty) {
            tkmessageBox(message="Left-hand side of model empty.", 
                icon="error", type="ok")
            tkgrab.release(top)
            tkdestroy(top)
            generalizedLinearModel()
            return()
            }
        check.empty <- gsub(" ", "", tclvalue(rhsVariable))
        if ("" == check.empty) {
            tkmessageBox(message="Right-hand side of model empty.", 
                icon="error", type="ok")
            }
        modelValue <- tclvalue(modelName)
        if (!is.valid.name(modelValue)){
            tkmessageBox(message=paste('"', modelValue, '" is not a valid name.', 
                sep=""), icon="error", type="ok")
            }
        y <- paste(tclvalue(lhsVariable), ">0", sep="")
        right <- tclvalue(rhsVariable)
        if (right == ".") right <- allvars
        formula <- paste(y, right, sep=" ~ ")
        subsetval <- tclvalue(subsetVariable)
        if (subsetval != "") {
            DataSet1 <- eval(parse(text=paste(.activeDataSet, sep="")), envir=.GlobalEnv)
            DataSet2 <- eval(parse(text=paste(.communityDataSet, sep="")), envir=.GlobalEnv)
            list <- (rownames(DataSet1) != subsetval)
            DataSet1 <- DataSet1[list,]
            DataSet2 <- DataSet2[list,]
            name1 <- paste(.activeDataSet,".m.", subsetval,sep="")
            name2 <- paste(.communityDataSet,".m.", subsetval,sep="")
            assign(name1,DataSet1, envir=.GlobalEnv)
            assign(name2,DataSet2, envir=.GlobalEnv)
            activeDataSet(name1)
            communityDataSet(name2)
        }
        stan <- tclvalue(standardVariable) == "1"
        if (stan==T) {
            DataSet1 <- eval(parse(text=paste(.activeDataSet, sep="")), envir=.GlobalEnv)
            standard <- paste(.activeDataSet, ".standard",sep="")
            for (j in 1:ncol(DataSet1)) {
                if (is.factor(DataSet1[,j]) == F) {DataSet1[,j] <- scale(DataSet1[,j])}
            }
            assign(standard,DataSet1, envir=.GlobalEnv)
            activeDataSet(standard)
        }
        justDoIt(paste(.activeDataSet, "$", tclvalue(lhsVariable), "<- ", .communityDataSet, "$",tclvalue(lhsVariable), sep=""))
        logger(paste(.activeDataSet, "$", tclvalue(lhsVariable), "<- ", .communityDataSet, "$",tclvalue(lhsVariable), sep=""))
        option <- options[as.numeric(tkcurselection(optionBox))+1]
        if (option=="gam model" || option=="gam quasi-binomial model") {
            justDoIt(paste("library(mgcv)"))
            logger(paste("library(mgcv)"))
        }
        if (option=="rpart") {
            justDoIt(paste("library(rpart)"))
            logger(paste("library(rpart)"))
        }
        if (option=="nnetrandom") {
            justDoIt(paste("library(nnet)"))
            logger(paste("library(nnet)"))
            justDoIt(paste(.activeDataSet, "$presence <- as.numeric(", .communityDataSet, "$",tclvalue(lhsVariable), ">0)", sep=""))
            logger(paste(.activeDataSet, "$presence <- as.numeric(", .communityDataSet, "$",tclvalue(lhsVariable), ">0)", sep=""))
            justDoIt(paste("attach(", .activeDataSet, ", pos=2)",sep=""))
            logger(paste("attach(", .activeDataSet, ", pos=2)",sep=""))
            activeDataSet(.activeDataSet)
            formula <- paste("presence", right, sep=" ~ ")            
        }
        if (option == "crosstab"){
            y <- tclvalue(lhsVariable)
            command <- paste("crosstabanalysis(na.omit(", .activeDataSet, "),'", y, "','", right, "')", sep="")
        }
        if (option == "binomial model"){
            command <- paste("glm(", formula, ", family=binomial(link=logit), data=",.activeDataSet, ", na.action=na.exclude)", sep="")
        }
        if (option == "quasi-binomial model"){
            command <- paste("glm(", formula, ", family=quasibinomial(link=logit), data=",.activeDataSet, ", na.action=na.exclude)", sep="")
        }
        if (option == "gam model"){
            command <- paste("gam(", formula, ", family=binomial(link=logit), data=",.activeDataSet, ", na.action=na.exclude)", sep="")
        }
        if (option == "gam quasi-binomial model"){
            command <- paste("gam(", formula, ", family=quasibinomial(link=logit), data=",.activeDataSet, ", na.action=na.exclude)", sep="")
        }
        if (option == "rpart"){
            command <- paste("rpart(", formula, ", data=",.activeDataSet, ", method='class', na.action=na.rpart)", sep="")
        }
        if (option == "nnetrandom"){
            command <- paste("nnetrandom(", formula, ", data=",.activeDataSet, ", size=2, skip=T, entropy=T, trace=F, maxit=1000, tries=500, leave.one.out=F)", sep="")
        }
        logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        sum <- tclvalue(summaryVariable) == "1"
        if (sum==T && option !="crosstab") {
            doItAndPrint(paste("summary(", modelValue, ")", sep=""))
            if (option=="binomial model") {
                doItAndPrint(paste("deviancepercentage(", modelValue, ", na.omit(", .activeDataSet, "), test='Chi', digits=2)", sep=""))
            }
            if (option=="quasi-binomial model") {
                doItAndPrint(paste("deviancepercentage(", modelValue, ", na.omit(", .activeDataSet, "), test='F', digits=2)", sep=""))
            }
        }
        if (sum==T && option=="crosstab") {        
            doItAndPrint(paste(modelValue))
            doItAndPrint(paste(modelValue, "$observed", sep=""))
            doItAndPrint(paste(modelValue, "$expected", sep=""))
        }
        anov <- tclvalue(anovaVariable) == "1"
        if (anov==T && (option=="binomial model" || option=="gam model")) {
            doItAndPrint(paste("anova(", modelValue, ",test='Chi')", sep=""))
                doItAndPrint(paste("car::vif(lm(", formula, ", data=na.omit(",.activeDataSet, ")))", sep=""))                  
                doItAndPrint(paste("drop1(", modelValue, ",test='Chi')", sep=""))
                doItAndPrint(paste("car::Anova(", modelValue, ",type='II', test='F', error.estimate='deviance')", sep=""))
                doItAndPrint(paste("car::Anova(", modelValue, ",type='II', test='Wald')", sep=""))
        } 
        if (anov==T && (option=="quasi-binomial model" || option=="gam quasi-binomial model")) {
            doItAndPrint(paste("anova(", modelValue, ",test='F')", sep=""))
                doItAndPrint(paste("car::vif(lm(", formula, ", data=na.omit(",.activeDataSet, ")))", sep=""))                  
                doItAndPrint(paste("drop1(", modelValue, ",test='F')", sep=""))
                doItAndPrint(paste("car::Anova(", modelValue, ",type='II', test='F', error.estimate='deviance')", sep=""))
                doItAndPrint(paste("car::Anova(", modelValue, ",type='II', test='Wald')", sep=""))
        }       
        data <- tclvalue(dataVariable) =="1"
        if (data==T) {
            if (option=="rpart") {
                 justDoIt(paste(.activeDataSet, "$", modelValue, ".fit <- predict(", modelValue, ", type='prob', na.action=na.fail)[,2]", sep=""))
                 logger(paste(.activeDataSet, "$", modelValue, ".fit <- predict(", modelValue, ", type='prob', na.action=na.fail)[,2]", sep=""))
            }
            if (option=="nnetrandom") {
                 justDoIt(paste(.activeDataSet, "$", modelValue, ".fit <- predict(", modelValue, ", newdata=", .activeDataSet, ", type='raw', na.action=na.fail)", sep=""))
                 logger(paste(.activeDataSet, "$", modelValue, ".fit <- predict(", modelValue, ", newdata=", .activeDataSet, ", type='raw', na.action=na.fail)", sep=""))
            }
            if (option!="rpart" && option!="nnetrandom" && option!="crosstab") {
                 justDoIt(paste(.activeDataSet, "$", modelValue, ".fit <- predict(", modelValue, ", type='response')", sep=""))
                 logger(paste(.activeDataSet, "$", modelValue, ".fit <- predict(", modelValue, ", type='response')", sep=""))
            }
            activeDataSet(.activeDataSet)
        }
    }
    onPlot <- function(){
        modelValue <- tclvalue(modelName)
        axisvar <- .variables[as.numeric(tkcurselection(axisBox))+1]
        varfactor <- eval(parse(text=paste("is.factor(",.activeDataSet, "$", axisvar, ")", sep="")), envir=.GlobalEnv)
        option <- options[as.numeric(tkcurselection(optionBox))+1]
        plottype <- types[as.numeric(tkcurselection(typeBox))+1]
        y <- paste(tclvalue(lhsVariable), ">0", sep="")
        if (plottype == "tabular") {
                doItAndPrint(paste("plot(as.factor(", y, ") ~", axisvar, ", na.omit(", .activeDataSet, "))", sep=""))
        }
        if (plottype == "diagnostic plots"){
            if (option=="gam model"  || option=="gam quasi-binomial model") {
                doItAndPrint(paste("gam.check(", modelValue, ")", sep=""))
            }
            if (option=="binomial model" || option=="quasi-binomial model"){
                doItAndPrint(paste("par(mfrow=c(2,2))"))
                doItAndPrint(paste("plot(", modelValue, ")", sep=""))
                doItAndPrint(paste("par(mfrow=c(1,1))"))
            }
            if (option=="rpart" || option=="nnetrandom"){
                doItAndPrint(paste("plot(as.factor(predict(", modelValue, ",newdata=na.omit(", .activeDataSet, "), type='class')) ~ as.factor(na.omit(", .activeDataSet, ")$", y, "), xlab='observed',ylab='predicted')", sep=""))
            }
        }
        if (plottype == "levene test (factor)" && option !="crosstab" && option !="rpart" && option !="nnetrandom" && varfactor==T) {
            doItAndPrint(paste("leveneTest(residuals(", modelValue, "), ", .activeDataSet ,"$", axisvar, ")", sep=""))
            justDoIt(paste("plot(residuals(", modelValue, ") ~ ", .activeDataSet ,"$", axisvar, ")", sep=""))
            logger(paste("plot(residuals(", modelValue, ") ~ ", .activeDataSet ,"$", axisvar, ")", sep=""))
            doItAndPrint(paste("points(", .activeDataSet ,"$", axisvar, ",residuals(", modelValue, "))", sep=""))
        }
        if (plottype == "term plot" && option !="crosstab" && option !="rpart" && option !="nnetrandom"){
            if (option == "gam model" || option=="gam quasi-binomial model"){
                doItAndPrint(paste("plot(", modelValue, ", se=T, rug=T)", sep=""))       
            }else{
                doItAndPrint(paste("termplot(", modelValue, ", se=T, partial.resid=T, rug=T, terms='", axisvar, "')", sep=""))
            }    
        }
        if (plottype == "effect plot" && option !="crosstab" && option !="rpart" && option !="nnetrandom") {
            justDoIt(paste("library(effects)", sep=""))
            logger(paste("library(effects)", sep=""))
            doItAndPrint(paste("as.data.frame(effect('", axisvar, "', ", modelValue, "))", sep=""))
            doItAndPrint(paste("plot(effect('", axisvar, "', ", modelValue, ", xlevels=500))", sep=""))
        }
        if (plottype == "qq plot" && option !="crosstab" && option !="rpart" && option !="nnetrandom") {
            doItAndPrint(paste("qqPlot(residuals(", modelValue, "))", sep=""))
            doItAndPrint(paste("shapiro.test(residuals(", modelValue, "))", sep=""))
            doItAndPrint(paste("ks.test(residuals(", modelValue, "), pnorm)", sep=""))
        }
        if (plottype == "result plot (new)" || plottype =="result plot (add)" || plottype == "result plot (interpolate)"){
            if (plottype == "result plot (new)"){
                if (varfactor==T){ 
                    justDoIt(paste("plot(rep(-9, nrow(", .activeDataSet, ")) ~ ", .activeDataSet, "$", axisvar, ", xlab='", axisvar, "', ylab='", tclvalue(lhsVariable), " (presence-absence)', type='n', ylim=c(0,1))", sep=""))
                    logger(paste("plot(rep(-9, nrow(", .activeDataSet, ")) ~ ", .activeDataSet, "$", axisvar, ", xlab='", axisvar, "', ylab='", tclvalue(lhsVariable), " (presence-absence)', type='n', ylim=c(0,1))", sep=""))
                }else{
                    justDoIt(paste("plot(", .activeDataSet, "$", y, "~ ", .activeDataSet, "$", axisvar, ", xlab='", axisvar, "', ylab='", tclvalue(lhsVariable), " (presence-absence)', ylim=c(0,1))", sep=""))
                    logger(paste("plot(", .activeDataSet, "$", y, "~ ", .activeDataSet, "$", axisvar, ", xlab='", axisvar, "', ylab='", tclvalue(lhsVariable), " (presence-absence)', ylim=c(0,1))", sep=""))
                }
                doItAndPrint(paste("abline(h=0,lty=3)"))
                doItAndPrint(paste("abline(h=0.5,lty=3)"))
                doItAndPrint(paste("abline(h=1,lty=3)"))
            }
            if (plottype=="result plot (interpolate)" && varfactor==F) {            
                varmin <- eval(parse(text=paste("min(",.activeDataSet, "$", axisvar, ")", sep="")), envir=.GlobalEnv)
                varmax <- eval(parse(text=paste("max(",.activeDataSet, "$", axisvar, ")", sep="")), envir=.GlobalEnv)
                prdata <- paste(.activeDataSet, ".pred", sep="")
                prdatacont <- data.frame(seq(varmin,varmax,length=1000))
                colnames(prdatacont) <- axisvar
                assign(prdata, prdatacont, envir=.GlobalEnv)
            }
            if (option=="rpart") {
                if (plottype=="result plot (interpolate)" && varfactor==F) {
                    doItAndPrint(paste("points(predict(", modelValue, ", newdata=", prdata, ", type='prob')[,2] ~ ", prdata, "$", axisvar, ", type='l', lwd=2, col='red')", sep=""))
                }else{
                    doItAndPrint(paste("points(predict(", modelValue, ", newdata=", .activeDataSet, ", type='prob')[,2] ~ ", .activeDataSet, "$", axisvar, ", col='red', cex=1.5)", sep=""))
                }
            }
            if (option=="nnetrandom") {
                if (plottype=="result plot (interpolate)" && varfactor==F) {
                    doItAndPrint(paste("points(predict(", modelValue, ",newdata=", prdata, ", type='raw') ~ ", prdata, "$", axisvar, ", col='red', type='l', lwd=2)", sep=""))
                }else{
                    doItAndPrint(paste("points(predict(", modelValue, ",newdata=", .activeDataSet, ", type='raw') ~ ", .activeDataSet, "$", axisvar, ", col='red', cex=1.5)", sep=""))
                }
            }
            if (option!="nnetrandom" && option!="rpart" && option!="crosstab" && plottype!="result plot (interpolate)") {
                prmodel <- paste(modelValue, ".pred", sep="")         
                logger(paste(prmodel, " <- predict(", modelValue, ", newdata=", .activeDataSet, ", type='response', se.fit=T)", sep=""))
                assign(prmodel, justDoIt(paste("predict(", modelValue, ", newdata=", .activeDataSet, ", type='response', se.fit=T)", sep="")), envir=.GlobalEnv)
                doItAndPrint(paste("points(", prmodel, "$fit ~ ", .activeDataSet, "$", axisvar, ", col='red', cex=1.5)", sep=""))
                doItAndPrint(paste("segments(as.numeric(", .activeDataSet, "$", axisvar, "),", prmodel, "$fit + 2*", prmodel, "$se.fit, as.numeric(", .activeDataSet, "$", axisvar, "),", prmodel, "$fit - 2*", prmodel, "$se.fit, lty=2, col='red')", sep=""))
            }
            if (option!="nnetrandom" && option!="rpart" && option!="crosstab" && plottype=="result plot (interpolate)" && varfactor==F) {
                prmodel <- paste(modelValue, ".pred", sep="")         
                logger(paste(prmodel, " <- predict(", modelValue, ", newdata=", prdata, ", type='response', se.fit=T)", sep=""))
                assign(prmodel, justDoIt(paste("predict(", modelValue, ", newdata=", prdata, ", type='response', se.fit=T)", sep="")), envir=.GlobalEnv)
                doItAndPrint(paste("points(", prmodel, "$fit ~ ", prdata, "$", axisvar, ", type='l', lwd=2, col='red')", sep=""))
                doItAndPrint(paste("points((", prmodel, "$fit + 2*", prmodel, "$se.fit) ~ ", prdata, "$", axisvar, ", type='l', lty=2, col='red')", sep=""))
                doItAndPrint(paste("points((", prmodel, "$fit - 2*", prmodel, "$se.fit) ~ ", prdata, "$", axisvar, ", type='l', lty=2, col='red')", sep=""))
            }
        }
        if (plottype == "cr plot" && option !="crosstab" && option !="rpart" && option !="nnetrandom") {
            doItAndPrint(paste("crPlots(", modelValue, ",'", axisvar, "')", sep=""))
        }
        if (plottype == "av plot" && option !="crosstab" && option !="rpart" && option !="nnetrandom") {
            doItAndPrint(paste("avPlots(", modelValue, ", ask=F, identify.points=F)", sep=""))
        }
        if (plottype == "influence plot" && option !="crosstab" && option !="rpart" && option !="nnetrandom") {
            doItAndPrint(paste("influencePlot(", modelValue, ", labels=F)", sep=""))
            doItAndPrint(paste("influence.measures(", modelValue, ")", sep=""))
        }
        if (plottype == "multcomp (factor)" && option !="crosstab" && option !="rpart" && option !="nnetrandom" && varfactor==T) {
            justDoIt(paste("library(multcomp)", sep=""))
            logger(paste("library(multcomp)", sep=""))
            doItAndPrint(paste("plot(print(confint(glht(", modelValue, ", linfct = mcp(", axisvar, "= 'Tukey')))))", sep=""))
        }
        if (plottype == "rpart" && option=="rpart") {
            justDoIt(paste("par(xpd=NA)"))
            logger(paste("par(xpd=NA)"))
            justDoIt(paste("plot(", modelValue, ", compress=T, uniform=F, branch=0.7)", sep=""))
            logger(paste("plot(", modelValue, ", compress=T, uniform=F, branch=0.7)", sep=""))
            doItAndPrint(paste("text(", modelValue, ", use.n=T, all=T, col='blue', cex=1, pretty=0, fancy=T, fwidth=0.99, fheight=0.99)", sep=""))
            justDoIt(paste("par(xpd=F)"))
            logger(paste("par(xpd=F)"))
        }
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
    }
    .operatorFont <- getRcmdr("operatorFont")
    plusButton <- tkbutton(x3Frame, text="+", width="3", command=onPlus, 
        font=.operatorFont)
    timesButton <- tkbutton(x3Frame, text="*", width="3", command=onTimes, 
        font=.operatorFont)
    colonButton <- tkbutton(x3Frame, text=":", width="3", command=onColon, 
        font=.operatorFont)
    slashButton <- tkbutton(x3Frame, text="/", width="3", command=onSlash, 
        font=.operatorFont)
    inButton <- tkbutton(xFrame, text="%in%", width="3", command=onIn,
        font=.operatorFont)
    minusButton <- tkbutton(x3Frame, text="I(", width="3", command=onMinus, 
        font=.operatorFont)
    minus2Button <- tkbutton(x3Frame, text="s(", width="3", command=onMinus2, 
        font=.operatorFont)
    powerButton <- tkbutton(x3Frame, text="^", width="3", command=onPower, 
        font=.operatorFont)
    leftParenButton <- tkbutton(x3Frame, text="(", width="3", command=onLeftParen, 
        font=.operatorFont)
    rightParenButton <- tkbutton(x3Frame, text=")", width="3", command=onRightParen, 
        font=.operatorFont)
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    plotButton <- tkbutton(buttonsFrame, text="Plot", width="12", command=onPlot)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        help(generalizedLinearModel)
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    tkgrid(tklabel(modelFrame, text="Save model as:               ", width=20), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(tklabel(option1Frame, text="Model options"), sticky="w")
    tkgrid(optionBox, optionScroll,sticky="w")
    tkgrid(standardCheckBox, tklabel(option2Frame, text="standardise"), sticky="w")
    tkgrid(summaryCheckBox, tklabel(option2Frame, text="print summary"), sticky="w")
    tkgrid(anovaCheckBox, tklabel(option2Frame, text="print anova"), sticky="w")
    tkgrid(dataCheckBox, tklabel(option2Frame, text="add predictions to dataframe"), sticky="w")
    tkgrid(option1Frame, tklabel(optionFrame, text="", width=1), option2Frame, sticky="w")
    tkgrid(optionFrame, sticky="w")
    tkgrid(tklabel(lhsFrame, text="Response"), sticky="w")
    tkgrid(lhsEntry, sticky="nw")
    tkgrid(yBox, yScroll, sticky="nw")
    tkgrid(lhsFrame,tklabel(resFrame, text="", width=1), yFrame)
    tkgrid(resFrame, sticky="w")
    tkgrid(rhsEntry, sticky="w")
    tkgrid(xBox, xScroll,sticky="w") 
    tkgrid(plusButton, timesButton, colonButton, slashButton, inButton, sticky="w")
    tkgrid(minusButton,powerButton, leftParenButton, rightParenButton, minus2Button, sticky="w")
    tkgrid(tklabel(xFrame, text="Explanatory"), sticky="w")
    tkgrid(x1Frame, sticky="w")
    tkgrid(x2Frame, tklabel(xFrame, text="", width=1), x3Frame, sticky="w")
    tkgrid(x4Frame, sticky="w")
    tkgrid(xFrame, sticky="w")
    tkgrid(tklabel(subsetFrame, text="Remove sites with name:     ", width=20), subsetEntry, sticky="w")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(tklabel(plot1Frame, text="Plot options"), sticky="w")
    tkgrid(typeBox, typeScroll, sticky="nw")
    tkgrid(tklabel(plot2Frame, text="Plot variable"), sticky="w")
    tkgrid(axisBox, axisScroll, sticky="nw")
    tkgrid(plot1Frame, tklabel(plotFrame, text="", width=1), plot2Frame, sticky="w")
    tkgrid(plotFrame, sticky="w")
    tkgrid(OKbutton, plotButton, cancelButton, tklabel(buttonsFrame, text="            "), 
        helpButton, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(xScroll, sticky="ns")
    tkgrid.configure(yScroll, sticky="ns")
    tkgrid.configure(typeScroll, sticky="ns")
    tkgrid.configure(axisScroll, sticky="ns")
    tkgrid.configure(optionScroll, sticky="ns")
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkselection.set(typeBox, 0)
    tkselection.set(optionBox, 0)
    tkselection.set(axisBox, 0)
    tkbind(top, "<Return>", onOK)
    tkbind(xBox, "<Double-ButtonPress-1>", onDoubleClick)
    tkbind(yBox, "<Double-ButtonPress-1>", onDoubleClick2)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(lhsEntry)
    tkwait.window(top)
} 

diversityresponse <- function(){
    .communityDataSet <- CommunityDataSet()
    .activeDataSet <- ActiveDataSet()
    justDoIt(paste(.activeDataSet, " <- diversityvariables(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
    logger(paste(.activeDataSet, " <- diversityvariables(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
    activeDataSet(.activeDataSet)
    communityDataSet(.activeDataSet)
    logger(paste("environmental data set (with added diversity variables) is now also the community data set", sep=""))
    logger(paste("use menu option of 'Species abundance as response...' now to analyze diversity as response", sep=""))
}

distmatrixGUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Distance matrix calculation")
    .activeDataSet <- ActiveDataSet()
    .communityDataSet <- CommunityDataSet()
    modelName <- tclVar("Distmatrix.1")
    modelFrame <- tkframe(top, relief="groove", borderwidth=2)
    model <- tkentry(modelFrame, width=50, textvariable=modelName)
    method2Frame <- tkframe(top, relief="groove", borderwidth=2)
    distBox <- tklistbox(method2Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    distScroll <- tkscrollbar(method2Frame, repeatinterval=5, command=function(...) tkyview(distBox, ...))
    treatasdistVariable <- tclVar("0")
    treatasdistCheckBox <- tkcheckbutton(method2Frame, variable=treatasdistVariable)
    tkconfigure(distBox, yscrollcommand=function(...) tkset(distScroll, ...))
    distances <- c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", 
        "chao", "cao", "mahalanobis",
        "w", "-1", "c", "wb", "r", "I", "e", "t", "me", "j", "sor", "m", "-2", "co", "cc", "g", "-3", "l", "19", "hk", "rlb", "sim", "gl", "z")
    for (x in distances) tkinsert(distBox, "end", x)
    onOK <- function(){
        dist <- distances[as.numeric(tkcurselection(distBox))+1]
        modelValue <- tclvalue(modelName)
        if (dist %in% c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", 
            "chao", "cao", "mahalanobis")) {
            logger(paste(modelValue, " <- vegdist(", .communityDataSet, ", method='", dist, "', na.rm=T)", sep=""))
            assign(modelValue, justDoIt(paste("vegdist(",.communityDataSet, ", method='",dist, "', na.rm=T)", sep="")), envir=.GlobalEnv)
        }else{
            logger(paste(modelValue, " <- betadiver(", .communityDataSet, ", method='", dist, "')", sep=""))
            assign(modelValue, justDoIt(paste("betadiver(",.communityDataSet, ", method='",dist, "')", sep="")), envir=.GlobalEnv)
        }
        doItAndPrint(paste(modelValue))
        doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist, "')", sep=""))
        treatasdist <- tclvalue(treatasdistVariable)==1
        if (treatasdist==T) {
            logger(paste(modelValue, " <- data.frame(as.matrix(", modelValue, "))", sep=""))
            assign(modelValue, justDoIt(paste("data.frame(as.matrix(", modelValue, "))", sep="")), envir=.GlobalEnv)
            communityDataSet(modelValue) 
        }
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
    }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(modelFrame, text="Save data as:", width=10), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(tklabel(method2Frame, text="Distance"), sticky="w")
    tkgrid(distBox, distScroll,sticky="w")
    tkgrid(tklabel(method2Frame, text="Make community dataset", width=25),treatasdistCheckBox, sticky="w")
    tkgrid(method2Frame, sticky="w")
    tkgrid(OKbutton, cancelButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(distScroll, sticky="ns")
    tkselection.set(distBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(distBox)
    tkwait.window(top)
}

unconordiGUI <- function(){
    contrasts <- c("contr.treatment", "contr.poly")
    .activeDataSet <- ActiveDataSet()
    .communityDataSet <- CommunityDataSet()
    .variables <- Variables()
    variables <- paste(.variables, ifelse(is.element(.variables, Factors()), "[factor]", ""))
    .cvariables <- CVariables()
    cvariables <- paste(.cvariables)
    top <- tktoplevel()
    tkwm.title(top, "Unconstrained ordination")
    modelName <- tclVar("Ordination.model1")
    modelFrame <- tkframe(top, relief="groove", borderwidth=2)
    model <- tkentry(modelFrame, width=40, textvariable=modelName)
    methodFrame <- tkframe(top, relief="groove", borderwidth=2)
    method1Frame <- tkframe(methodFrame)
    method2Frame <- tkframe(methodFrame)
    method3Frame <- tkframe(methodFrame)
    method4Frame <- tkframe(methodFrame)
    methodBox <- tklistbox(method1Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    methodScroll <- tkscrollbar(method1Frame, repeatinterval=5, command=function(...) tkyview(methodBox, ...))
    tkconfigure(methodBox, yscrollcommand=function(...) tkset(methodScroll, ...))
    methods <- c("PCA","PCA (prcomp)", "PCoA","PCoA (Caillez)","CA","DCA","metaMDS","wcmdscale","pcnm","NMS (standard)","isomap")
    for (x in methods) tkinsert(methodBox, "end", x)
    distBox <- tklistbox(method2Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    distScroll <- tkscrollbar(method2Frame, repeatinterval=5, command=function(...) tkyview(distBox, ...))
    tkconfigure(distBox, yscrollcommand=function(...) tkset(distScroll, ...))
    distances <- c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", 
            "chao", "cao", "mahalanobis")
    for (x in distances) tkinsert(distBox, "end", x)
    summaryVariable <- tclVar("1")
    summaryCheckBox <- tkcheckbutton(method4Frame, variable=summaryVariable)
    scalingVariable <- tclVar("1")
    scale <- tkentry(method4Frame, width=10, textvariable=scalingVariable)
    NMSVariable <- tclVar("2")
    NMSa <- tkentry(method3Frame, width=10, textvariable=NMSVariable)
    NMSpermVariable <- tclVar("1")
    NMSperm <- tkentry(method3Frame, width=10, textvariable=NMSpermVariable)
    addspecVariable <- tclVar("0")
    addspecCheckBox <- tkcheckbutton(method3Frame, variable=addspecVariable)
    treatasdistVariable <- tclVar("0")
    treatasdistCheckBox <- tkcheckbutton(method4Frame, variable=treatasdistVariable)
    plotFrame <- tkframe(top, relief="groove", borderwidth=2)
    plot1Frame <- tkframe(plotFrame)
    plot2Frame <- tkframe(plotFrame)
    plot3Frame <- tkframe(plotFrame)
    plot4Frame <- tkframe(plotFrame)
    typeBox <- tklistbox(plot1Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    typeScroll <- tkscrollbar(plot1Frame, repeatinterval=5, command=function(...) tkyview(typeBox, ...))
    tkconfigure(typeBox, yscrollcommand=function(...) tkset(typeScroll, ...))
    types <- c("plot","ordiplot","ordiplot empty","origin axes","identify sites","identify species","text sites","text species","points sites","points species",
        "label sites","label species","orditorp sites","orditorp species",
        "envfit","ordihull (factor)", "ordiarrows (factor)","ordisegments (factor)","ordispider (factor)","ordiellipse (factor)","ordisurf (continuous)",
        "ordibubble (continuous)","ordisymbol (factor)","ordisymbol (factor, legend)","ordivector (species)","ordivector interpretation",
        "ordicluster","ordicluster2","ordispantree","ordinearest","ordiequilibriumcircle","screeplot",
        "distance displayed","coenocline","stressplot",
        "orditkplot sites","orditkplot species","orditkplot pointlabel")
    for (x in types) tkinsert(typeBox, "end", x)
    choicesVariable <- tclVar("1,2")
    choice <- tkentry(plot3Frame, width=10, textvariable=choicesVariable)
    dataVariable <- tclVar("0")
    dataCheckBox <- tkcheckbutton(plot3Frame, variable=dataVariable)
    axisBox <- tklistbox(plot2Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    axisScroll <- tkscrollbar(plot2Frame, repeatinterval=5, command=function(...) tkyview(axisBox, ...))
    tkconfigure(axisBox, yscrollcommand=function(...) tkset(axisScroll, ...))
    for (x in variables) tkinsert(axisBox, "end", x)
    cexVariable <- tclVar("1")
    cexa <- tkentry(plot4Frame, width=10, textvariable=cexVariable)
    colVariable <- tclVar("blue")
    cola <- tkentry(plot4Frame, width=10, textvariable=colVariable)
    onOK <- function(){
        doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
        method <- methods[as.numeric(tkcurselection(methodBox))+1]
        dist <- distances[as.numeric(tkcurselection(distBox))+1]
        k <- tclvalue(NMSVariable)
        perm <- tclvalue(NMSpermVariable)
        treatasdist <- tclvalue(treatasdistVariable)==1
        addspec <- tclvalue(addspecVariable) == "1"
        if (method=="PCA") {
            command <- paste("rda(", .communityDataSet, ")", sep="")
            doItAndPrint(paste("dist.eval(", .communityDataSet, ",'euclidean')", sep=""))
        }
        if (method=="PCA (prcomp)") {
            command <- paste("prcomp(", .communityDataSet, ")", sep="")
            doItAndPrint(paste("dist.eval(", .communityDataSet, ",'euclidean')", sep=""))
        }
        if (method=="CA") {
            command <- paste("cca(", .communityDataSet, ")", sep="")
        }
        if (method=="DCA") {
            command <- paste("decorana(", .communityDataSet, ")", sep="")
        }
        if (method=="PCoA") {
            if(treatasdist==F){
                logger(paste("distmatrix <- vegdist(", .communityDataSet, ",method='", dist, "', na.rm=T)", sep=""))
                assign("distmatrix", justDoIt(paste("vegdist(",.communityDataSet, ",method='",dist, "', na.rm=T)", sep="")), envir=.GlobalEnv)
                doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist, "')", sep=""))
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            command <- paste("cmdscale(distmatrix, k=", k, ", eig=T, add=F)", sep="")
        }
        if (method=="PCoA (Caillez)") {
            if(treatasdist==F){
                logger(paste("distmatrix <- vegdist(", .communityDataSet, ",method='", dist, "', na.rm=T)", sep=""))
                assign("distmatrix", justDoIt(paste("vegdist(",.communityDataSet, ",method='",dist, "', na.rm=T)", sep="")), envir=.GlobalEnv)
                doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist, "')", sep=""))
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            command <- paste("cmdscale(distmatrix, k=", k, ", eig=T, add=T)", sep="")
        }
        if (method=="metaMDS") { 
            addspec <- F          
            command <- paste("metaMDS(", .communityDataSet, ",distance='", dist, "', k=", k, ", trymax=", perm, ", autotransform=T, noshare=0.1, expand=T, trace=1, plot=F)", sep="")
            doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist, "')", sep=""))
        }
        if (method=="wcmdscale") {
            if(treatasdist==F){
                logger(paste("distmatrix <- vegdist(", .communityDataSet, ", method='", dist, "', na.rm=T)", sep=""))
                assign("distmatrix", justDoIt(paste("vegdist(",.communityDataSet, ", method='",dist, "', na.rm=T)", sep="")), envir=.GlobalEnv)
                doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist, "')", sep=""))
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            command <- paste("wcmdscale(distmatrix, k=", k, ", eig=T)", sep="")          
        }
        if (method=="pcnm") {
            if(treatasdist==F){
                logger(paste("distmatrix <- vegdist(", .communityDataSet, ", method='", dist, "', na.rm=T)", sep=""))
                assign("distmatrix", justDoIt(paste("vegdist(",.communityDataSet, ", method='",dist, "', na.rm=T)", sep="")), envir=.GlobalEnv)
                doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist, "')", sep=""))
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            command <- paste("pcnm(distmatrix)", sep="")          
        }
        if (method=="NMS (standard)") {
            if(treatasdist==F){
                logger(paste("distmatrix <- vegdist(", .communityDataSet, ", method='", dist, "', na.rm=T)", sep=""))
                assign("distmatrix", justDoIt(paste("vegdist(",.communityDataSet, ", method='",dist, "', na.rm=T)", sep="")), envir=.GlobalEnv)
                doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist, "')", sep=""))
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            command <- paste("NMSrandom(distmatrix, perm=", perm,", k=", k, ")", sep="")          
        }
        if (method=="isomap") {
            if(treatasdist==F){
                logger(paste("distmatrix <- vegdist(", .communityDataSet, ", method='", dist, "', na.rm=T)", sep=""))
                assign("distmatrix", justDoIt(paste("vegdist(",.communityDataSet, ", method='",dist, "', na.rm=T)", sep="")), envir=.GlobalEnv)
                doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist, "')", sep=""))
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            command <- paste("isomap(distmatrix, k=4)", sep="")          
        }
        modelValue <- tclvalue(modelName)
#        if (!is.valid.name(modelValue)){
#            tkmessageBox(message=paste('"', modelValue, '" is not a valid name.', 
#                sep=""), icon="error", type="ok")
#        }        
        logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        sum <- tclvalue(summaryVariable) == "1"
        scaling <- tclvalue(scalingVariable)
        if (method == "PCoA"  || method == "PCoA (Caillez)") {
            doItAndPrint(paste("rownames(", modelValue, "$points) <- rownames(", .communityDataSet, ")", sep=""))
        }
        if (addspec==T) {
            if (method=="PCoA") {
                doItAndPrint(paste(modelValue, "<- add.spec.scores(", modelValue, ",", .communityDataSet, ", method='pcoa.scores', Rscale=T, scaling=", scaling, ", multi=1)", sep=""))
            }
            if (method=="PCoA (Caillez)") {
                doItAndPrint(paste(modelValue, "<- add.spec.scores(", modelValue, ",", .communityDataSet, ", method='pcoa.scores', Rscale=T, scaling=", scaling, ", multi=1)", sep=""))
            }
            if (method=="NMS (standard)") {
                doItAndPrint(paste(modelValue, "<- add.spec.scores(", modelValue, ",", .communityDataSet, ", method='wa.scores')", sep=""))
            }
        }
        if (method == "PCA" || method == "PCA (prcomp)"  || method == "CA" || method == "DCA" || addspec == T) {
            doItAndPrint(paste("check.ordiscores(", .communityDataSet, ",", modelValue, ", check.species=T)", sep=""))
        }else{
            doItAndPrint(paste("check.ordiscores(", .communityDataSet, ",", modelValue, ", check.species=F)", sep=""))
        }
        if (sum==T) {
            if (method %in% c("PCA", "CA", "DCA")) {
                doItAndPrint(paste("summary(", modelValue, ", scaling=", scaling, ")", sep=""))
                if (method=="PCA") {doItAndPrint(paste("PCAsignificance(", modelValue, ")", sep=""))}
                doItAndPrint(paste("goodness(", modelValue, ", display='sites', choices=c(1:4), statistic='explained')", sep=""))
                doItAndPrint(paste("inertcomp(", modelValue, ", display='sites', statistic='explained', proportional=F)", sep=""))
            }         
            if (method %in% c("PCA (prcomp)", "PCoA","PCoA (Caillez)","metaMDS","wcmdscale","pcnm","NMS (standard)")) {
                doItAndPrint(paste(modelValue, sep=""))
                if (method=="metaMDS") {doItAndPrint(paste("goodness(", modelValue, ")", sep=""))}
            }
            if (method=="isomap") {doItAndPrint(paste("summary(", modelValue, ")", sep=""))}
        }
    }
    onPlot <- function(){
        method <- methods[as.numeric(tkcurselection(methodBox))+1]
        modelValue <- tclvalue(modelName)
        plottype <- types[as.numeric(tkcurselection(typeBox))+1]
        scaling <- tclvalue(scalingVariable)
        perm <- tclvalue(NMSpermVariable)
        choices <- tclvalue(choicesVariable)
        dist <- distances[as.numeric(tkcurselection(distBox))+1]
        col <- tclvalue(colVariable)
        cex <- tclvalue(cexVariable)
        treatasdist <- tclvalue(treatasdistVariable)==1
        addspec <- tclvalue(addspecVariable) == "1"
        justDoIt(paste("attach(", .activeDataSet, ")", sep=""))
        if (plottype == "plot"){
            justDoIt(paste("par(cex=",cex,")", sep=""))
            logger(paste("par(cex=",cex,")", sep=""))
            if (method=="PCA" || method=="CA"  || method=="DCA") {
                justDoIt(paste("plot1 <- plot(", modelValue, ", choices=c(", choices, "), scaling=", scaling, ")", sep=""))
                logger(paste("plot1 <- plot(", modelValue, ", choices=c(", choices, "), scaling=", scaling, ")", sep=""))
            }
            if (method=="metaMDS") {
                justDoIt(paste("plot1 <- plot(", modelValue, ", choices=c(", choices, "))", sep=""))
                logger(paste("plot1 <- plot(", modelValue, ", choices=c(", choices, "))", sep=""))
            }
            if (method=="PCoA" || method=="PCoA (Caillez)"  || method=="NMS (standard)" || method=="wcmdscale" || method=="pcnm") {
                justDoIt(paste("plot1 <- plot(scores(", modelValue, ", display='sites', choices=c(", choices, ")))", sep=""))
                logger(paste("plot1 <- plot(scores(", modelValue, ", display='sites', choices=c(", choices, ")))", sep=""))
                if (addspec==T) {
                    justDoIt(paste("plot1 <- plot(scores(", modelValue, ", display='species', choices=c(", choices, ")), pch='+', col='red')", sep=""))
                    logger(paste("plot1 <- plot(scores(", modelValue, ", display='species', choices=c(", choices, ")), pch='+', col='red')", sep=""))
                }
                justDoIt(paste("text(scores(", modelValue, ", display='sites', choices=c(", choices, ")), rownames(", .communityDataSet, "), pos=3)", sep=""))
                logger(paste("text(scores(", modelValue, ", display='sites', choices=c(", choices, ")), rownames(", .communityDataSet, "), pos=3)", sep=""))
            }
        }
        if (plottype == "ordiplot"){
            justDoIt(paste("par(cex=",cex,")", sep=""))
            logger(paste("par(cex=",cex,")", sep=""))
            if (method=="PCA" || method=="CA"  || method=="DCA") {
                justDoIt(paste("plot1 <- ordiplot(", modelValue, ", choices=c(", choices, "), scaling=", scaling, ")", sep=""))
                logger(paste("plot1 <- ordiplot(", modelValue, ", choices=c(", choices, "), scaling=", scaling, ")", sep=""))
            }else{
                justDoIt(paste("plot1 <- ordiplot(", modelValue, ", choices=c(", choices, "))", sep="")) 
                logger(paste("plot1 <- ordiplot(", modelValue, ", choices=c(", choices, "))", sep="")) 
            }
        }
        if (plottype == "ordiplot empty"){
            justDoIt(paste("par(cex=",cex,")", sep=""))
            logger(paste("par(cex=",cex,")", sep=""))
            if (method=="PCA" || method=="CA"  || method=="DCA") {
                justDoIt(paste("plot1 <- ordiplot(", modelValue, ", type='none', choices=c(", choices, "), scaling=", scaling, ")", sep=""))
                logger(paste("plot1 <- ordiplot(", modelValue, ", type='none', choices=c(", choices, "), scaling=", scaling, ")", sep=""))
            }else{
                justDoIt(paste("plot1 <- ordiplot(", modelValue, ", type='none', choices=c(", choices, "))", sep="")) 
                logger(paste("plot1 <- ordiplot(", modelValue, ", type='none', choices=c(", choices, "))", sep="")) 
            }
        }
        if (plottype == "identify sites"){
            doItAndPrint(paste("identify(plot1, 'sites', col='", col,"', cex=", cex, ")", sep=""))             
        }
        if (plottype == "identify species"){
            doItAndPrint(paste("identify(plot1, 'species', col='", col,"', cex=", cex, ")", sep=""))             
        }
        if (plottype == "text sites"){
            doItAndPrint(paste("text(plot1, 'sites', col='", col,"', cex=", cex, ")",sep=""))             
        }
        if (plottype == "text species"){
            doItAndPrint(paste("text(plot1, 'species', col='", col,"', cex=", cex, ")",sep=""))             
        }
        if (plottype == "points sites"){
            doItAndPrint(paste("points(plot1, 'sites', col='", col,"', cex=", cex, ")",sep=""))             
        }
        if (plottype == "points species"){
            doItAndPrint(paste("points(plot1, 'species', col='", col,"', cex=", cex, ")",sep=""))             
        }
        if (plottype == "label sites"){
            doItAndPrint(paste("ordilabel(plot1, 'sites', col='", col,"', cex=", cex, ")",sep=""))             
        }
        if (plottype == "label species"){
            doItAndPrint(paste("ordilabel(plot1, 'species', col='", col,"', cex=", cex, ")",sep=""))             
        }
        if (plottype == "orditorp sites"){
            doItAndPrint(paste("orditorp(plot1, 'sites', col='", col,"', cex=", cex, ")",sep=""))             
        }
        if (plottype == "orditorp species"){
            doItAndPrint(paste("orditorp(plot1, 'species', col='", col,"', cex=", cex, ")",sep=""))             
        }
        if (plottype == "origin axes"){
            doItAndPrint(paste("abline(h = 0, lty = 3)", sep=""))
            doItAndPrint(paste("abline(v = 0, lty = 3)", sep=""))
        }
        if (plottype == "screeplot"){
            justDoIt(paste("par(cex=",cex,")", sep=""))
            logger(paste("par(cex=",cex,")", sep=""))
#           if (method=="PCA" || method=="PCA (prcomp)") {
                doItAndPrint(paste("plot1 <- screeplot(", modelValue, ", bstick=T)", sep=""))
#            }
        }
        axisvar <- .variables[as.numeric(tkcurselection(axisBox))+1]
        varfactor <- eval(parse(text=paste("is.factor(",.activeDataSet, "$", axisvar, ")", sep="")), envir=.GlobalEnv)
        if (plottype == "envfit"){
            doItAndPrint(paste("fitted <- envfit(plot1, data.frame(", axisvar, "), permutations=", perm, ")", sep=""))
            doItAndPrint(paste("plot(fitted, col='", col,"', cex=", cex, ")", sep=""))
        }
        if (plottype == "ordihull (factor)" && varfactor==T){
            doItAndPrint(paste("ordihull(plot1,", axisvar, ", col='", col, "')", sep=""))
        }
        if (plottype == "ordiarrows (factor)" && varfactor==T){
            doItAndPrint(paste("ordiarrows(plot1,", axisvar, ", col='", col, "')", sep=""))
        }
        if (plottype == "ordisegments (factor)" && varfactor==T){
            doItAndPrint(paste("ordisegments(plot1,", axisvar, ", col='", col, "')", sep=""))
        }
        if (plottype == "ordispider (factor)" && varfactor==T){
            doItAndPrint(paste("ordispider(plot1,", axisvar, ", col='", col, "')", sep=""))
            }
        if (plottype == "ordiellipse (factor)" && varfactor==T){
            doItAndPrint(paste("ordiellipse(plot1,", axisvar, ", col='", col, "', conf=0.9)", sep=""))
            }
        if (plottype == "ordisurf (continuous)" && varfactor==F){
            doItAndPrint(paste("ordisurf(plot1,", axisvar, ", add=T)", sep=""))
            }
        if (plottype == "ordibubble (continuous)" && varfactor==F){
            doItAndPrint(paste("ordibubble(plot1,", axisvar, ", fg='", col, "')", sep=""))
            }
        if (plottype == "ordisymbol (factor)" && varfactor==T){
            justDoIt(paste("ordisymbol(plot1,", .activeDataSet, ",'", axisvar, "', legend=F, rainbow=T, cex=", cex, ")", sep=""))
            logger(paste("ordisymbol(plot1,", .activeDataSet, ",'", axisvar, "', legend=F, rainbow=T, cex=", cex, ")", sep=""))
        }
        if (plottype == "ordisymbol (factor, legend)" && varfactor==T){
            justDoIt(paste("ordisymbol(plot1,", .activeDataSet, ",'", axisvar, "', legend=T, legend.x='topleft', legend.ncol=1, rainbow=T, cex=", cex, ")", sep=""))
            logger(paste("ordisymbol(plot1,", .activeDataSet, ",'", axisvar, "', legend=T, legend.x='topleft', legend.ncol=1, rainbow=T, cex=", cex, ")", sep=""))
        }
        if (plottype == "ordivector (species)"){
            realspecies <- eval(parse(text=paste("any(colnames(", .communityDataSet, ")=='", axisvar, "')", sep="")), envir=.GlobalEnv)
            if (realspecies == T) {
                doItAndPrint(paste("ordivector(plot1,'", axisvar, "',lty=0, angle=5, length=0.5)", sep=""))
            }
        }
        if (plottype == "ordivector interpretation"){
            realspecies <- eval(parse(text=paste("any(colnames(", .communityDataSet, ")=='", axisvar, "')", sep="")), envir=.GlobalEnv)
            if (realspecies == T) {
                doItAndPrint(paste("ordivector(plot1,'", axisvar, "',lty=2)", sep=""))
            }
        }
        if (plottype == "ordiequilibriumcircle" && method == "PCA"){
            doItAndPrint(paste("ordiequilibriumcircle(", modelValue, ",plot1, col='", col, "')", sep=""))
            }            
        if (plottype == "ordicluster"){
            if(treatasdist==F){
                logger(paste("distmatrix <- vegdist(", .communityDataSet, ",method='", dist, "', na.rm=T)", sep=""))
                assign("distmatrix", justDoIt(paste("vegdist(",.communityDataSet, ",method='",dist, "', na.rm=T)", sep="")), envir=.GlobalEnv)
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            logger(paste("cluster <- hclust(distmatrix, method='single')", sep=""))
            assign("cluster", justDoIt(paste("hclust(distmatrix, method='single')", sep="")), envir=.GlobalEnv)
            doItAndPrint(paste("ordicluster(plot1, cluster, prune=0, col='", col, "')", sep=""))
            }
        if (plottype == "ordicluster2"){
            if(treatasdist==F){
                logger(paste("distmatrix <- vegdist(", .communityDataSet, ",method='", dist, "', na.rm=T)", sep=""))
                assign("distmatrix", justDoIt(paste("vegdist(",.communityDataSet, ",method='",dist, "', na.rm=T)", sep="")), envir=.GlobalEnv)
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            logger(paste("cluster <- hclust(distmatrix, method='single')", sep=""))
            assign("cluster", justDoIt(paste("hclust(distmatrix, method='single')", sep="")), envir=.GlobalEnv)
            doItAndPrint(paste("ordicluster2(plot1, cluster, mingroups=1, col='", col, "')", sep=""))
            }
        if (plottype == "ordinearest"){
            if(treatasdist==F){
                logger(paste("distmatrix <- vegdist(", .communityDataSet, ",method='", dist, "', na.rm=T)", sep=""))
                assign("distmatrix", justDoIt(paste("vegdist(",.communityDataSet, ",method='",dist, "', na.rm=T)", sep="")), envir=.GlobalEnv)
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            doItAndPrint(paste("ordinearest(plot1, distmatrix,col='", col, "')", sep=""))
            }
        if (plottype == "ordispantree"){
            if(treatasdist==F){
                logger(paste("distmatrix <- vegdist(", .communityDataSet, ",method='", dist, "', na.rm=T)", sep=""))
                assign("distmatrix", justDoIt(paste("vegdist(",.communityDataSet, ",method='",dist, "', na.rm=T)", sep="")), envir=.GlobalEnv)
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            doItAndPrint(paste("lines(spantree(distmatrix,toolong=0),plot1,col='", col, "')", sep=""))
            }
        if (plottype == "distance displayed"){
            if(treatasdist==F){
                doItAndPrint(paste("distdisplayed(", .communityDataSet ,",plot1, distx='", dist, "',plotit=T)", sep=""))
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
                doItAndPrint(paste("distdisplayed(distmatrix, plot1, plotit=T)", sep=""))
            }
        }
        if (plottype == "coenocline"){
            doItAndPrint(paste("ordicoeno(", .communityDataSet ,", ordiplot=plot1, axis=1, legend=T, cex=0.8, ncol=4)", sep=""))
            }
        if (plottype == "stressplot"  && method == "metaMDS"){
            doItAndPrint(paste("stressplot(", modelValue ,")", sep=""))
            }
        if (plottype == "orditkplot sites"){
            justDoIt(paste("par(cex=",cex,")", sep=""))
            logger(paste("par(cex=",cex,")", sep=""))
            if (method=="PCA" || method=="CA"  || method=="DCA") {
                justDoIt(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", display='sites', choices=c(", choices, "), scaling=", scaling, "))", sep=""))
                logger(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", display='sites', choices=c(", choices, "), scaling=", scaling, "))", sep=""))
            }else{
                justDoIt(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", display='sites', choices=c(", choices, ")))", sep="")) 
                logger(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", display='sites', choices=c(", choices, ")))", sep="")) 
            }
        }
        if (plottype == "orditkplot species"){
            justDoIt(paste("par(cex=",cex,")", sep=""))
            logger(paste("par(cex=",cex,")", sep=""))
            if (method=="PCA" || method=="CA"  || method=="DCA") {
                justDoIt(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", display='species', choices=c(", choices, "), scaling=", scaling, "))", sep=""))
                logger(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", display='species', choices=c(", choices, "), scaling=", scaling, "))", sep=""))
            }else{
                justDoIt(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", display='species', choices=c(", choices, ")))", sep="")) 
                logger(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", display='species', choices=c(", choices, ")))", sep="")) 
            }
        }
        if (plottype == "orditkplot pointlabel"){
            justDoIt(paste("par(cex=",cex,")", sep=""))
            logger(paste("par(cex=",cex,")", sep=""))
            if (method=="PCA" || method=="CA"  || method=="DCA") {
                justDoIt(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", choices=c(", choices, "), scaling=", scaling, "))", sep=""))
                logger(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", choices=c(", choices, "), scaling=", scaling, "))", sep=""))
            }else{
                justDoIt(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", choices=c(", choices, ")))", sep="")) 
                logger(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", choices=c(", choices, ")))", sep="")) 
            }
        }

        data <- tclvalue(dataVariable) =="1"
        if (data==T) {
            justDoIt(paste(.activeDataSet, "$", modelValue, ".ax1 <- scores(plot1,display='sites')[,1]", sep=""))
            logger(paste(.activeDataSet, "$", modelValue, ".ax1 <- scores(plot1,display='sites')[,1]", sep=""))
            justDoIt(paste(.activeDataSet, "$", modelValue, ".ax2 <- scores(plot1,display='sites')[,2]", sep=""))
            logger(paste(.activeDataSet, "$", modelValue, ".ax2 <- scores(plot1,display='sites')[,2]", sep=""))
            activeDataSet(.activeDataSet)
        }
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    plotButton <- tkbutton(buttonsFrame, text="Plot", width="12", command=onPlot)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(modelFrame, text="Save model as:               ", width=20), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(tklabel(method1Frame, text="Ordination method"), sticky="w")
    tkgrid(methodBox, methodScroll,sticky="w")
    tkgrid(tklabel(method2Frame, text="Distance"), sticky="w")
    tkgrid(distBox, distScroll,sticky="w")
    tkgrid(tklabel(method3Frame, text="PCoA/NMS axes", width=15), NMSa, sticky="w")
    tkgrid(tklabel(method3Frame, text="NMS permutations", width=15), NMSperm, sticky="w")
    tkgrid(addspecCheckBox, tklabel(method3Frame, text="PCoa/NMS species", width=15), sticky="w")
    tkgrid(summaryCheckBox, tklabel(method4Frame, text="model summary"), sticky="w")
    tkgrid(tklabel(method4Frame, text="scaling", width=10), scale, sticky="w")
    tkgrid(treatasdistCheckBox, tklabel(method4Frame, text="as.dist(Community)", width=15), sticky="w")
    tkgrid(method1Frame, tklabel(methodFrame, text="", width=1), method2Frame, sticky="w")
    tkgrid(method3Frame, tklabel(methodFrame, text="", width=1), method4Frame, sticky="w")
    tkgrid(methodFrame, sticky="w")
    tkgrid(tklabel(plot1Frame, text="Plot method"), sticky="w")
    tkgrid(typeBox, typeScroll, sticky="nw")
    tkgrid(tklabel(plot2Frame, text="Plot variable"), sticky="w")
    tkgrid(axisBox, axisScroll, sticky="nw")
    tkgrid(tklabel(plot3Frame, text="axes", width=10), choice, sticky="w")
    tkgrid(dataCheckBox, tklabel(plot3Frame, text="add scores to dataframe"), sticky="w")
    tkgrid(tklabel(plot4Frame, text="cex", width=10), cexa, sticky="w")
    tkgrid(tklabel(plot4Frame, text="colour", width=10), cola, sticky="w")
    tkgrid(plot1Frame, tklabel(plotFrame, text="", width=1), plot2Frame, sticky="w")
    tkgrid(plot3Frame, tklabel(plotFrame, text="", width=1), plot4Frame, sticky="w")
    tkgrid(plotFrame, sticky="w")
    tkgrid(OKbutton, plotButton, cancelButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(typeScroll, sticky="ns")
    tkgrid.configure(axisScroll, sticky="ns")
    tkgrid.configure(methodScroll, sticky="ns")
    tkgrid.configure(distScroll, sticky="ns")
    tkselection.set(typeBox, 0)
    tkselection.set(methodBox, 0)
    tkselection.set(axisBox, 0)
    tkselection.set(distBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(methodBox)
    tkwait.window(top)
    }


conordiGUI <- function(){
    contrasts <- c("contr.treatment", "contr.poly")
    checkAddOperator <- function(rhs){
        rhs.chars <- rev(strsplit(rhs, "")[[1]])
        if (length(rhs.chars) < 1) return(FALSE)
        check.char <- if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1)) 
                rhs.chars[1] else rhs.chars[2]
        !is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%"))
        }
    top <- tktoplevel()
    tkwm.title(top, "Constrained ordination")
    .activeDataSet <- ActiveDataSet()
    .communityDataSet <- CommunityDataSet()
    .variables <- Variables()
    variables <- paste(.variables, ifelse(is.element(.variables, Factors()), "[factor]", ""))
    .cvariables <- CVariables()
    cvariables <- paste(.cvariables)
    modelName <- tclVar("Ordination.model1")
    modelFrame <- tkframe(top, relief="groove", borderwidth=2)
    model <- tkentry(modelFrame, width=40, textvariable=modelName)
    methodFrame <- tkframe(top, relief="groove", borderwidth=2)
    method1Frame <- tkframe(methodFrame)
    method2Frame <- tkframe(methodFrame)
    method3Frame <- tkframe(methodFrame)
    method4Frame <- tkframe(methodFrame)
    methodBox <- tklistbox(method1Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    methodScroll <- tkscrollbar(method1Frame, repeatinterval=5, command=function(...) tkyview(methodBox, ...))
    tkconfigure(methodBox, yscrollcommand=function(...) tkset(methodScroll, ...))
    methods <- c("RDA","CCA","capscale","capscale(add)","CAPdiscrim","prc","multiconstrained (RDA)", "multiconstrained (CCA)", "multiconstrained (capscale)", "multiconstrained (capscale add)")
    for (x in methods) tkinsert(methodBox, "end", x)
    distBox <- tklistbox(method2Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    distScroll <- tkscrollbar(method2Frame, repeatinterval=5, command=function(...) tkyview(distBox, ...))
    tkconfigure(distBox, yscrollcommand=function(...) tkset(distScroll, ...))
    distances <- c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", 
            "chao", "cao", "mahalanobis")
    for (x in distances) tkinsert(distBox, "end", x)
    summaryVariable <- tclVar("1")
    summaryCheckBox <- tkcheckbutton(method3Frame, variable=summaryVariable)
    treatasdistVariable <- tclVar("0")
    treatasdistCheckBox <- tkcheckbutton(method3Frame, variable=treatasdistVariable)
    scalingVariable <- tclVar("1")
    scale <- tkentry(method4Frame, width=10, textvariable=scalingVariable)
    permVariable <- tclVar("100")
    permutation <- tkentry(method4Frame, width=10, textvariable=permVariable)
    xFrame <- tkframe(top, relief="groove", borderwidth=2)
    x1Frame <- tkframe(xFrame)
    x4Frame <- tkframe(xFrame)
    x2Frame <- tkframe(x4Frame)
    x3Frame <- tkframe(x4Frame)
    xBox <- tklistbox(x2Frame, width=28, height=5,
        selectmode="single", background="white", exportselection="FALSE")
    xScroll <- tkscrollbar(x2Frame, repeatinterval=5, command=function(...) tkyview(xBox, ...))
    tkconfigure(xBox, yscrollcommand=function(...) tkset(xScroll, ...))
    for (x in variables) tkinsert(xBox, "end", x)
    rhsVariable <- tclVar("")
    rhsEntry <- tkentry(x1Frame, width=60, textvariable=rhsVariable)
    plotFrame <- tkframe(top, relief="groove", borderwidth=2)
    plot1Frame <- tkframe(plotFrame)
    plot2Frame <- tkframe(plotFrame)
    plot3Frame <- tkframe(plotFrame)
    plot4Frame <- tkframe(plotFrame)
    typeBox <- tklistbox(plot1Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    typeScroll <- tkscrollbar(plot1Frame, repeatinterval=5, command=function(...) tkyview(typeBox, ...))
    tkconfigure(typeBox, yscrollcommand=function(...) tkset(typeScroll, ...))
    types <- c("plot","ordiplot","ordiplot empty","origin axes","identify sites","identify species","identify centroids","text sites","text species","text centroids",
        "points sites","points species","points centroids",
        "label sites","label species","label centroids","orditorp sites","orditorp species","orditorp centroids",
        "envfit","ordihull (factor)","ordiarrows (factor)","ordisegments (factor)","ordispider (factor)","ordiellipse (factor)","ordisurf (continuous)",
        "ordibubble (continuous)","ordisymbol (factor)","ordisymbol (factor, legend)","ordivector (species)","ordivector interpretation","ordicluster","ordicluster2",
        "ordinearest","ordispantree","ordiresids","distance displayed","coenocline", "screeplot", "orditkplot sites","orditkplot species","orditkplot pointlabel")
    for (x in types) tkinsert(typeBox, "end", x)
    choicesVariable <- tclVar("1,2")
    choice <- tkentry(plot3Frame, width=10, textvariable=choicesVariable)
    dataVariable <- tclVar("0")
    dataCheckBox <- tkcheckbutton(plot3Frame, variable=dataVariable)
    axisBox <- tklistbox(plot2Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    axisScroll <- tkscrollbar(plot2Frame, repeatinterval=5, command=function(...) tkyview(axisBox, ...))
    tkconfigure(axisBox, yscrollcommand=function(...) tkset(axisScroll, ...))
    for (x in variables) tkinsert(axisBox, "end", x)
    cexVariable <- tclVar("1")
    cexa <- tkentry(plot4Frame, width=10, textvariable=cexVariable)
    colVariable <- tclVar("blue")
    cola <- tkentry(plot4Frame, width=10, textvariable=colVariable)
    onDoubleClick <- function(){
        var <- as.character(tkget(xBox, "active"))[1]
            tkfocus(rhsEntry)
            rhs <- tclvalue(rhsVariable)
            rhs.chars <- rev(strsplit(rhs, "")[[1]])
            check.char <- if (length(rhs.chars) > 0){
                if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1)) 
                    rhs.chars[1] else rhs.chars[2]
                }
                else ""
            tclvalue(rhsVariable) <- if (rhs == "" || 
                is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%")))
                    paste(rhs, var, sep="")
                else paste(rhs, "+", var)
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onPlus <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "+ ")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onTimes <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "*", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onColon <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, ":", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onSlash <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "/",  sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onIn <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "%in% ")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onMinus <- function(){
        rhs <- tclvalue(rhsVariable)
        tclvalue(rhsVariable) <- paste(rhs, "+Condition(")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onPower <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, "^", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onLeftParen <- function(){
        tkfocus(rhsEntry)
        rhs <- tclvalue(rhsVariable)
        tclvalue(rhsVariable) <- paste(rhs, "(", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onRightParen <- function(){
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs)) return()
        tclvalue(rhsVariable) <- paste(rhs, ")", sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
        }
    onOK <- function(){
        doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
        method <- methods[as.numeric(tkcurselection(methodBox))+1]
        perm <- as.numeric(tclvalue(permVariable))
        dist <- distances[as.numeric(tkcurselection(distBox))+1]
        treatasdist <- tclvalue(treatasdistVariable)==1
        check.empty <- gsub(" ", "", tclvalue(rhsVariable))
        if ("" == check.empty) {
            tkmessageBox(message="Right-hand side of model empty.", 
                icon="error", type="ok")
            }
        formula <- paste(.communityDataSet, tclvalue(rhsVariable), sep=" ~ ")
        if (method=="RDA") {
            command <- paste("rda(", formula, ",", .activeDataSet, ")", sep="")
            doItAndPrint(paste("dist.eval(", .communityDataSet, ",'euclidean')", sep=""))
        }
        if (method=="CCA") {
            command <- paste("cca(", formula, ",", .activeDataSet, ")", sep="")
        }
        if (method=="capscale") {
            if(treatasdist==T){
                logger(paste(.communityDataSet, " <- as.dist(", .communityDataSet, ")", sep=""))
                assign(.communityDataSet, justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            command <- paste("capscale(", formula, ",", .activeDataSet, ",dist='", dist, "', add=F)", sep="")
            if(treatasdist==F){
                doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist, "')", sep=""))
            }
            doItAndPrint(paste("adonis(", formula, "," , .activeDataSet, ", method='", dist, "', permutations=", perm, ")", sep=""))
        }
        if (method=="capscale(add)") {
            if(treatasdist==T){
                logger(paste(.communityDataSet, " <- as.dist(", .communityDataSet, ")", sep=""))
                assign(.communityDataSet, justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            command <- paste("capscale(", formula, ",", .activeDataSet, ",dist='", dist, "', add=T)", sep="")
            if(treatasdist==F){
                doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist, "')", sep=""))
            }
            doItAndPrint(paste("adonis(", formula, "," , .activeDataSet, ", method='", dist, "', permutations=", perm, ")", sep=""))
        }
        if (method=="prc") {
            command <- paste("prc(", .communityDataSet, "," ,tclvalue(rhsVariable), ")", sep="")
            doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist, "')", sep=""))
        }
        if (method=="CAPdiscrim") {
            if(treatasdist==T){
                logger(paste(.communityDataSet, " <- as.dist(", .communityDataSet, ")", sep=""))
                assign(.communityDataSet, justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            command <- paste("CAPdiscrim(", formula, ",", .activeDataSet, ", dist='", dist, "', permutations=", perm,")", sep="")
            doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist, "')", sep=""))
        }
        if (method=="multiconstrained (RDA)") {
            command <- paste("multiconstrained(method='rda',", formula, ",", .activeDataSet, ", contrast=0, step=", perm, ")", sep="")
            doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist, "')", sep=""))
        }
        if (method=="multiconstrained (CCA)") {
            command <- paste("multiconstrained(method='cca',", formula, ",", .activeDataSet, ", contrast=0, step=", perm, ")", sep="")
        }
        if (method=="multiconstrained (capscale)") {
            if(treatasdist==T){
                logger(paste(.communityDataSet, " <- as.dist(", .communityDataSet, ")", sep=""))
                assign(.communityDataSet, justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            command <- paste("multiconstrained(method='capscale',", formula, ",", .activeDataSet, ",dist='", dist, "', add=F, contrast=0, step=", perm, ")", sep="")
            if(treatasdist==F){
                doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist, "')", sep=""))
            }
        }
        if (method=="multiconstrained (capscale add)") {
            if(treatasdist==T){
                logger(paste(.communityDataSet, " <- as.dist(", .communityDataSet, ")", sep=""))
                assign(.communityDataSet, justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            command <- paste("multiconstrained(method='capscale',", formula, ",", .activeDataSet, ",dist='", dist, "', add=T, contrast=0, step=", perm, ")", sep="")
            if(treatasdist==F){
                doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist, "')", sep=""))
            }
        }
        modelValue <- tclvalue(modelName)
#        if (!is.valid.name(modelValue)){
#            tkmessageBox(message=paste('"', modelValue, '" is not a valid name.', 
#                sep=""), icon="error", type="ok")
#            }        
        logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        if (method == "RDA" || method == "CCA"  || method == "capscale" || method == "capscale(add)") {
            doItAndPrint(paste("check.ordiscores(", .communityDataSet, ",", modelValue, ", check.species=T)", sep=""))
        }
        if (method == "CAPdiscrim") {
            doItAndPrint(paste("check.ordiscores(", .communityDataSet, ",", modelValue, ", check.species=F)", sep=""))
        }
        sum <- tclvalue(summaryVariable) == "1"
        scaling <- tclvalue(scalingVariable)
        if (sum==T) {
            if (method=="CAPdiscrim" || method=="multiconstrained (RDA)" || method=="multiconstrained (CCA)" || method=="multiconstrained (capscale)" || method=="multiconstrained (capscale add)") {
                doItAndPrint(paste(modelValue, sep=""))
            }else{
                doItAndPrint(paste("summary(", modelValue, ", scaling=", scaling, ")", sep=""))
            }
            if (method=="RDA" || method=="CCA"  || method=="capscale"  || method=="capscale(add)") {
                doItAndPrint(paste("RsquareAdj(", modelValue, ")", sep=""))
                doItAndPrint(paste("deviance(", modelValue, ")", sep=""))
                doItAndPrint(paste("vif.cca(", modelValue, ")", sep=""))
            }
            if (method=="RDA" || method=="CCA"  || method=="capscale"  || method=="capscale(add)") {            
                doItAndPrint(paste("goodness(", modelValue, ", display='sites', statistic='explained')", sep=""))            
                doItAndPrint(paste("inertcomp(", modelValue, ", display='sites', statistic='explained', proportional=T)", sep="")) 
            }
        }
        if (perm>0  && method !="CAPdiscrim" && method!="multiconstrained (RDA)" && method!="multiconstrained (CCA)" && method!="multiconstrained (capscale)" && method!="multiconstrained (capscale add)") {
            doItAndPrint(paste("permutest(", modelValue, ", permutations=", perm, ")", sep=""))
            doItAndPrint(paste("permutest(", modelValue, ", permutations=", perm, ", first=T)", sep=""))
            if (method !="prc") {doItAndPrint(paste("anova.cca(", modelValue, ", step=", perm, ", by='terms')", sep=""))}
            if (method !="prc") {doItAndPrint(paste("anova.cca(", modelValue, ", step=", perm, ", by='margin')", sep=""))}
            }
        if(treatasdist==T  && method!="RDA" && method!="CCA"  && method!="prc" && method!="multiconstrained (RDA)" && method!="multiconstrained (CCA)"){
            logger(paste(.communityDataSet, " <- data.frame(as.matrix(", .communityDataSet, "))", sep=""))
            assign(.communityDataSet, justDoIt(paste("data.frame(as.matrix(",.communityDataSet, "))", sep="")), envir=.GlobalEnv)
        }
    }
    onPlot <- function(){
        method <- methods[as.numeric(tkcurselection(methodBox))+1]
        modelValue <- tclvalue(modelName)
        perm <- as.numeric(tclvalue(permVariable))
        plottype <- types[as.numeric(tkcurselection(typeBox))+1]
        scaling <- tclvalue(scalingVariable)
        choices <- tclvalue(choicesVariable)
        dist <- distances[as.numeric(tkcurselection(distBox))+1]
        cex <- tclvalue(cexVariable)
        col <- tclvalue(colVariable)
        treatasdist <- tclvalue(treatasdistVariable)==1
        justDoIt(paste("attach(", .activeDataSet, ")", sep=""))
        if (plottype == "plot"  && method != "CAPdiscrim" && method!="multiconstrained (RDA)" && method!="multiconstrained (CCA)" && method!="multiconstrained (capscale)" && method!="multiconstrained (capscale add)"){
            justDoIt(paste("par(cex=",cex,")", sep=""))
            logger(paste("par(cex=",cex,")", sep=""))
            justDoIt(paste("plot1 <- plot(", modelValue, ", choices=c(", choices, "), scaling=", scaling, ")", sep=""))
            logger(paste("plot1 <- plot(", modelValue, ", choices=c(", choices, "), scaling=", scaling, ")", sep=""))
            }
        if (plottype == "ordiplot"){
            if (method!="multiconstrained (RDA)" && method!="multiconstrained (CCA)" && method!="multiconstrained (capscale)" && method!="multiconstrained (capscale add)"){
                justDoIt(paste("par(cex=",cex,")", sep=""))
                logger(paste("par(cex=",cex,")", sep=""))
                if (method == "CAPdiscrim") {
                    justDoIt(paste("plot1 <- ordiplot(", modelValue, ", choices=c(", choices, "))", sep=""))              
                    logger(paste("plot1 <- ordiplot(", modelValue, ", choices=c(", choices, "))", sep=""))              
                }else{
                    justDoIt(paste("plot1 <- ordiplot(", modelValue, ", choices=c(", choices, "), scaling=", scaling, ")", sep=""))              
                    logger(paste("plot1 <- ordiplot(", modelValue, ", choices=c(", choices, "), scaling=", scaling, ")", sep=""))              
                }
            }
        }
        if (plottype == "ordiplot empty"){
            if (method!="multiconstrained (RDA)" && method!="multiconstrained (CCA)" && method!="multiconstrained (capscale)" && method!="multiconstrained (capscale add)"){
                justDoIt(paste("par(cex=",cex,")", sep=""))
                logger(paste("par(cex=",cex,")", sep=""))
                if (method == "CAPdiscrim") {
                    justDoIt(paste("plot1 <- ordiplot(", modelValue, ", type='none', choices=c(", choices, "))", sep=""))              
                    logger(paste("plot1 <- ordiplot(", modelValue, ", type='none', choices=c(", choices, "))", sep=""))              
                }else{
                    justDoIt(paste("plot1 <- ordiplot(", modelValue, ", type='none',choices=c(", choices, "), scaling=", scaling, ")", sep="")) 
                    logger(paste("plot1 <- ordiplot(", modelValue, ", type='none',choices=c(", choices, "), scaling=", scaling, ")", sep="")) 
                }
            }
        }
        if (plottype == "identify sites"){
            doItAndPrint(paste("identify(plot1, 'sites', col='", col,"', cex=", cex, ")", sep=""))             
            }
        if (plottype == "identify species"){
            doItAndPrint(paste("identify(plot1, 'species', col='", col,"', cex=", cex, ")", sep=""))             
            }
        if (plottype == "identify centroids"){
            doItAndPrint(paste("identify(plot1, 'centroids', col='", col,"', cex=", cex, ")",sep=""))             
            }
        if (plottype == "text sites"){
            doItAndPrint(paste("text(plot1, 'sites', col='", col,"', cex=", cex, ")",sep=""))             
            }
        if (plottype == "text species"){
            doItAndPrint(paste("text(plot1, 'species', col='", col,"', cex=", cex, ")",sep=""))             
            }
        if (plottype == "text centroids"){
            doItAndPrint(paste("text(plot1, 'centroids', col='", col,"', cex=", cex, ")",sep=""))             
            }
        if (plottype == "points sites"){
            doItAndPrint(paste("points(plot1, 'sites', col='", col,"', cex=", cex, ")",sep=""))             
        }
        if (plottype == "points species"){
            doItAndPrint(paste("points(plot1, 'species', col='", col,"', cex=", cex, ")",sep=""))             
        }
        if (plottype == "points centroids"){
            doItAndPrint(paste("points(plot1, 'centroids', col='", col,"', cex=", cex, ")",sep=""))             
            }
        if (plottype == "label sites"){
            doItAndPrint(paste("ordilabel(plot1, 'sites', col='", col,"', cex=", cex, ")",sep=""))             
            }
        if (plottype == "label species"){
            doItAndPrint(paste("ordilabel(plot1, 'species', col='", col,"', cex=", cex, ")",sep=""))             
            }
        if (plottype == "label centroids"){
            doItAndPrint(paste("ordilabel(plot1, 'centroids', col='", col,"', cex=", cex, ")",sep=""))             
            }
        if (plottype == "orditorp sites"){
            doItAndPrint(paste("orditorp(plot1, 'sites', col='", col,"', cex=", cex, ")",sep=""))             
            }
        if (plottype == "orditorp species"){
            doItAndPrint(paste("orditorp(plot1, 'species', col='", col,"', cex=", cex, ")",sep=""))             
            }
        if (plottype == "orditorp centroids"){
            doItAndPrint(paste("orditorp(plot1, 'centroids', col='", col,"', cex=", cex, ")",sep=""))             
            }
        if (plottype == "origin axes"){
            doItAndPrint(paste("abline(h = 0, lty = 3)", sep=""))
            doItAndPrint(paste("abline(v = 0, lty = 3)", sep=""))
        }
        if (plottype == "ordiresids"){
            if (method != "CAPdiscrim"  && method!="multiconstrained (RDA)" && method!="multiconstrained (CCA)" && method!="multiconstrained (capscale)" && method!="multiconstrained (capscale add)"){
                justDoIt(paste("par(cex=",cex,")", sep=""))
                logger(paste("par(cex=",cex,")", sep=""))
                justDoIt(paste("ordiresids(", modelValue, ", kind='residuals')", sep=""))              
                logger(paste("ordiresids(", modelValue, ", kind='residuals')", sep=""))              
            }
        }
        if (plottype == "orditkplot sites"){
            if (method!="multiconstrained (RDA)" && method!="multiconstrained (CCA)" && method!="multiconstrained (capscale)" && method!="multiconstrained (capscale add)"){
                justDoIt(paste("par(cex=",cex,")", sep=""))
                logger(paste("par(cex=",cex,")", sep=""))
                if (method != "CAPdiscrim") {
                    justDoIt(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", display='sites', choices=c(", choices, "), scaling=", scaling, "))", sep=""))
                    logger(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", display='sites', choices=c(", choices, "), scaling=", scaling, "))", sep=""))
                }else{
                    justDoIt(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", display='sites', choices=c(", choices, ")))", sep="")) 
                    logger(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", display='sites', choices=c(", choices, ")))", sep="")) 
                }
            }
        }
        if (plottype == "orditkplot species"){
            if (method!="multiconstrained (RDA)" && method!="multiconstrained (CCA)" && method!="multiconstrained (capscale)" && method!="multiconstrained (capscale add)"){
                justDoIt(paste("par(cex=",cex,")", sep=""))
                logger(paste("par(cex=",cex,")", sep=""))
                if (method != "CAPdiscrim") {
                    justDoIt(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", display='species', choices=c(", choices, "), scaling=", scaling, "))", sep=""))
                    logger(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", display='species', choices=c(", choices, "), scaling=", scaling, "))", sep=""))
                }else{
                    justDoIt(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", display='species', choices=c(", choices, ")))", sep="")) 
                    logger(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", display='species', choices=c(", choices, ")))", sep="")) 
                }
            }
        }
        if (plottype == "orditkplot pointlabel"){
            if (method!="multiconstrained (RDA)" && method!="multiconstrained (CCA)" && method!="multiconstrained (capscale)" && method!="multiconstrained (capscale add)"){
                justDoIt(paste("par(cex=",cex,")", sep=""))
                logger(paste("par(cex=",cex,")", sep=""))
                if (method != "CAPdiscrim") {
                    justDoIt(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", choices=c(", choices, "), scaling=", scaling, "))", sep=""))
                    logger(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", choices=c(", choices, "), scaling=", scaling, "))", sep=""))
                }else{
                    justDoIt(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", choices=c(", choices, ")))", sep="")) 
                    logger(paste("plot1 <- orditkplot(ordipointlabel(", modelValue, ", choices=c(", choices, ")))", sep="")) 
                }
            }
        }
        axisvar <- .variables[as.numeric(tkcurselection(axisBox))+1]
        varfactor <- eval(parse(text=paste("is.factor(",.activeDataSet, "$", axisvar, ")", sep="")), envir=.GlobalEnv)
        if (plottype == "envfit"){
            doItAndPrint(paste("fitted <- envfit(plot1, data.frame(", axisvar, "), permutations=", perm, ")", sep=""))
            doItAndPrint(paste("plot(fitted, col='", col,"', cex=", cex, ")", sep=""))
        }
        if (plottype == "ordihull (factor)" && varfactor==T){
            doItAndPrint(paste("ordihull(plot1,", axisvar, ", col='", col, "')", sep=""))
            }
        if (plottype == "ordiarrows (factor)" && varfactor==T){
            doItAndPrint(paste("ordiarrows(plot1,", axisvar, ", col='", col, "')", sep=""))
            }
        if (plottype == "ordisegments (factor)" && varfactor==T){
            doItAndPrint(paste("ordisegments(plot1,", axisvar, ", col='", col, "')", sep=""))
            }
        if (plottype == "ordispider (factor)" && varfactor==T){
            doItAndPrint(paste("ordispider(plot1,", axisvar, ", col='", col, "')", sep=""))
            }
        if (plottype == "ordiellipse (factor)" && varfactor==T){
            doItAndPrint(paste("ordiellipse(plot1,", axisvar, ", col='", col, "', conf=0.9)", sep=""))
            }
        if (plottype == "ordisurf (continuous)" && varfactor==F){
            doItAndPrint(paste("ordisurf(plot1,", axisvar, ", add=T)", sep=""))
            }
        if (plottype == "ordibubble (continuous)" && varfactor==F){
            doItAndPrint(paste("ordibubble(plot1,", axisvar, ", fg='", col, "')", sep=""))
            } 
        if (plottype == "ordisymbol (factor)" && varfactor==T){
            justDoIt(paste("ordisymbol(plot1,", .activeDataSet, ",'", axisvar, "', rainbow=T, legend=F, cex=", cex, ")", sep=""))
            logger(paste("ordisymbol(plot1,", .activeDataSet, ",'", axisvar, "', rainbow=T, legend=F, cex=", cex, ")", sep=""))
            }
        if (plottype == "ordisymbol (factor, legend)" && varfactor==T){
            justDoIt(paste("ordisymbol(plot1,", .activeDataSet, ",'", axisvar, "', rainbow=T, legend=T, legend.x='topleft', legend.ncol=1, cex=", cex, ")", sep=""))
            logger(paste("ordisymbol(plot1,", .activeDataSet, ",'", axisvar, "', rainbow=T, legend=T, legend.x='topleft', legend.ncol=1, cex=", cex, ")", sep=""))
            }
        if (plottype == "ordivector (species)"){
            realspecies <- eval(parse(text=paste("any(colnames(", .communityDataSet, ")=='", axisvar, "')", sep="")), envir=.GlobalEnv)
            if (realspecies == T) {
                doItAndPrint(paste("ordivector(plot1,'", axisvar, "',lty=0, angle=5, length=0.5)", sep=""))
                }
            }
        if (plottype == "ordivector interpretation"){
            realspecies <- eval(parse(text=paste("any(colnames(", .communityDataSet, ")=='", axisvar, "')", sep="")), envir=.GlobalEnv)
            if (realspecies == T) {
                doItAndPrint(paste("ordivector(plot1,'", axisvar, "',lty=2)", sep=""))
            }
        }         
        if (plottype == "ordicluster"){
            if(treatasdist==F){
                logger(paste("distmatrix <- vegdist(", .communityDataSet, ",method='", dist, "', na.rm=T)", sep=""))
                assign("distmatrix", justDoIt(paste("vegdist(",.communityDataSet, ",method='",dist, "', na.rm=T)", sep="")), envir=.GlobalEnv)
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            logger(paste("cluster <- hclust(distmatrix, method='single')", sep=""))
            assign("cluster", justDoIt(paste("hclust(distmatrix, method='single')", sep="")), envir=.GlobalEnv)
            doItAndPrint(paste("ordicluster(plot1, cluster, prune=0,col='", col, "')", sep=""))
            }
        if (plottype == "ordicluster2"){
            if(treatasdist==F){
                logger(paste("distmatrix <- vegdist(", .communityDataSet, ",method='", dist, "', na.rm=T)", sep=""))
                assign("distmatrix", justDoIt(paste("vegdist(",.communityDataSet, ",method='",dist, "', na.rm=T)", sep="")), envir=.GlobalEnv)
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            logger(paste("cluster <- hclust(distmatrix, method='single')", sep=""))
            assign("cluster", justDoIt(paste("hclust(distmatrix, method='single')", sep="")), envir=.GlobalEnv)
            doItAndPrint(paste("ordicluster2(plot1, cluster, mingroups=1, col='", col, "')", sep=""))
            }
        if (plottype == "ordinearest"){
            if(treatasdist==F){
                logger(paste("distmatrix <- vegdist(", .communityDataSet, ",method='", dist, "', na.rm=T)", sep=""))
                assign("distmatrix", justDoIt(paste("vegdist(",.communityDataSet, ",method='",dist, "', na.rm=T)", sep="")), envir=.GlobalEnv)
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            doItAndPrint(paste("ordinearest(plot1, distmatrix,col='", col, "')", sep=""))
            }
        if (plottype == "ordispantree"){
            if(treatasdist==F){
                logger(paste("distmatrix <- vegdist(", .communityDataSet, ",method='", dist, "', na.rm=T)", sep=""))
                assign("distmatrix", justDoIt(paste("vegdist(",.communityDataSet, ",method='",dist, "', na.rm=T)", sep="")), envir=.GlobalEnv)
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            doItAndPrint(paste("lines(spantree(distmatrix,toolong=0),plot1,col='", col, "')", sep=""))
            }
        if (plottype == "distance displayed"){
            if(treatasdist==F){
                doItAndPrint(paste("distdisplayed(", .communityDataSet ,",plot1, distx='", dist, "',plotit=T)", sep=""))
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
                doItAndPrint(paste("distdisplayed(distmatrix, plot1, plotit=T)", sep=""))
            }
            }
        if (plottype == "coenocline"){
            doItAndPrint(paste("ordicoeno(", .communityDataSet, ", ordiplot=plot1, axis=1, legend=T, cex=0.8, ncol=4)", sep=""))
            }
        if (plottype == "screeplot"){
            justDoIt(paste("par(cex=",cex,")", sep=""))
            logger(paste("par(cex=",cex,")", sep=""))
            justDoIt(paste("plot1 <- screeplot(", modelValue, ")", sep=""))
            logger(paste("plot1 <- screeplot(", modelValue, ")", sep=""))
        }              
        data <- tclvalue(dataVariable) =="1"
        if (data==T) {
            justDoIt(paste(.activeDataSet, "$", modelValue, ".ax1 <- scores(plot1,display='sites')[,1]", sep=""))
            logger(paste(.activeDataSet, "$", modelValue, ".ax1 <- scores(plot1,display='sites')[,1]", sep=""))
            justDoIt(paste(.activeDataSet, "$", modelValue, ".ax2 <- scores(plot1,display='sites')[,2]", sep=""))
            logger(paste(.activeDataSet, "$", modelValue, ".ax2 <- scores(plot1,display='sites')[,2]", sep=""))
            activeDataSet(.activeDataSet)
        }
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    .operatorFont <- getRcmdr("operatorFont")
    plusButton <- tkbutton(x3Frame, text="+", width="3", command=onPlus, 
        font=.operatorFont)
    timesButton <- tkbutton(x3Frame, text="*", width="3", command=onTimes, 
        font=.operatorFont)
    colonButton <- tkbutton(x3Frame, text=":", width="3", command=onColon, 
        font=.operatorFont)
    slashButton <- tkbutton(x3Frame, text="/", width="3", command=onSlash, 
        font=.operatorFont)
    inButton <- tkbutton(xFrame, text="%in%", width="3", command=onIn,
        font=.operatorFont)
    minusButton <- tkbutton(x3Frame, text="Cond", width="3", command=onMinus, 
        font=.operatorFont)
    powerButton <- tkbutton(x3Frame, text="^", width="3", command=onPower, 
        font=.operatorFont)
    leftParenButton <- tkbutton(x3Frame, text="(", width="3", command=onLeftParen, 
        font=.operatorFont)
    rightParenButton <- tkbutton(x3Frame, text=")", width="3", command=onRightParen, 
        font=.operatorFont)
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    plotButton <- tkbutton(buttonsFrame, text="Plot", width="12", command=onPlot)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(modelFrame, text="Save model as:               ", width=20), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(tklabel(method1Frame, text="Ordination method"), sticky="w")
    tkgrid(methodBox, methodScroll,sticky="w")
    tkgrid(tklabel(method2Frame, text="Distance"), sticky="w")
    tkgrid(distBox, distScroll,sticky="w")
    tkgrid(summaryCheckBox, tklabel(method3Frame, text="model summary"), sticky="w")
    tkgrid(treatasdistCheckBox, tklabel(method3Frame, text="as.dist(Community)", width=15), sticky="w")
    tkgrid(tklabel(method4Frame, text="scaling", width=10), scale, sticky="w")
    tkgrid(tklabel(method4Frame, text="permutations", width=10), permutation, sticky="w")
    tkgrid(method1Frame, tklabel(methodFrame, text="", width=1), method2Frame, sticky="w")
    tkgrid(method3Frame, tklabel(methodFrame, text="", width=1), method4Frame, sticky="w")
    tkgrid(methodFrame, sticky="w")
    tkgrid(rhsEntry, sticky="w")
    tkgrid(xBox, xScroll,sticky="w") 
    tkgrid(plusButton, timesButton, colonButton, slashButton, inButton, sticky="w")
    tkgrid(minusButton,powerButton, leftParenButton, rightParenButton, sticky="w")
    tkgrid(tklabel(xFrame, text="Explanatory"), sticky="w")
    tkgrid(x1Frame, sticky="w")
    tkgrid(x2Frame, tklabel(xFrame, text="", width=1), x3Frame, sticky="w")
    tkgrid(x4Frame, sticky="w")
    tkgrid(xFrame, sticky="w")
    tkgrid(tklabel(plot1Frame, text="Plot method"), sticky="w")
    tkgrid(typeBox, typeScroll, sticky="nw")
    tkgrid(tklabel(plot2Frame, text="Plot variable"), sticky="w")
    tkgrid(axisBox, axisScroll, sticky="nw")
    tkgrid(tklabel(plot3Frame, text="axes", width=10), choice, sticky="w")
    tkgrid(dataCheckBox, tklabel(plot3Frame, text="add scores to dataframe"), sticky="w")
    tkgrid(tklabel(plot4Frame, text="cex", width=10), cexa, sticky="w")
    tkgrid(tklabel(plot4Frame, text="colour", width=10), cola, sticky="w")
    tkgrid(plot1Frame, tklabel(plotFrame, text="", width=1), plot2Frame, sticky="w")
    tkgrid(plot3Frame, tklabel(plotFrame, text="", width=1), plot4Frame, sticky="w")
    tkgrid(plotFrame, sticky="w")
    tkgrid(OKbutton, plotButton, cancelButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(xScroll, sticky="ns")
    tkgrid.configure(typeScroll, sticky="ns")
    tkgrid.configure(axisScroll, sticky="ns")
    tkgrid.configure(methodScroll, sticky="ns")
    tkselection.set(xBox, 0)
    tkselection.set(typeBox, 0)
    tkselection.set(distBox, 0)
    tkselection.set(methodBox, 0)
    tkselection.set(axisBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkbind(xBox, "<Double-ButtonPress-1>", onDoubleClick)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(methodBox)
    tkwait.window(top)
    }


clusterGUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Cluster analysis")
    .activeDataSet <- ActiveDataSet()
    .communityDataSet <- CommunityDataSet()
    .variables <- Variables()
    variables <- paste(.variables, ifelse(is.element(.variables, Factors()), "[factor]", ""))
    modelName <- tclVar("Cluster.1")
    modelFrame <- tkframe(top, relief="groove", borderwidth=2)
    model <- tkentry(modelFrame, width=40, textvariable=modelName)
    methodFrame <- tkframe(top, relief="groove", borderwidth=2)
    method1Frame <- tkframe(methodFrame)
    method2Frame <- tkframe(methodFrame)
    method3Frame <- tkframe(methodFrame)
    method4Frame <- tkframe(methodFrame)
    methodBox <- tklistbox(method1Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    methodScroll <- tkscrollbar(method1Frame, repeatinterval=5, command=function(...) tkyview(methodBox, ...))
    tkconfigure(methodBox, yscrollcommand=function(...) tkset(methodScroll, ...))
    methods <- c("hclust","agnes","diana","kmeans","cascadeKM","pam","clara","fanny")
    for (x in methods) tkinsert(methodBox, "end", x)
    distBox <- tklistbox(method2Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    distScroll <- tkscrollbar(method2Frame, repeatinterval=5, command=function(...) tkyview(distBox, ...))
    tkconfigure(distBox, yscrollcommand=function(...) tkset(distScroll, ...))
    distances <- c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", 
        "chao", "cao", "mahalanobis")
    for (x in distances) tkinsert(distBox, "end", x)
    treatasdistVariable <- tclVar("0")
    treatasdistCheckBox <- tkcheckbutton(method3Frame, variable=treatasdistVariable)
    summaryVariable <- tclVar("1")
    summaryCheckBox <- tkcheckbutton(method3Frame, variable=summaryVariable)
    copheneticVariable <- tclVar("0")
    copheneticCheckBox <- tkcheckbutton(method3Frame, variable=copheneticVariable)
    clustersVariable <- tclVar("5")
    clustersa <- tkentry(method3Frame, width=10, textvariable=clustersVariable)
    dataVariable <- tclVar("0")
    dataCheckBox <- tkcheckbutton(method3Frame, variable=dataVariable)
    algoBox <- tklistbox(method4Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    algoScroll <- tkscrollbar(method4Frame, repeatinterval=5, command=function(...) tkyview(algoBox, ...))
    tkconfigure(algoBox, yscrollcommand=function(...) tkset(algoScroll, ...))
    algos <- c("average","single","complete","ward","weighted","median","centroid")
    for (x in algos) tkinsert(algoBox, "end", x)
    plotFrame <- tkframe(top, relief="groove", borderwidth=2)
    plot1Frame <- tkframe(plotFrame)
    plot2Frame <- tkframe(plotFrame)
    plot3Frame <- tkframe(plotFrame)
    plot4Frame <- tkframe(plotFrame)
    typeBox <- tklistbox(plot1Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    typeScroll <- tkscrollbar(plot1Frame, repeatinterval=5, command=function(...) tkyview(typeBox, ...))
    tkconfigure(typeBox, yscrollcommand=function(...) tkset(typeScroll, ...))
    types <- c("dendrogram1 (hang = -1)","dendrogram2 (hang = 0.1)", "dendrogram3 (horizontal)", 
        "phylogram (ape package)", "cladogram (ape package)", "fan (ape package)", "unrooted (ape package)",    
        "rectangles","pruned dendrogram","silhouette","kgs","cophenetic","cascadeKM","reorder (variable)", "tiplabels (variable size)", "tiplabels (factor)")
    for (x in types) tkinsert(typeBox, "end", x)
    cexVariable <- tclVar("1")
    cexa <- tkentry(plot3Frame, width=8, textvariable=cexVariable)
    colVariable <- tclVar("blue")
    cola <- tkentry(plot4Frame, width=8, textvariable=colVariable)
    axisBox <- tklistbox(plot2Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    axisScroll <- tkscrollbar(plot2Frame, repeatinterval=5, command=function(...) tkyview(axisBox, ...))
    tkconfigure(axisBox, yscrollcommand=function(...) tkset(axisScroll, ...))
    for (x in variables) tkinsert(axisBox, "end", x)
    onOK <- function(){
        doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
        method <- methods[as.numeric(tkcurselection(methodBox))+1]
        dist <- distances[as.numeric(tkcurselection(distBox))+1]
        algo <- algos[as.numeric(tkcurselection(algoBox))+1]
        treatasdist <- tclvalue(treatasdistVariable)==1
        clusters <- tclvalue(clustersVariable)
        axisvar <- .variables[as.numeric(tkcurselection(axisBox))+1]
        if (method=="agnes" || method=="diana" || method=="pam" || method=="clara" || method=="fanny") {
            justDoIt(paste("library(cluster)"))
            logger(paste("library(cluster)"))
        }
        if (method != "kmeans"  && method != "cascadeKM") {
            if(treatasdist==F){
                logger(paste("distmatrix <- vegdist(", .communityDataSet, ",method='", dist, "', na.rm=T)", sep=""))
                assign("distmatrix", justDoIt(paste("vegdist(",.communityDataSet, ",method='",dist, "', na.rm=T)", sep="")), envir=.GlobalEnv)
                doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist, "')", sep=""))
            }else{
                logger(paste("distmatrix <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
        }
        if (method=="hclust") {
            command <- paste("hclust(distmatrix, method='", algo, "')", sep="")
        }
        if (method=="agnes") {
            command <- paste("agnes(distmatrix, method='", algo, "')", sep="")
        }
        if (method=="diana") {
            command <- paste("diana(distmatrix)", sep="")
        }
        if (method=="kmeans") {
            command <- paste("kmeans(", .communityDataSet, ", centers=", clusters, ", iter.max=100)", sep="")
        }
        if (method=="cascadeKM") {
            command <- paste("cascadeKM(", .communityDataSet, ", inf.gr=2, sup.gr=", clusters, ", iter = 100, criterion = 'calinski')", sep="")
        }
        if (method=="pam") {
            command <- paste("pam(distmatrix, k=", clusters, ")", sep="")
        }
        if (method=="clara") {
            size <- as.numeric(clusters)+1
            command <- paste("clara(distmatrix, k=", clusters, ", sampsize=", size,")", sep="")
        }
        if (method=="fanny") {
            command <- paste("fanny(distmatrix, k=", clusters, ")", sep="")
        }
        modelValue <- tclvalue(modelName)      
        logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        sum <- tclvalue(summaryVariable) == "1"
        if (sum==T) {
            if (method=="kmeans" || method=="hclust" || method=="cascadeKM") {
                doItAndPrint(paste(modelValue))
                doItAndPrint(paste("attributes(", modelValue, ")", sep=""))
                if (method=="hclust") {doItAndPrint(paste("treeheight(", modelValue, ")", sep=""))}
            }else{
                doItAndPrint(paste("summary(", modelValue, ")", sep=""))
            }
        }
        coph <- tclvalue(copheneticVariable) == "1"
        if (coph==T && method != "kmeans" && method != "pam" && method != "clara" && method != "fanny") {
            logger(paste("copheneticdist <- cophenetic(", modelValue, ")", sep=""))
            assign("copheneticdist", justDoIt(paste("cophenetic(", modelValue, ")", sep="")), envir=.GlobalEnv)
            doItAndPrint(paste("mantel(distmatrix,copheneticdist,permutations=100)",sep=""))
        }
        data <- tclvalue(dataVariable) == "1"
        if (data==T && method!="cascadeKM") {
            if (method =="kmeans" || method== "pam" || method=="clara" || method=="fanny") {
                justDoIt(paste(.activeDataSet, "$", modelValue, ".cluster <- as.factor(", modelValue, "$cluster)", sep=""))
                logger(paste(.activeDataSet, "$", modelValue, ".cluster <- as.factor(", modelValue, "$cluster)", sep=""))
            }else{
                justDoIt(paste(.activeDataSet, "$", modelValue, ".cluster <- as.factor(cutree(", modelValue, ", k=",clusters, "))", sep=""))
                logger(paste(.activeDataSet, "$", modelValue, ".cluster <- as.factor(cutree(", modelValue, ", k=",clusters, "))", sep=""))
            }
            activeDataSet(.activeDataSet)
        }
    }
    onPlot <- function(){
        method <- methods[as.numeric(tkcurselection(methodBox))+1]
        modelValue <- tclvalue(modelName)
        plottype <- types[as.numeric(tkcurselection(typeBox))+1]
        dist <- distances[as.numeric(tkcurselection(distBox))+1]
        clusters <- tclvalue(clustersVariable)
        axisvar <- .variables[as.numeric(tkcurselection(axisBox))+1]
        col <- tclvalue(colVariable)
        cex <- tclvalue(cexVariable)
        justDoIt(paste("par(cex=",cex,")", sep=""))
        logger(paste("par(cex=",cex,")", sep=""))
        if (plottype == "dendrogram1 (hang = -1)"){
            doItAndPrint(paste("plot(as.dendrogram(", modelValue, ", hang = -1), horiz=F, edgePar=list(col='", col, "'), nodePar=list(pch=NA, lab.col='", col, "'), main='', sub='', xlab='', ylab='')", sep=""))              
        }
        if (plottype == "dendrogram2 (hang = 0.1)"){
            doItAndPrint(paste("plot(as.dendrogram(", modelValue, ", hang = 0.1), horiz=F, edgePar=list(col='", col, "'), nodePar=list(pch=NA, lab.col='", col, "'), main='', sub='', xlab='', ylab='')", sep=""))              
        }
        if (plottype == "dendrogram3 (horizontal)"){
            doItAndPrint(paste("plot(as.dendrogram(", modelValue, "), horiz=T, edgePar=list(col='", col, "'), nodePar=list(pch=NA, lab.col='", col, "'), main='', sub='', xlab='', ylab='')", sep=""))              
        }
        if (plottype == "phylogram (ape package)"){
            justDoIt(paste("library(ape)", sep=""))
            logger(paste("library(ape)", sep=""))
            doItAndPrint(paste("plot(as.phylo(as.hclust(", modelValue, ")), type='phylogram', direction='rightwards', edge.color='", col, "', tip.color='", col, "', font=1)", sep=""))              
        }
        if (plottype == "cladogram (ape package)"){
            justDoIt(paste("library(ape)", sep=""))
            logger(paste("library(ape)", sep=""))
            doItAndPrint(paste("plot(as.phylo(as.hclust(", modelValue, ")), type='cladogram', edge.color='", col, "', tip.color='", col, "', font=1)", sep=""))              
        }
        if (plottype == "fan (ape package)"){
            justDoIt(paste("library(ape)", sep=""))
            logger(paste("library(ape)", sep=""))
            doItAndPrint(paste("plot(as.phylo(as.hclust(", modelValue, ")), type='fan', edge.color='", col, "', tip.color='", col, "', font=1)", sep=""))              
        }
        if (plottype == "unrooted (ape package)"){
            justDoIt(paste("library(ape)", sep=""))
            logger(paste("library(ape)", sep=""))
            doItAndPrint(paste("plot(as.phylo(as.hclust(", modelValue, ")), type='unrooted', edge.color='", col, "', tip.color='", col, "', font=1)", sep=""))              
        }
        if (plottype == "pruned dendrogram"){
            justDoIt(paste("library(maptree)", sep=""))
            logger(paste("library(maptree)", sep=""))
            doItAndPrint(paste("plot(clip.clust(as.hclust(", modelValue,"), data=", .communityDataSet, ", k=", clusters, "))", sep=""))              
        }
        if (plottype == "silhouette"){
            doItAndPrint(paste("plot(silhouette(cutree(as.hclust(", modelValue,"), k=", clusters, "), distmatrix))", sep=""))
        }  
        if (plottype == "kgs" && method != "kmeans" && method != "cascadeKM" && method != "pam" && method != "clara" && method != "fanny"){
            justDoIt(paste("library(maptree)", sep=""))
            logger(paste("library(maptree)", sep=""))
            doItAndPrint(paste("plot(kgs(", modelValue, ", distmatrix, maxclust=min(20, nrow(", .communityDataSet, "))))", sep=""))             
        }        
        if (plottype == "rectangles" && method != "kmeans" && method != "cascadeKM" && method != "pam" && method != "clara" && method != "fanny"){
            justDoIt(paste("rect.hclust(", modelValue, ", k=", clusters, ", border='", col, "')", sep=""))             
            logger(paste("rect.hclust(", modelValue, ", k=", clusters, ", border='", col, "')", sep=""))  
        } 
        if (plottype == "cophenetic" && method != "kmeans" && method != "cascadeKM" && method != "pam" && method != "clara" && method != "fanny"){
            logger(paste("copheneticdist <- cophenetic(", modelValue, ")", sep=""))
            assign("copheneticdist", justDoIt(paste("cophenetic(", modelValue, ")", sep="")), envir=.GlobalEnv)
            doItAndPrint(paste("plot(distmatrix, copheneticdist, col='", col, "')", sep=""))
            doItAndPrint(paste("abline(0,1)", sep="")) 
        }
        if (plottype == "cascadeKM"){
            doItAndPrint(paste("plot(", modelValue, ")", sep=""))             
        }
        if (plottype == "reorder (variable)"){
            command <- paste("reorder(as.hclust(", modelValue, "), wts=as.numeric(", .activeDataSet, "$", axisvar, "))", sep="")
            logger(paste(modelValue, " <- ", command, sep=""))
            assign(modelValue, justDoIt(command), envir=.GlobalEnv)
            doItAndPrint(paste("plot(as.dendrogram(", modelValue, ", hang = 0.1), horiz=F, edgePar=list(col='", col, "'), nodePar=list(pch=NA, lab.col='", col, "'), main='', sub='', xlab='', ylab='')", sep=""))              
        }
        if (plottype == "tiplabels (variable size)"){
            justDoIt(paste("library(ape)", sep=""))
            logger(paste("library(ape)", sep=""))
            doItAndPrint(paste("plot(as.phylo(as.hclust(", modelValue, ")), type='phylogram', direction='rightwards', edge.color='", col, "', tip.color='", col, "', font=1, label.offset=0.05)", sep=""))
            doItAndPrint(paste("tiplabels(pch=19, col='", col, "', cex=3*as.numeric(", .activeDataSet, "$", axisvar, ")/max(as.numeric(", .activeDataSet, "$", axisvar, ")))", sep=""))
        }
        if (plottype == "tiplabels (factor)"){
            justDoIt(paste("library(ape)", sep=""))
            logger(paste("library(ape)", sep=""))
            doItAndPrint(paste("plot(as.phylo(as.hclust(", modelValue, ")), type='phylogram', direction='rightwards', edge.color='", col, "', tip.color='", col, "', font=1, label.offset=0.05)", sep=""))
            doItAndPrint(paste("tiplabels(pch=19, cex=2, col=as.numeric(", .activeDataSet, "$", axisvar, "))", sep=""))
        }
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
    }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    plotButton <- tkbutton(buttonsFrame, text="Plot", width="12", command=onPlot)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(modelFrame, text="Save cluster as:               ", width=20), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(tklabel(method1Frame, text="Cluster method"), sticky="w")
    tkgrid(methodBox, methodScroll,sticky="w")
    tkgrid(tklabel(method2Frame, text="Distance"), sticky="w")
    tkgrid(distBox, distScroll,sticky="w")
    tkgrid(summaryCheckBox, tklabel(method3Frame, text="cluster summary"), sticky="w")
    tkgrid(treatasdistCheckBox, tklabel(method3Frame, text="as.dist(Community)", width=15), sticky="w")
    tkgrid(copheneticCheckBox, tklabel(method3Frame, text="cophenetic correlation"), sticky="w")
    tkgrid(tklabel(method3Frame, text="clusters", width=8), clustersa, sticky="w")
    tkgrid(dataCheckBox, tklabel(method3Frame, text="save cluster membership"), sticky="w")
    tkgrid(tklabel(method4Frame, text="Cluster options"), sticky="w")
    tkgrid(algoBox, algoScroll, sticky="nw")
    tkgrid(method1Frame, tklabel(methodFrame, text="", width=1), method2Frame, sticky="w")
    tkgrid(method3Frame, tklabel(methodFrame, text="", width=1), method4Frame, sticky="w")
    tkgrid(methodFrame, sticky="w")
    tkgrid(tklabel(plot1Frame, text="Plot options"), sticky="w")
    tkgrid(typeBox, typeScroll, sticky="nw")
    tkgrid(tklabel(plot2Frame, text="Plot variable"), sticky="w")
    tkgrid(axisBox, axisScroll, sticky="nw")
    tkgrid(tklabel(plot3Frame, text="cex", width=5), cexa, sticky="w")
    tkgrid(tklabel(plot4Frame, text="colour", width=5), cola, sticky="w")
    tkgrid(plot1Frame, tklabel(plotFrame, text="", width=1), plot2Frame, sticky="w")
    tkgrid(plot3Frame, tklabel(plotFrame, text="", width=1), plot4Frame, sticky="w")
    tkgrid(plotFrame, sticky="w")
    tkgrid(OKbutton, plotButton, cancelButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(typeScroll, sticky="ns")
    tkgrid.configure(algoScroll, sticky="ns")
    tkgrid.configure(methodScroll, sticky="ns")
    tkgrid.configure(distScroll, sticky="ns")
    tkgrid.configure(axisScroll, sticky="ns")
    tkselection.set(typeBox, 0)
    tkselection.set(methodBox, 0)
    tkselection.set(algoBox, 0)
    tkselection.set(distBox, 0)
    tkselection.set(axisBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(methodBox)
    tkwait.window(top)
    }


mantelGUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Compare distance matrices")
    .activeDataSet <- ActiveDataSet()
    .communityDataSet <- CommunityDataSet()
    .variables <- Variables()
    variables <- paste(.variables, ifelse(is.element(.variables, Factors()), "[factor]", ""))
    methodFrame <- tkframe(top, relief="groove", borderwidth=2)
    method1Frame <- tkframe(methodFrame)
    method2Frame <- tkframe(methodFrame)
    method3Frame <- tkframe(methodFrame)
    method4Frame <- tkframe(methodFrame)
    method5Frame <- tkframe(methodFrame)
    method6Frame <- tkframe(methodFrame)
    testBox <- tklistbox(method1Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    testScroll <- tkscrollbar(method1Frame, repeatinterval=5, command=function(...) tkyview(testBox, ...))
    tkconfigure(testBox, yscrollcommand=function(...) tkset(testScroll, ...))
    tests <- c("mantel","anosim (factor)","mrpp (factor)","rankindex", "bioenv", "betadisper (factor)", "meandist (factor)", "simper (factor)")
    for (x in tests) tkinsert(testBox, "end", x)
    dist1Box <- tklistbox(method3Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    dist1Scroll <- tkscrollbar(method3Frame, repeatinterval=5, command=function(...) tkyview(dist1Box, ...))
    tkconfigure(dist1Box, yscrollcommand=function(...) tkset(dist1Scroll, ...))
    distances <- c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", 
        "chao", "cao", "mahalanobis")
    for (x in distances) tkinsert(dist1Box, "end", x)
    dist2Box <- tklistbox(method2Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    dist2Scroll <- tkscrollbar(method2Frame, repeatinterval=5, command=function(...) tkyview(dist2Box, ...))
    tkconfigure(dist2Box, yscrollcommand=function(...) tkset(dist2Scroll, ...))
    distances2 <- c("daisy (factor)",distances)
    for (x in distances2) tkinsert(dist2Box, "end", x)   
    scaleBox <- tklistbox(method4Frame, width=27, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    scaleScroll <- tkscrollbar(method4Frame, repeatinterval=5, command=function(...) tkyview(scaleBox, ...))
    tkconfigure(scaleBox, yscrollcommand=function(...) tkset(scaleScroll, ...))
    svariables <- c("all",variables)
    for (x in svariables) tkinsert(scaleBox, "end", x)
    treatasdistVariable <- tclVar("0")
    treatasdistCheckBox <- tkcheckbutton(method5Frame, variable=treatasdistVariable)
    plotVariable <- tclVar("0")
    plotCheckBox <- tkcheckbutton(method5Frame, variable=plotVariable)
    permVariable <- tclVar("100")
    perma <- tkentry(method5Frame, width=10, textvariable=permVariable)
    methodBox <- tklistbox(method6Frame, width=27, height=3,
        selectmode="single", background="white", exportselection="FALSE") 
    methodScroll <- tkscrollbar(method6Frame, repeatinterval=3, command=function(...) tkyview(methodBox, ...))
    tkconfigure(methodBox, yscrollcommand=function(...) tkset(methodScroll, ...))
    methods <- c("pearson","spearman","kendall")
    for (x in methods) tkinsert(methodBox, "end", x)
    onOK <- function(){
        test <- tests[as.numeric(tkcurselection(testBox))+1]
        dist1 <- distances[as.numeric(tkcurselection(dist1Box))+1]
        dist2 <- distances2[as.numeric(tkcurselection(dist2Box))+1]
        method <- methods[as.numeric(tkcurselection(methodBox))+1]
        permutations <- tclvalue(permVariable)
        treatasdist <- tclvalue(treatasdistVariable)==1
        var2 <- svariables[as.numeric(tkcurselection(scaleBox))+1]
        if (test == "mantel") {
            doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
            if(treatasdist==F){
                logger(paste("distmatrix1 <- vegdist(", .communityDataSet, ",method='", dist1, "', na.rm=T)", sep=""))
                assign("distmatrix1", justDoIt(paste("vegdist(",.communityDataSet, ",method='",dist1, "', na.rm=T)", sep="")), envir=.GlobalEnv)
                doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist1, "')", sep=""))
            }else{
                logger(paste("distmatrix1 <- as.dist(", .communityDataSet, ")", sep=""))
                assign("distmatrix1", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
            }
            if (var2 == "all") {
                if(dist2 == "daisy (factor)") {
                    justDoIt(paste("library(cluster)"))
                    logger(paste("library(cluster)"))
                    logger(paste("distmatrix2 <- as.dist(as.matrix(daisy(", .activeDataSet, ")))", sep=""))
                    assign("distmatrix2", justDoIt(paste("as.dist(as.matrix(daisy(", .activeDataSet, ")))", sep="")), envir=.GlobalEnv)
                }else{
                    logger(paste("distmatrix2 <- vegdist(", .activeDataSet, ",method='",dist2, "', na.rm=T)", sep=""))
                    assign("distmatrix2", justDoIt(paste("vegdist(",.activeDataSet, ",method='",dist2, "', na.rm=T)", sep="")), envir=.GlobalEnv)                
                }
            }else{
                var2 <- .variables[as.numeric(tkcurselection(scaleBox))]
                varfactor <- eval(parse(text=paste("is.factor(",.activeDataSet, "$", var2, ")", sep="")), envir=.GlobalEnv)
                if (varfactor==T) {
                    justDoIt(paste("library(cluster)"))
                    logger(paste("library(cluster)"))
                    logger(paste("distmatrix2 <- as.dist(as.matrix(daisy(", .activeDataSet, "[,'", var2, "',drop=F])))", sep=""))
                    assign("distmatrix2", justDoIt(paste("distmatrix2 <- daisy(", .activeDataSet, "[,'", var2, "',drop=F])", sep="")), envir=.GlobalEnv)
                }else{
                    logger(paste("distmatrix2 <- vegdist(", .activeDataSet, "$", var2, ",method='",dist2, "')", sep=""))
                    assign("distmatrix2", justDoIt(paste("vegdist(", .activeDataSet, "$", var2,", method='",dist2, "')", sep="")), envir=.GlobalEnv)
                }
            }
            doItAndPrint(paste("mantel(distmatrix1, distmatrix2, method='", method, "', permutations=", permutations, ")",sep=""))
        }
        if (test == "anosim (factor)" && var2 != "all") {
            var2 <- .variables[as.numeric(tkcurselection(scaleBox))]
            varfactor <- eval(parse(text=paste("is.factor(",.activeDataSet, "$", var2, ")", sep="")), envir=.GlobalEnv)
            if (varfactor==T) {
                doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
                if(treatasdist==F){
                    logger(paste("distmatrix1 <- vegdist(", .communityDataSet, ", method='", dist1, "', na.rm=T)", sep=""))
                    assign("distmatrix1", justDoIt(paste("vegdist(", .communityDataSet, ", method='",dist1, "', na.rm=T)", sep="")), envir=.GlobalEnv)
                    doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist1, "')", sep=""))
                }else{
                    logger(paste("distmatrix1 <- as.dist(", .communityDataSet, ")", sep=""))
                    assign("distmatrix1", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
                }
                doItAndPrint(paste("summary(anosim(distmatrix1, grouping=", .activeDataSet, "$", var2, ", permutations=", permutations, "))",sep=""))
            }
        }
        if (test == "mrpp (factor)" && var2 != "all") {
            var2 <- .variables[as.numeric(tkcurselection(scaleBox))]
            varfactor <- eval(parse(text=paste("is.factor(",.activeDataSet, "$", var2, ")", sep="")), envir=.GlobalEnv)
            if (varfactor==T) {
                doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
                if(treatasdist==F){
                    logger(paste("distmatrix1 <- vegdist(", .communityDataSet, ",method='", dist1, "', na.rm=T)", sep=""))
                    assign("distmatrix1", justDoIt(paste("vegdist(",.communityDataSet, ", method='", dist1, "', na.rm=T)", sep="")), envir=.GlobalEnv)
                    doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist1, "')", sep=""))
                }else{
                    logger(paste("distmatrix1 <- as.dist(", .communityDataSet, ")", sep=""))
                    assign("distmatrix1", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
                }
                doItAndPrint(paste("mrpp(distmatrix1, grouping=", .activeDataSet, "$", var2, ", permutations=", permutations, ")",sep=""))
            }
        }
        if (test == "rankindex") {
            doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
            if (var2 == "all") {
                doItAndPrint(paste("rankindex(", .activeDataSet, ", ", .communityDataSet, ", method='",method, "')", sep=""))
            }else{
                var2 <- .variables[as.numeric(tkcurselection(scaleBox))]
                doItAndPrint(paste("rankindex(", .activeDataSet, "$", var2, ", ", .communityDataSet, ", method='", method, "')", sep=""))
            }
        }
        if (test == "bioenv") {
            doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
            doItAndPrint(paste("bioenv(", .communityDataSet, ", prepare.bioenv(", .activeDataSet, ", as.numeric=c()), method='", method, "', index='", dist1, "'", sep=""))
        }
        if (test == "betadisper (factor)" && var2 != "all") {
            var2 <- .variables[as.numeric(tkcurselection(scaleBox))]
            varfactor <- eval(parse(text=paste("is.factor(",.activeDataSet, "$", var2, ")", sep="")), envir=.GlobalEnv)
            if (varfactor==T) {
                doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
                if(treatasdist==F){
                    logger(paste("distmatrix1 <- vegdist(", .communityDataSet, ", method='", dist1, "', na.rm=T)", sep=""))
                    assign("distmatrix1", justDoIt(paste("vegdist(", .communityDataSet, ", method='", dist1, "', na.rm=T)", sep="")), envir=.GlobalEnv)
                    doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist1, "')", sep=""))
                }else{
                    logger(paste("distmatrix1 <- as.dist(", .communityDataSet, ")", sep=""))
                    assign("distmatrix1", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
                }
                doItAndPrint(paste("betadisper(distmatrix1, group=", .activeDataSet, "$", var2, ")",sep=""))
                doItAndPrint(paste("anova(betadisper(distmatrix1, group=", .activeDataSet, "$", var2, "))",sep=""))
            }
        }
        if (test == "meandist (factor)" && var2 != "all") {
            var2 <- .variables[as.numeric(tkcurselection(scaleBox))]
            varfactor <- eval(parse(text=paste("is.factor(",.activeDataSet, "$", var2, ")", sep="")), envir=.GlobalEnv)
            if (varfactor==T) {
                doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
                if(treatasdist==F){
                    logger(paste("distmatrix1 <- vegdist(", .communityDataSet, ", method='", dist1, "', na.rm=T)", sep=""))
                    assign("distmatrix1", justDoIt(paste("vegdist(", .communityDataSet, ", method='", dist1, "', na.rm=T)", sep="")), envir=.GlobalEnv)
                    doItAndPrint(paste("dist.eval(", .communityDataSet, ",'", dist1, "')", sep=""))
                }else{
                    logger(paste("distmatrix1 <- as.dist(", .communityDataSet, ")", sep=""))
                    assign("distmatrix1", justDoIt(paste("as.dist(",.communityDataSet, ")", sep="")), envir=.GlobalEnv)
                }
                doItAndPrint(paste("meandist(distmatrix1, grouping=", .activeDataSet, "$", var2, ")",sep=""))
                doItAndPrint(paste("summary(meandist(distmatrix1, grouping=", .activeDataSet, "$", var2, "))",sep=""))
            }
        }
        if (test == "simper (factor)" && var2 != "all") {
            var2 <- .variables[as.numeric(tkcurselection(scaleBox))]
            varfactor <- eval(parse(text=paste("is.factor(",.activeDataSet, "$", var2, ")", sep="")), envir=.GlobalEnv)
            if (varfactor==T) {
                doItAndPrint(paste("check.datasets(", .communityDataSet, ", ", .activeDataSet, ")", sep=""))
                doItAndPrint(paste("simper(", .communityDataSet, ", group=", .activeDataSet, "$", var2, ", permutations=", permutations, ")",sep=""))
                doItAndPrint(paste("summary(simper(", .communityDataSet, ", group=", .activeDataSet, "$", var2, ", permutations=", permutations, "))",sep=""))
            }
        }
        plotit <- tclvalue(plotVariable) == "1"
        var2 <- svariables[as.numeric(tkcurselection(scaleBox))+1]
        if (plotit==T && test=="mantel"  && var2=="all") {
            justDoIt(paste("plot(distmatrix2, distmatrix1,xlab='environmental distance',ylab='ecological distance')", sep=""))
            logger(paste("plot(distmatrix2, distmatrix1,xlab='environmental distance',ylab='ecological distance')", sep=""))
        }
        if (plotit==T && test=="mantel"  && var2!="all") {
            var2 <- .variables[as.numeric(tkcurselection(scaleBox))]
            varfactor <- eval(parse(text=paste("is.factor(",.activeDataSet, "$", var2, ")", sep="")), envir=.GlobalEnv)
            if (varfactor==T) {
                justDoIt(paste("plot(distmatrix1~as.factor(distmatrix2),xlab='environmental distance',ylab='ecological distance')", sep="")) 
                logger(paste("plot(distmatrix1~as.factor(distmatrix2),xlab='environmental distance',ylab='ecological distance')", sep="")) 
            }else{
                justDoIt(paste("plot(distmatrix2, distmatrix1,xlab='environmental distance',ylab='ecological distance')", sep=""))
                logger(paste("plot(distmatrix2, distmatrix1,xlab='environmental distance',ylab='ecological distance')", sep=""))
            }
        }
        if (plotit==T && test!="mantel" && test!="rankindex" && test!="bioenv" && var2!="all") {
            var2 <- .variables[as.numeric(tkcurselection(scaleBox))]
            varfactor <- eval(parse(text=paste("is.factor(",.activeDataSet, "$", var2, ")", sep="")), envir=.GlobalEnv)
            if (varfactor==T) {
                justDoIt(paste("library(cluster)"))
                logger(paste("library(cluster)"))
                logger(paste("distmatrix2 <- as.dist(as.matrix(daisy(", .activeDataSet, "[,'", var2, "',drop=F])))", sep=""))
                assign("distmatrix2", justDoIt(paste("distmatrix2 <- daisy(", .activeDataSet, "[,'", var2, "',drop=F])", sep="")), envir=.GlobalEnv)
                justDoIt(paste("plot(distmatrix1~as.factor(distmatrix2),xlab='environmental distance',ylab='ecological distance')", sep="")) 
                logger(paste("plot(distmatrix1~as.factor(distmatrix2),xlab='environmental distance',ylab='ecological distance')", sep="")) 
            }
        }
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
    }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(method1Frame, text="Type of test"), sticky="nw")
    tkgrid(testBox, testScroll,sticky="nw")
    tkgrid(tklabel(method3Frame, text="Community distance"), sticky="nw")
    tkgrid(dist1Box, dist1Scroll,sticky="nw")
    tkgrid(tklabel(method2Frame, text="Environmental distance"), sticky="w")
    tkgrid(dist2Box, dist2Scroll,sticky="w")
    tkgrid(tklabel(method4Frame, text="Environmental variable"), sticky="w")
    tkgrid(scaleBox, scaleScroll, sticky="w")
    tkgrid(treatasdistCheckBox, tklabel(method5Frame, text="as.dist(Community)", width=15), sticky="w")
    tkgrid(plotCheckBox, tklabel(method5Frame, text="plot results"), sticky="w")
    tkgrid(tklabel(method5Frame, text="permutations", width=10), perma, sticky="w")
    tkgrid(tklabel(method6Frame, text="correlation"), sticky="w")
    tkgrid(methodBox, methodScroll, sticky="nw")
    tkgrid(method1Frame, tklabel(methodFrame, text="", width=1), method2Frame, sticky="w")
    tkgrid(method3Frame, tklabel(methodFrame, text="", width=1), method4Frame, sticky="w")
    tkgrid(method5Frame, tklabel(methodFrame, text="", width=1), method6Frame, sticky="w")
    tkgrid(methodFrame, sticky="w")
    tkgrid(OKbutton, cancelButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(dist1Scroll, sticky="ns")
    tkgrid.configure(dist2Scroll, sticky="ns")
    tkgrid.configure(methodScroll, sticky="ns")
    tkgrid.configure(scaleScroll, sticky="ns")
    tkgrid.configure(testScroll, sticky="ns")
    tkselection.set(dist1Box, 0)
    tkselection.set(dist2Box, 0)
    tkselection.set(methodBox, 0)
    tkselection.set(scaleBox, 0)
    tkselection.set(testBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(methodBox)
    tkwait.window(top)
}


cepNamesCommunity <- function() {
    .communityDataSet <- CommunityDataSet()
    justDoIt(paste("colnames(", .communityDataSet, ") <- make.cepnames(colnames(", .communityDataSet, "))", sep=""))
    communityDataSet(.communityDataSet)
}

helpBiodiversityR <- function() {
    print(help("BiodiversityRGUI",help_type="html"))
}

allcitations <- function() {
    doItAndPrint(paste("loaded.citations()"))
}

browseTDAwebsite <- function() {
    browseURL("http://www.worldagroforestry.org/treesandmarkets/tree_diversity_analysis.asp")
}

browseTDAmanual <- function() {
    browseURL("http://www.worldagroforestry.org/downloads/publications/PDFs/kindt%20b2005.pdf")
}

###################
# GUI for ensemble suitability modelling
# March 2016


ens.start <- function() {
    browseURL(paste(file.path(path.package(package="BiodiversityR")[1], "doc"),
        "/", gettextRcmdr("Getting-Started-ensemble-GUI"), ".txt", sep=""))
}

dismo.pdf <- function() {
    browseURL(paste(system.file(package="dismo"), "/doc/sdm.pdf", sep=""))
}

ens.directory <- function() {
    logger(paste("Select directory with raster layers and (ideally) presence locations", sep=""))
    logger(paste("Results will be saved in subfolders of this directory", sep=""))
    logger(paste("As alternative, use: R Commander > File > Change working directory...", sep=""))
    doItAndPrint(paste("setwd(choose.dir())", sep=""))
    doItAndPrint(paste("getwd()", sep=""))
}

ens.workspace <- function(){
    logger(paste("Select .RData file", sep=""))
    doItAndPrint(paste("load(choose.files(default='*.RData', multi=F))", sep=""))
    putRcmdr("dialog.values", list())
    activateMenus()
}

ens.grd <- function(s.files) {
    if(is.null(s.files) == T) {stop("no files provided")}
    if(length(s.files) < 1) {stop("no files provided")}
    for (i in 1:length(s.files)) {
        working.raster <- raster::raster(s.files[i])
        print(working.raster)
        raster::writeRaster(working.raster, filename=working.raster@file@name, format="raster", overwrite=T)
    }
}

ens.grd.menu <- function() {
    logger(paste("Select files to be converted to grd format of raster package", sep=""))
    doItAndPrint(paste("selected.files <- choose.files()", sep=""))
    logger(paste("As alternative use: selected.files <- list.files(path=getwd(), pattern='.tif', full.names=TRUE)", sep=""))
    justDoIt(paste("selected.files <- normalizePath(selected.files, mustWork=F)", sep=""))
    doItAndPrint(paste("selected.files", sep=""))
    doItAndPrint(paste("ens.grd(selected.files)", sep=""))
}

if (exists("stack.list") == F) {assign("stack.list", NULL, envir=.GlobalEnv)}
stack.listP <- function() {return(!is.null(stack.list))}
if (exists("stack.focal") == F) {assign("stack.focal", NULL, envir=.GlobalEnv)}
stack.focalP <- function() {return(!is.null(stack.focal))}

update.stacklist <- function(){
    all.objects <- ls(envir=.GlobalEnv)
    stack.list2 <- stack.list
    for (i in 1:length(all.objects)) {
        eval1 <- eval(parse(text=paste(all.objects[i])), envir=.GlobalEnv)
        if(class(eval1) == "RasterStack") {stack.list2 <- c(stack.list2, all.objects[i])}
    }
    stack.list2 <- unique(stack.list2)
    stack.list2 <- stack.list2[stack.list2 != "eval1"]
    assign("eval1", NULL, envir=.GlobalEnv)
    assign("stack.list", stack.list2, envir=.GlobalEnv)
}

update.stacklist()

update.stacklist.menu <- function(){
    justDoIt(paste("update.stacklist()", sep=""))
    doItAndPrint(paste("levels(as.factor(stack.list))", sep=""))    
}

stack.create.GUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Create stack of environmental layers")
    modelName <- tclVar("calibration1")
    modelFrame <- tkframe(top, relief="groove", borderwidth=2)
    model <- tkentry(modelFrame, width=40, textvariable=modelName)
    all.varsFrame <- tkframe(top, relief="groove", borderwidth=2)
    all.vars1Frame <- tkframe(all.varsFrame)
    all.vars2Frame <- tkframe(all.varsFrame)
    all.varsVariable <- tclVar("0")
    all.varsCheckBox <- tkcheckbutton(all.vars2Frame, variable=all.varsVariable)
    dismoVariable <- tclVar("0")
    dismoCheckBox <- tkcheckbutton(all.vars2Frame, variable=dismoVariable)
    onOK <- function(){
        all.vars <- tclvalue(all.varsVariable) == "1"
        dismoValue <- tclvalue(dismoVariable) == "1"
        if (dismoValue == T) {all.vars <- T}
        if (all.vars==T) {
            if (dismoValue == T) {
                dismo.dir <- normalizePath(paste(system.file(package="dismo"), '/ex', sep=''), mustWork=F)
                assign("dismo.ex", dismo.dir, envir=.GlobalEnv)
                doItAndPrint(paste("selected.files <- list.files(path=dismo.ex, pattern='.grd', full.names=TRUE)", sep=""))
            }else{
                doItAndPrint(paste("selected.files <- list.files(path=getwd(), pattern='.grd', full.names=TRUE)", sep=""))
            }
        }else{
            doItAndPrint(paste("selected.files <- choose.files(default='*.grd')", sep=""))
        }
        if(length(selected.files) > 0) {
            justDoIt(paste("selected.files <- normalizePath(selected.files, mustWork=F)", sep="")) 
            doItAndPrint(paste("selected.files", sep=""))
            modelValue <- tclvalue(modelName)
            modelValue <- gsub(".", "_", modelValue, fixed=T)
            logger(paste(modelValue, " <- raster::stack(selected.files)", sep=""))
            stack.eval <- eval(parse(text=paste("raster::stack(selected.files)", sep="")), envir=.GlobalEnv)
            assign(modelValue, stack.eval, envir=.GlobalEnv)
            doItAndPrint(paste(modelValue, "@title <- '", modelValue, "'", sep="")) 
            stack.title <- as.character(eval(parse(text=paste(modelValue, "@title", sep="")), envir=.GlobalEnv))
            doItAndPrint(paste(modelValue, sep=""))
            stack.list <- c(stack.list, stack.title)
            stack.list <- unique(stack.list)
            assign("stack.list", stack.list, envir=.GlobalEnv)
            if (length(stack.list) == 1) {assign("stack.focal", stack.list[1], envir=.GlobalEnv)}
        }
        putRcmdr("dialog.values", list())
        activateMenus()
	closeDialog()
	tkfocus(CommanderWindow())
        }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('stack', package='raster', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(modelFrame, text="Save result as:    ", width=15), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(tklabel(all.vars1Frame, text="Use <OK> button to select (<CTRL>-select) layers"), sticky="w")
    tkgrid(all.varsCheckBox, tklabel(all.vars2Frame, text="select all 'grd' raster layers in the working directory"), sticky="w")
    tkgrid(dismoCheckBox, tklabel(all.vars2Frame, text="select all 'grd' raster layers from the dismo package examples"), sticky="w")
    tkgrid(all.vars1Frame, sticky="w")
    tkgrid(all.vars2Frame, sticky="w")
    tkgrid(all.varsFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
#    tkfocus(modelName)
    tkwait.window(top)
}

stack.select.GUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Select calibration stack")
    update.stacklist()
    variables <- stack.list
    varFrame <- tkframe(top, relief="groove", borderwidth=2)
    subsetBox <- tklistbox(varFrame, width=40, height=7,
        selectmode="single", background="white", exportselection="FALSE")
    subsetScroll <- tkscrollbar(varFrame, repeatinterval=5, command=function(...) tkyview(subsetBox, ...))
    tkconfigure(subsetBox, yscrollcommand=function(...) tkset(subsetScroll, ...))
    for (x in variables) tkinsert(subsetBox, "end", x)
    onOK <- function(){
        var <- variables[as.numeric(tkcurselection(subsetBox))+1]
        assign("stack.focal", var, envir=.GlobalEnv)
        doItAndPrint(paste(stack.focal, sep=""))
        doItAndPrint(paste(stack.focal, "@title", sep=""))
        putRcmdr("dialog.values", list())
        activateMenus()
	closeDialog()
	tkfocus(CommanderWindow())
        }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(varFrame, text="Select calibration stack"), sticky="w")
    tkgrid(subsetBox, subsetScroll, sticky="w")
    tkgrid(varFrame, sticky="w")
    tkgrid(OKbutton, tklabel(buttonsFrame, text="    "), cancelButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(subsetScroll, sticky="ns")
    if (is.null(stack.focal) == F) {
        tkselection.set(subsetBox, which(stack.focal == variables)-1)
    }else{
        tkselection.set(subsetBox, 0)
    }
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(subsetBox)
    tkwait.window(top)
}

viewstack <- function(){
    doItAndPrint(paste(stack.focal, sep=""))
    doItAndPrint(paste(stack.focal, "@title", sep=""))
}

stack.delete.GUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Delete (drop) layers from calibration stack")
    lvariables <- eval(parse(text=paste("names(", stack.focal, ")", sep="")), envir=.GlobalEnv)
    layerFrame <- tkframe(top, relief="groove", borderwidth=2)
    layerBox <- tklistbox(layerFrame, width=40, height=8,
        selectmode="multiple", background="white", exportselection="FALSE") 
    layerScroll <- tkscrollbar(layerFrame, repeatinterval=5, command=function(...) tkyview(layerBox, ...))
    tkconfigure(layerBox, yscrollcommand=function(...) tkset(layerScroll, ...))
    for (x in lvariables) tkinsert(layerBox, "end", x)
    onOK <- function(){
        layers <- lvariables[as.numeric(tkcurselection(layerBox))+1]
        for (i in 1:length(layers)) {
            logger(paste(stack.focal, " <- raster::dropLayer(", stack.focal, ", which(names(", stack.focal, ")=='", layers[i], "'))", sep=""))
            stack.eval <- eval(parse(text=paste("raster::dropLayer(", stack.focal, ", which(names(", stack.focal, ")=='", layers[i], "'))", sep="")), envir=.GlobalEnv)
            assign(stack.focal, stack.eval, envir=.GlobalEnv)
        }
        activateMenus()
	closeDialog()
	tkfocus(CommanderWindow())
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('dropLayer', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(layerFrame, text="Layers to delete (drop)"), sticky="w")
    tkgrid(layerBox, layerScroll, sticky="w")
    tkgrid(layerFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(layerScroll, sticky="ns")
#    tkselection.set(layerBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(layerBox)
    tkwait.window(top)
}

stack.rename.GUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Rename layers of calibration stack")
    lvariables <- eval(parse(text=paste("names(", stack.focal, ")", sep="")), envir=.GlobalEnv)
    layerFrame <- tkframe(top, relief="groove", borderwidth=2)
    layerBox <- tklistbox(layerFrame, width=50, height=15,
        selectmode="single", background="white", exportselection="FALSE") 
    layerScroll <- tkscrollbar(layerFrame, repeatinterval=5, command=function(...) tkyview(layerBox, ...))
    tkconfigure(layerBox, yscrollcommand=function(...) tkset(layerScroll, ...))
    for (x in lvariables) tkinsert(layerBox, "end", x)
    modelName <- tclVar("calibration1")
    modelFrame <- tkframe(top, relief="groove", borderwidth=2)
    model <- tkentry(modelFrame, width=40, textvariable=modelName)
    varName <- tclVar("bio0x")
    varFrame <- tkframe(top, relief="groove", borderwidth=2)
    var <- tkentry(varFrame, width=40, textvariable=varName)
    onOK <- function(){
        layer <- lvariables[as.numeric(tkcurselection(layerBox))+1]
        varValue <- tclvalue(varName)
        doItAndPrint(paste("names(", stack.focal, ")[which(names(", stack.focal, ")=='", layer, "')] <- '", varValue, "'", sep=""))
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('names', package='raster', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(layerFrame, text="Select layer to rename"), sticky="w")
    tkgrid(layerBox, layerScroll, sticky="w")
    tkgrid(layerFrame, sticky="w")
    tkgrid(tklabel(varFrame, text="New name:", width=10), var, sticky="w")
    tkgrid(varFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(layerScroll, sticky="ns")
#    tkselection.set(layerBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(layerBox)
    tkwait.window(top)
}

if (exists("stack.factors") == F) {assign("stack.factors", NULL, envir=.GlobalEnv)}

stack.factors.GUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Identify factors for calibration stack")
    lvariables <- eval(parse(text=paste("names(", stack.focal, ")", sep="")), envir=.GlobalEnv)
    layerFrame <- tkframe(top, relief="groove", borderwidth=2)
    layerBox <- tklistbox(layerFrame, width=40, height=8,
        selectmode="multiple", background="white", exportselection="FALSE") 
    layerScroll <- tkscrollbar(layerFrame, repeatinterval=5, command=function(...) tkyview(layerBox, ...))
    tkconfigure(layerBox, yscrollcommand=function(...) tkset(layerScroll, ...))
    for (x in lvariables) tkinsert(layerBox, "end", x)
    onOK <- function(){
        layers <- lvariables[as.numeric(tkcurselection(layerBox))+1]
        if (length(layers) < 1) {
            stack.factors0 <- NULL
        }else{
            stack.factors0 <- character(length=length(layers))
            for (i in 1:length(layers)) {stack.factors0[i] <- paste(layers[i], sep="")}
        }
        assign("stack.factors", stack.factors0, envir=.GlobalEnv)
        if (is.null(stack.factors) == F) {
            for (i in 1:length(stack.factors)) {
                if (i==1) {
                    factor.string <- paste("c('", stack.factors[i], "'", sep="")
                }else{
                    factor.string <- paste(factor.string, ", '", stack.factors[i], "'", sep="")
                }
            }
            factor.string <- paste(factor.string, ")", sep="")
            doItAndPrint(paste("Selected factors: ", factor.string, sep=""))
        }
        activateMenus()
	closeDialog()
	tkfocus(CommanderWindow())
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(layerFrame, text="Select factors"), sticky="w")
    tkgrid(layerBox, layerScroll, sticky="w")
    tkgrid(layerFrame, sticky="w")
    tkgrid(OKbutton, cancelButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(layerScroll, sticky="ns")
#    tkselection.set(layerBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(layerBox)
    tkwait.window(top)
}

if (exists("stack.dummies") == F) {assign("stack.dummies", NULL, envir=.GlobalEnv)}

stack.dummies.GUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Identify dummy variables for calibration stack")
    lvariables <- eval(parse(text=paste("names(", stack.focal, ")", sep="")), envir=.GlobalEnv)
    layerFrame <- tkframe(top, relief="groove", borderwidth=2)
    layerBox <- tklistbox(layerFrame, width=40, height=8,
        selectmode="multiple", background="white", exportselection="FALSE") 
    layerScroll <- tkscrollbar(layerFrame, repeatinterval=5, command=function(...) tkyview(layerBox, ...))
    tkconfigure(layerBox, yscrollcommand=function(...) tkset(layerScroll, ...))
    for (x in lvariables) tkinsert(layerBox, "end", x)
    onOK <- function(){
        layers <- lvariables[as.numeric(tkcurselection(layerBox))+1]
        if (length(layers) < 1) {
            stack.dummies0 <- NULL
        }else{
            stack.dummies0 <- character(length=length(layers))
            for (i in 1:length(layers)) {stack.dummies0[i] <- paste(layers[i], sep="")}
        }
        assign("stack.dummies", stack.dummies0, envir=.GlobalEnv)
        if (is.null(stack.dummies) == F) {
            for (i in 1:length(stack.dummies)) {
                if (i==1) {
                    dummy.string <- paste("c('", stack.dummies[i], "'", sep="")
                }else{
                    dummy.string <- paste(dummy.string, ", '", stack.dummies[i], "'", sep="")
                }
            }
            dummy.string <- paste(dummy.string, ")", sep="")
            doItAndPrint(paste("Selected dummy variables: ", dummy.string, sep=""))
        }
        activateMenus()
	closeDialog()
	tkfocus(CommanderWindow())
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(layerFrame, text="Select dummy variables"), sticky="w")
    tkgrid(layerBox, layerScroll, sticky="w")
    tkgrid(layerFrame, sticky="w")
    tkgrid(OKbutton, cancelButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(layerScroll, sticky="ns")
#    tkselection.set(layerBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(layerBox)
    tkwait.window(top)
}


presence_file <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
BradypusPresence <- read.table(presence_file, header=TRUE, sep=',')
BradypusPresence[, 'species'] <- as.factor(gsub(' ', '_', BradypusPresence[, 'species']))

if (exists("presence.focal") == F) {assign("presence.focal", NULL, envir=.GlobalEnv)}
presence.focalP <- function() {return(!is.null(presence.focal))}

make.presence.GUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Make presence data set")
    .activeDataSet <- ActiveDataSet()
    .fvariables <- Factors()
    fvariables <- paste(.fvariables)
    .nvariables <- Numeric()
    nvariables <- paste(.nvariables)
    modelName <- tclVar("Presence1")
    modelFrame <- tkframe(top, relief="groove", borderwidth=2)
    model <- tkentry(modelFrame, width=24, textvariable=modelName)
    speciesFrame <- tkframe(top, relief="groove", borderwidth=2)
    speciesBox <- tklistbox(speciesFrame, width=40, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    speciesScroll <- tkscrollbar(speciesFrame, repeatinterval=5, command=function(...) tkyview(speciesBox, ...))
    tkconfigure(speciesBox, yscrollcommand=function(...) tkset(speciesScroll, ...))
    for (x in fvariables) tkinsert(speciesBox, "end", x)
    xlonFrame <- tkframe(top, relief="groove", borderwidth=2)
    xlonBox <- tklistbox(xlonFrame, width=40, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    xlonScroll <- tkscrollbar(xlonFrame, repeatinterval=5, command=function(...) tkyview(xlonBox, ...))
    tkconfigure(xlonBox, yscrollcommand=function(...) tkset(xlonScroll, ...))
    for (x in nvariables) tkinsert(xlonBox, "end", x)
    ylatFrame <- tkframe(top, relief="groove", borderwidth=2)
    ylatBox <- tklistbox(ylatFrame, width=40, height=5,
        selectmode="single", background="white", exportselection="FALSE") 
    ylatScroll <- tkscrollbar(ylatFrame, repeatinterval=5, command=function(...) tkyview(ylatBox, ...))
    tkconfigure(ylatBox, yscrollcommand=function(...) tkset(ylatScroll, ...))
    for (x in nvariables) tkinsert(ylatBox, "end", x)
    onOK <- function(){
        modelValue <- tclvalue(modelName)
        species <- .fvariables[as.numeric(tkcurselection(speciesBox))+1]
        xlon <- .nvariables[as.numeric(tkcurselection(xlonBox))+1]
        ylat <- .nvariables[as.numeric(tkcurselection(ylatBox))+1]
        command <- paste(.activeDataSet, "[as.numeric(na.omit(match(c('", species, "', '", xlon, "', '", ylat, "'), names(", .activeDataSet, "))))]", sep="")
        logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        doItAndPrint(paste("names(", modelValue, ") <- c('species', 'x', 'y')", sep=""))
        assign("presence.focal", modelValue, envir=.GlobalEnv)
        doItAndPrint(paste(presence.focal, "[, 'species'] <- as.factor(gsub(' ', '_', ", presence.focal, "[, 'species']))", sep=""))
        doItAndPrint(paste("summary(", presence.focal, ")", sep=""))
        putRcmdr("dialog.values", list())
        activateMenus()
	closeDialog()
	tkfocus(CommanderWindow())
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(modelFrame, text="Save result as:    ", width=15), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(tklabel(speciesFrame, text="Species variable"), sticky="w")
    tkgrid(speciesBox, speciesScroll, sticky="w")
    tkgrid(speciesFrame, sticky="w")
    tkgrid(tklabel(xlonFrame, text="x variable (eg, longitude)"), sticky="w")
    tkgrid(xlonBox, xlonScroll, sticky="w")
    tkgrid(xlonFrame, sticky="w")
    tkgrid(tklabel(ylatFrame, text="y variable (eg, latitude)"), sticky="w")
    tkgrid(ylatBox, ylatScroll, sticky="w")
    tkgrid(ylatFrame, sticky="w")
    tkgrid(OKbutton, cancelButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(speciesScroll, sticky="ns")
    tkgrid.configure(xlonScroll, sticky="ns")
    tkgrid.configure(ylatScroll, sticky="ns")
    tkselection.set(speciesBox, 0)
    tkselection.set(xlonBox, 0)
    tkselection.set(ylatBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(speciesBox)
    tkwait.window(top)
}

presence.select.GUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Select presence data set")
    variables <- listDataSets()
    varFrame <- tkframe(top, relief="groove", borderwidth=2)
    subsetBox <- tklistbox(varFrame, width=40, height=7,
        selectmode="single", background="white", exportselection="FALSE")
    subsetScroll <- tkscrollbar(varFrame, repeatinterval=5, command=function(...) tkyview(subsetBox, ...))
    tkconfigure(subsetBox, yscrollcommand=function(...) tkset(subsetScroll, ...))
    for (x in variables) tkinsert(subsetBox, "end", x)
    onOK <- function(){
        var <- variables[as.numeric(tkcurselection(subsetBox))+1]
        assign("presence.focal", var, envir=.GlobalEnv)
        presence.names <- eval(parse(text=paste("names(", presence.focal, ")", sep="")), envir=.GlobalEnv)
        presence.names <- presence.names[1:3]
        if (all.equal(presence.names, c('species', 'lon', 'lat')) == T) {presence.names <- c('species', 'x', 'y')}        
        if (all.equal(presence.names, c('species', 'x', 'y')) == F) {doItAndPrint(paste("WARNING: variable names are not 'species', 'x' and 'y'", sep=""))}
        doItAndPrint(paste(presence.focal, "[, 'species'] <- as.factor(gsub(' ', '_', ", presence.focal, "[, 'species']))", sep=""))
        doItAndPrint(paste("summary(", presence.focal, ")", sep=""))
        putRcmdr("dialog.values", list())
        activateMenus()
	closeDialog()
	tkfocus(CommanderWindow())
        }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(varFrame, text="Select presence data set"), sticky="w")
    tkgrid(subsetBox, subsetScroll, sticky="w")
    tkgrid(varFrame, sticky="w")
    tkgrid(OKbutton, tklabel(buttonsFrame, text="    "), cancelButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(subsetScroll, sticky="ns")
    if (is.null(presence.focal) == F) {
        tkselection.set(subsetBox, which(presence.focal == variables)-1)
    }else{
        tkselection.set(subsetBox, 0)
    }
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(subsetBox)
    tkwait.window(top)
}

viewpresence <- function(){
    command <- justDoIt(paste("invisible(edit(", presence.focal, "))", sep=""))
}

make.absences <- function(an=1000, x, excludep=F, presence.data=NULL) {
    if (excludep == F) {
        as <- data.frame(dismo::randomPoints(x[[1]], n=an, p=NULL, excludep=F))
        names(as) <- c("x", "y")
        return(as)
    }else{
        species.names <- levels(droplevels(factor(presence.data[, 1])))
        for (s in 1:length(species.names)) {
            focal.species <- species.names[s]
            ps <- presence.data[presence.data[, 1]==focal.species, c(2:3)]
            as1 <- data.frame(dismo::randomPoints(x[[1]], n=an, p=ps, excludep=T))
            as2 <- cbind(rep(focal.species, length=nrow(as1)), as1)
            names(as2) <- c("species", "x", "y")
            if(s == 1) {
                as3 <- as2
            }else{
                as3 <- rbind(as3, as2)
            }
        } 
        return(as3)
    }
}

predictor.files <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''),
    pattern='grd', full.names=TRUE)
predictors <- raster::stack(predictor.files)
BradypusAbsence <- make.absences(x=predictors, excludep=T, presence.data=BradypusPresence)
assign("predictors", NULL, envir=.GlobalEnv)

if (exists("absence.focal") == F) {assign("absence.focal", NULL, envir=.GlobalEnv)}
absence.focalP <- function() {return(!is.null(absence.focal))}

make.absence.GUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Create absence data set")
    modelName <- tclVar("Absence1")
    modelFrame <- tkframe(top, relief="groove", borderwidth=2)
    model <- tkentry(modelFrame, width=24, textvariable=modelName)
    excludeFrame <- tkframe(top, relief="groove", borderwidth=2)
    exclude1Frame <- tkframe(excludeFrame)
    exclude2Frame <- tkframe(excludeFrame)
    n.abs <- tclVar("1000")
    n.absEntry <- tkentry(exclude1Frame, width=8, textvariable=n.abs)
    excludeVariable <- tclVar("0")
    excludeCheckBox <- tkcheckbutton(exclude2Frame, variable=excludeVariable)
    onOK <- function(){
        modelValue <- tclvalue(modelName)
        n.absValue <- tclvalue(n.abs)
        exclude <- tclvalue(excludeVariable) == "1"
        if (exclude==T) {
            command <- paste("make.absences(an=", n.absValue, ", x=", stack.focal, ", excludep=T, presence.data=", presence.focal, ")", sep="")
        }else{
            command <- paste("make.absences(an=", n.absValue, ", x=", stack.focal, ", excludep=F)", sep="")
        }
        logger(paste(modelValue, " <- ", command, sep=""))
        assign(modelValue, justDoIt(command), envir=.GlobalEnv)
        if (exclude==T) {
            doItAndPrint(paste("names(", modelValue, ") <- c('species', 'x', 'y')", sep=""))
        }else{
            doItAndPrint(paste("names(", modelValue, ") <- c('x', 'y')", sep=""))
        }
        if (is.null(absence.focal) == T) {assign("absence.focal", modelValue, envir=.GlobalEnv)}
        doItAndPrint(paste("summary(", absence.focal, ")", sep="")) 
        putRcmdr("dialog.values", list())
        activateMenus()
	closeDialog()
	tkfocus(CommanderWindow())
        }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('randomPoints', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(modelFrame, text="Save result as:    ", width=15), model, sticky="w")
    tkgrid(modelFrame, sticky="w")
    tkgrid(tklabel(exclude1Frame, text="number of absence points: ", width=26), n.absEntry, sticky="w")
    tkgrid(excludeCheckBox, tklabel(exclude2Frame, text="exclude raster cells with presence        "), sticky="w")
    tkgrid(exclude1Frame, sticky="w")
    tkgrid(exclude2Frame, sticky="w")
    tkgrid(excludeFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
#    tkfocus(modelName)
    tkwait.window(top)
}

absence.select.GUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Select absence data set")
    variables <- listDataSets()
    varFrame <- tkframe(top, relief="groove", borderwidth=2)
    subsetBox <- tklistbox(varFrame, width=40, height=7,
        selectmode="single", background="white", exportselection="FALSE")
    subsetScroll <- tkscrollbar(varFrame, repeatinterval=5, command=function(...) tkyview(subsetBox, ...))
    tkconfigure(subsetBox, yscrollcommand=function(...) tkset(subsetScroll, ...))
    for (x in variables) tkinsert(subsetBox, "end", x)
    onOK <- function(){
        var <- variables[as.numeric(tkcurselection(subsetBox))+1]
        assign("absence.focal", var, envir=.GlobalEnv)
        absence.names <- eval(parse(text=paste("names(", absence.focal, ")", sep="")), envir=.GlobalEnv)
        if (length(absence.names) == 3) {
            if (all.equal(absence.names, c('species', 'lon', 'lat')) == T) {absence.names <- c('species', 'x', 'y')}        
            if (all.equal(absence.names, c('species', 'x', 'y')) == F) {doItAndPrint(paste("WARNING: variable names are not 'species', 'x' and 'y'", sep=""))}
            doItAndPrint(paste(absence.focal, "[, 'species'] <- as.factor(gsub(' ', '_', ", absence.focal, "[, 'species']))", sep=""))
        }
        doItAndPrint(paste("summary(", absence.focal, ")", sep="")) 
        putRcmdr("dialog.values", list())
        activateMenus()
	closeDialog()
	tkfocus(CommanderWindow())
        }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(varFrame, text="Select absence data set"), sticky="w")
    tkgrid(subsetBox, subsetScroll, sticky="w")
    tkgrid(varFrame, sticky="w")
    tkgrid(OKbutton, tklabel(buttonsFrame, text="    "), cancelButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(subsetScroll, sticky="ns")
    if (is.null(absence.focal) == F) {
        tkselection.set(subsetBox, which(absence.focal == variables)-1)
    }else{
        tkselection.set(subsetBox, 0)
    }
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(subsetBox)
    tkwait.window(top)
}

viewabsence <- function(){
    command <- justDoIt(paste("invisible(edit(", absence.focal, "))", sep=""))
}

if (exists("ensmodels.file") == F) {assign("ensmodels.file", NULL, envir=.GlobalEnv)}
if (exists("ensmodels.focal") == F) {assign("ensmodels.focal", NULL, envir=.GlobalEnv)}
ensmodelsP <- function() {return(!is.null(ensmodels.focal))}

batch.GUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Ensemble suitability modelling")
    firstFrame <- tkframe(top, relief="groove", borderwidth=2)
    leftFrame <- tkframe(firstFrame)
    left1Frame <- tkframe(leftFrame)
    left2Frame <- tkframe(leftFrame) 
    n.ensemblesVariable <- tclVar("1")
    n.ensembleEntry <- tkentry(left1Frame, width=15, textvariable=n.ensemblesVariable)
    k.splitsVariable <- tclVar("4")
    k.splitsEntry <- tkentry(left1Frame, width=15, textvariable=k.splitsVariable)
    ENSEMBLE.bestVariable <- tclVar("0")
    ENSEMBLE.bestEntry <- tkentry(left1Frame, width=15, textvariable=ENSEMBLE.bestVariable)
    ENSEMBLE.minVariable <- tclVar("0.7")
    ENSEMBLE.minEntry <- tkentry(left1Frame, width=15, textvariable=ENSEMBLE.minVariable)
    ENSEMBLE.exponentVariable <- tclVar("c(1, 2)")
    ENSEMBLE.exponentEntry <- tkentry(left1Frame, width=15, textvariable=ENSEMBLE.exponentVariable)
    PROBITVariable <- tclVar("0")
    PROBITCheckBox <- tkcheckbutton(left2Frame, variable=PROBITVariable)
    AUC.weightsVariable <- tclVar("1")
    AUC.weightsCheckBox <- tkcheckbutton(left2Frame, variable=AUC.weightsVariable)
    infoFrame <- tkframe(firstFrame)
    info1Frame <- tkframe(infoFrame)
    if (is.null(stack.focal) == F) {stack.focalVariable <- tclVar(stack.focal)}
    if (is.null(stack.focal) == T) {stack.focalVariable <- tclVar("(insert calibration stack)")}
    stack.focalEntry <- tkentry(info1Frame, width=25, textvariable=stack.focalVariable)
    if (is.null(presence.focal) == F) {presence.focalVariable <- tclVar(presence.focal)}
    if (is.null(presence.focal) == T) {presence.focalVariable <- tclVar("(insert presence)")}
    presence.focalEntry <- tkentry(info1Frame, width=25, textvariable=presence.focalVariable)
    if (is.null(absence.focal) == F) {absence.focalVariable <- tclVar(absence.focal)}
    if (is.null(absence.focal) == T) {absence.focalVariable <- tclVar("(insert absence)")}
    absence.focalEntry <- tkentry(info1Frame, width=25, textvariable=absence.focalVariable)
    info2Frame <- tkframe(infoFrame)
    secondFrame <- tkframe(top, relief="groove", borderwidth=2)
    stackFrame <- tkframe(secondFrame, relief="groove", borderwidth=2)
    update.stacklist()
    if (is.null(stack.list) == F) {s.variables <- stack.list}
    if (is.null(stack.list) == T) {s.variables <- c("(none available)")}
    stackBox <- tklistbox(stackFrame, width=40, height=8,
        selectmode="multiple", background="white", exportselection="FALSE") 
    stackScroll <- tkscrollbar(stackFrame, repeatinterval=5, command=function(...) tkyview(stackBox, ...))
    tkconfigure(stackBox, yscrollcommand=function(...) tkset(stackScroll, ...))
    for (x in s.variables) tkinsert(stackBox, "end", x)
    rightFrame <- tkframe(secondFrame)
    right1Frame <- tkframe(rightFrame)
    right2Frame <- tkframe(rightFrame)
    right3Frame <- tkframe(rightFrame)
    variables <- c("spec_sens", "equal_sens_spec", "sensitivity", "no_omission", "prevalence", "kappa", "threshold.min", "threshold.mean", 
        "Sens=Spec", "MaxSens+Spec", "MaxKappa", "MaxPCC", "PredPrev=Obs", "ObsPrev", "MeanProb", "MinROCdist", "ReqSens")
    varFrame <- tkframe(right1Frame, relief="groove", borderwidth=2)
    subsetBox <- tklistbox(varFrame, width=40, height=5,
        selectmode="single", background="white", exportselection="FALSE")
    subsetScroll <- tkscrollbar(varFrame, repeatinterval=5, command=function(...) tkyview(subsetBox, ...))
    tkconfigure(subsetBox, yscrollcommand=function(...) tkset(subsetScroll, ...))
    for (x in variables) tkinsert(subsetBox, "end", x)
    sensitivityVariable <- tclVar("0.9")
    sensitivityEntry <- tkentry(right2Frame, width=15, textvariable=sensitivityVariable)
    PRESABSVariable <- tclVar("0")
    PRESABSCheckBox <- tkcheckbutton(right3Frame, variable=PRESABSVariable)
    algoFrame <- tkframe(top, relief="groove", borderwidth=2)
    algo1Frame <- tkframe(algoFrame)
    algo2Frame <- tkframe(algoFrame)
    algo3Frame <- tkframe(algoFrame)
    algo4Frame <- tkframe(algoFrame)
    MAXENTVariable <- tclVar("0")
    MAXENTCheckBox <- tkcheckbutton(algo1Frame, variable=MAXENTVariable)
    GLMVariable <- tclVar("0")
    GLMCheckBox <- tkcheckbutton(algo1Frame, variable=GLMVariable)
    MGCVVariable <- tclVar("0")
    MGCVCheckBox <- tkcheckbutton(algo1Frame, variable=MGCVVariable)
    NNETVariable <- tclVar("0")
    NNETCheckBox <- tkcheckbutton(algo1Frame, variable=NNETVariable)
    BIOCLIMVariable <- tclVar("0")
    BIOCLIMCheckBox <- tkcheckbutton(algo1Frame, variable=BIOCLIMVariable)
    GBMVariable <- tclVar("0")
    GBMCheckBox <- tkcheckbutton(algo2Frame, variable=GBMVariable)
    GLMSTEPVariable <- tclVar("0")
    GLMSTEPCheckBox <- tkcheckbutton(algo2Frame, variable=GLMSTEPVariable)
    MGCVFIXVariable <- tclVar("0")
    MGCVFIXCheckBox <- tkcheckbutton(algo2Frame, variable=MGCVFIXVariable)
    FDAVariable <- tclVar("0")
    FDACheckBox <- tkcheckbutton(algo2Frame, variable=FDAVariable)
    DOMAINVariable <- tclVar("0")
    DOMAINCheckBox <- tkcheckbutton(algo2Frame, variable=DOMAINVariable)
    GBMSTEPVariable <- tclVar("0")
    GBMSTEPCheckBox <- tkcheckbutton(algo3Frame, variable=GBMSTEPVariable)
    GAMVariable <- tclVar("0")
    GAMCheckBox <- tkcheckbutton(algo3Frame, variable=GAMVariable)
    EARTHVariable <- tclVar("0")
    EARTHCheckBox <- tkcheckbutton(algo3Frame, variable=EARTHVariable)
    SVMVariable <- tclVar("0")
    SVMCheckBox <- tkcheckbutton(algo3Frame, variable=SVMVariable)
    MAHALVariable <- tclVar("0")
    MAHALCheckBox <- tkcheckbutton(algo3Frame, variable=MAHALVariable)
    RFVariable <- tclVar("0")
    RFCheckBox <- tkcheckbutton(algo4Frame, variable=RFVariable)
    GAMSTEPVariable <- tclVar("0")
    GAMSTEPCheckBox <- tkcheckbutton(algo4Frame, variable=GAMSTEPVariable)
    RPARTVariable <- tclVar("0")
    RPARTCheckBox <- tkcheckbutton(algo4Frame, variable=RPARTVariable)
    SVMEVariable <- tclVar("0")
    SVMECheckBox <- tkcheckbutton(algo4Frame, variable=SVMEVariable)
    if (is.null(presence.focal) == F) {species.names <- eval(parse(text=paste("levels(droplevels(factor(", presence.focal, "[, 1])))", sep="")), envir=.GlobalEnv)}
    if (is.null(presence.focal) == F) {species.last <- species.names[length(species.names)]}
    onOK <- function(){
        n.ensembles <- tclvalue(n.ensemblesVariable)
        k.splits <- tclvalue(k.splitsVariable)
        ENSEMBLE.best <- tclvalue(ENSEMBLE.bestVariable)
        ENSEMBLE.min <- tclvalue(ENSEMBLE.minVariable)
        ENSEMBLE.exponent <- tclvalue(ENSEMBLE.exponentVariable)
        PROBIT <- tclvalue(PROBITVariable) == "1"
        AUC.weights <- tclvalue(AUC.weightsVariable) == "1"
        stack.focalValue <- tclvalue(stack.focalVariable)
        presence.focalValue <- tclvalue(presence.focalVariable)
        absence.focalValue <- tclvalue(absence.focalVariable)
        var <- variables[as.numeric(tkcurselection(subsetBox))+1]
        sensitivity <- tclvalue(sensitivityVariable)
        PRESABS <- tclvalue(PRESABSVariable) == "1"
        if (var == "Sens=Spec"){PRESABS <- TRUE}
        if (var == "MaxSens+Spec"){PRESABS <- TRUE}
        if (var == "MaxKappa"){PRESABS <- TRUE}
        if (var == "MaxPCC"){PRESABS <- TRUE}
        if (var == "PredPrev=Obs"){PRESABS <- TRUE}
        if (var == "ObsPrev"){PRESABS <- TRUE}
        if (var == "MeanProb"){PRESABS <- TRUE}
        if (var == "MinROCdist"){PRESABS <- TRUE}
        if (var == "ReqSens"){PRESABS <- TRUE}
        MAXENT <- tclvalue(MAXENTVariable)
        GBM <- tclvalue(GBMVariable)
        GBMSTEP <- tclvalue(GBMSTEPVariable)
        RF <- tclvalue(RFVariable)
        GLM <- tclvalue(GLMVariable)
        GLMSTEP <- tclvalue(GLMSTEPVariable)
        GAM <- tclvalue(GAMVariable)
        GAMSTEP <- tclvalue(GAMSTEPVariable)
        MGCV <- tclvalue(MGCVVariable)
        MGCVFIX <- tclvalue(MGCVFIXVariable)
        EARTH <- tclvalue(EARTHVariable)
        RPART <- tclvalue(RPARTVariable)
        NNET <- tclvalue(NNETVariable)
        FDA <- tclvalue(FDAVariable)
        SVM <- tclvalue(SVMVariable)
        SVME <- tclvalue(SVMEVariable)
        BIOCLIM <- tclvalue(BIOCLIMVariable)
        DOMAIN <- tclvalue(DOMAINVariable)
        MAHAL <- tclvalue(MAHALVariable)
        if (is.null(stack.factors) == F) {
            for (i in 1:length(stack.factors)) {
                if (i==1) {
                    factor.string <- paste("c('", stack.factors[i], "'", sep="")
                }else{
                    factor.string <- paste(factor.string, ", '", stack.factors[i], "'", sep="")
                }
            }
            factor.string <- paste(factor.string, ")", sep="")
        }else{
            factor.string <- paste("c()", sep="")
        }
        if (is.null(stack.dummies) == F) {
            for (i in 1:length(stack.dummies)) {
                if (i==1) {
                    dummy.string <- paste("c('", stack.dummies[i], "'", sep="")
                }else{
                    dummy.string <- paste(dummy.string, ", '", stack.dummies[i], "'", sep="")
                }
            }
            dummy.string <- paste(dummy.string, ")", sep="")
        }else{
            dummy.string <- paste("c()", sep="")
        }
        stacks <- s.variables[as.numeric(tkcurselection(stackBox))+1]
        if (is.null(stack.list) == T) {stacks <- tclvalue(stack.focalVariable)}
        if (length(stacks) > 0) {
            for (i in 1:length(stacks)) {
                if (i==1) {
                    stack.string <- paste("c(", stacks[i], sep="")
                }else{
                    stack.string <- paste(stack.string, ", ", stacks[i], sep="")
                }
            }
            stack.string <- paste(stack.string, ")", sep="")
        }else{
            stack.string <- paste("c()", sep="")
        }
        logger(paste("Note that it can take a while before results will be shown", sep=""))
        logger(paste("When calculations and projections are finalized, the window interface will close", sep=""))
        logger(paste("You can also monitor progress in the 'ensembles' and 'models' subfolders of the working directory", sep=""))
        doItAndPrint(paste("ensemble.batch(x=", stack.focalValue, ", xn=", stack.string, 
            ", species.presence=", presence.focalValue, ", species.absence=", absence.focalValue,
            ", presence.min=20, k.splits=", k.splits,
            ", n.ensembles=", n.ensembles, ", KML.out=TRUE, models.save=TRUE",
            ", threshold.method='", var, "', threshold.sensitivity=", sensitivity, 
            ", threshold.PresenceAbsence=", PRESABS,
            ", factors=", factor.string, ", dummy.vars=", dummy.string,
            ", ENSEMBLE.best=", ENSEMBLE.best, ", ENSEMBLE.min=", ENSEMBLE.min, ", ENSEMBLE.exponent=", ENSEMBLE.exponent,
            ", MAXENT=", MAXENT, ", GBM=", GBM, ", GBMSTEP=", GBMSTEP, ", RF=", RF, ", GLM=", GLM, ", GLMSTEP=", GLMSTEP, ", GAM=", GAM,
            ", GAMSTEP=", GAMSTEP, ", MGCV=", MGCV, ", MGCVFIX=", MGCVFIX, ", EARTH=", EARTH, ", RPART=", RPART, ", NNET=", NNET,
            ", FDA=", FDA, ", SVM=", SVM, ", SVME=", SVME, ", BIOCLIM=", BIOCLIM, ", DOMAIN=", DOMAIN, ", MAHAL=", MAHAL,
            ", PROBIT=", PROBIT, ", AUC.weights=", AUC.weights, ")", sep=""))
        if (is.null(presence.focal) == F) {
        if (n.ensembles > 1) {
            models.file <- paste(getwd(), "//models//", species.last, "_ENSEMBLE_", n.ensembles, "_models", sep="")
            models.file <- normalizePath(models.file, mustWork=F)
            if (file.exists(models.file) == T) {
                assign("ensmodels.file", models.file, envir=.GlobalEnv)
                load(ensmodels.file)
                assign("ensmodels.focal", ensemble.models, envir=.GlobalEnv)          
                logger(paste("Focal ensemble models (object ensmodels.focal) loaded from: ", models.file, sep=""))
            }
        }else{
            models.file <- paste(getwd(), "//models//", species.last, "_models", sep="")
            models.file <- normalizePath(models.file, mustWork=F)
            if (file.exists(models.file) == T) {
                assign("ensmodels.file", models.file, envir=.GlobalEnv)
                load(ensmodels.file)
                assign("ensmodels.focal", ensemble.models, envir=.GlobalEnv)          
                logger(paste("Focal ensemble models (object ensmodels.focal) loaded from: ", models.file, sep=""))
            }
        }
        }
        putRcmdr("dialog.values", list())
        activateMenus()
	closeDialog()
	tkfocus(CommanderWindow())
        }
    onCapture <- function(){
        n.ensembles <- tclvalue(n.ensemblesVariable)
        k.splits <- tclvalue(k.splitsVariable)
        ENSEMBLE.best <- tclvalue(ENSEMBLE.bestVariable)
        ENSEMBLE.min <- tclvalue(ENSEMBLE.minVariable)
        ENSEMBLE.exponent <- tclvalue(ENSEMBLE.exponentVariable)
        PROBIT <- tclvalue(PROBITVariable) == "1"
        AUC.weights <- tclvalue(AUC.weightsVariable) == "1"
        stack.focalValue <- tclvalue(stack.focalVariable)
        presence.focalValue <- tclvalue(presence.focalVariable)
        absence.focalValue <- tclvalue(absence.focalVariable)
        var <- variables[as.numeric(tkcurselection(subsetBox))+1]
        sensitivity <- tclvalue(sensitivityVariable)
        PRESABS <- tclvalue(PRESABSVariable) == "1"
        if (var == "Sens=Spec"){PRESABS <- TRUE}
        if (var == "MaxSens+Spec"){PRESABS <- TRUE}
        if (var == "MaxKappa"){PRESABS <- TRUE}
        if (var == "MaxPCC"){PRESABS <- TRUE}
        if (var == "PredPrev=Obs"){PRESABS <- TRUE}
        if (var == "ObsPrev"){PRESABS <- TRUE}
        if (var == "MeanProb"){PRESABS <- TRUE}
        if (var == "MinROCdist"){PRESABS <- TRUE}
        if (var == "ReqSens"){PRESABS <- TRUE}
        MAXENT <- tclvalue(MAXENTVariable)
        GBM <- tclvalue(GBMVariable)
        GBMSTEP <- tclvalue(GBMSTEPVariable)
        RF <- tclvalue(RFVariable)
        GLM <- tclvalue(GLMVariable)
        GLMSTEP <- tclvalue(GLMSTEPVariable)
        GAM <- tclvalue(GAMVariable)
        GAMSTEP <- tclvalue(GAMSTEPVariable)
        MGCV <- tclvalue(MGCVVariable)
        MGCVFIX <- tclvalue(MGCVFIXVariable)
        EARTH <- tclvalue(EARTHVariable)
        RPART <- tclvalue(RPARTVariable)
        NNET <- tclvalue(NNETVariable)
        FDA <- tclvalue(FDAVariable)
        SVM <- tclvalue(SVMVariable)
        SVME <- tclvalue(SVMEVariable)
        BIOCLIM <- tclvalue(BIOCLIMVariable)
        DOMAIN <- tclvalue(DOMAINVariable)
        MAHAL <- tclvalue(MAHALVariable)
        if (is.null(stack.factors) == F) {
            for (i in 1:length(stack.factors)) {
                if (i==1) {
                    factor.string <- paste("c('", stack.factors[i], "'", sep="")
                }else{
                    factor.string <- paste(factor.string, ", '", stack.factors[i], "'", sep="")
                }
            }
            factor.string <- paste(factor.string, ")", sep="")
        }else{
            factor.string <- paste("c()", sep="")
        }
        if (is.null(stack.dummies) == F) {
            for (i in 1:length(stack.dummies)) {
                if (i==1) {
                    dummy.string <- paste("c('", stack.dummies[i], "'", sep="")
                }else{
                    dummy.string <- paste(dummy.string, ", '", stack.dummies[i], "'", sep="")
                }
            }
            dummy.string <- paste(dummy.string, ")", sep="")
        }else{
            dummy.string <- paste("c()", sep="")
        }
        stacks <- s.variables[as.numeric(tkcurselection(stackBox))+1]
        if (is.null(stack.list) == T) {stacks <- tclvalue(stack.focalVariable)}
        if (length(stacks) > 0) {
            for (i in 1:length(stacks)) {
                if (i==1) {
                    stack.string <- paste("c(", stacks[i], sep="")
                }else{
                    stack.string <- paste(stack.string, ", ", stacks[i], sep="")
                }
            }
            stack.string <- paste(stack.string, ")", sep="")
        }else{
            stack.string <- paste("c()", sep="")
        }
        dir.create("outputs", showWarnings = F)
        filename1 <- paste(path=getwd(), "//outputs//", presence.focal, "_output.txt", sep="")
        filename1 <- normalizePath(filename1, mustWork=F)
        assign("output.filename", filename1, envir=.GlobalEnv)
        logger(paste("Results will be written to file: ", output.filename, sep=""))
        logger(paste("When calculations and projections are finalized, the window interface will close", sep=""))
        logger(paste("You can also monitor progress in the 'ensembles' and 'models' subfolders of the working directory", sep=""))
        doItAndPrint(paste("capture.output(ensemble.batch(x=", stack.focalValue, ", xn=", stack.string, 
            ", species.presence=", presence.focalValue, ", species.absence=", absence.focalValue,
            ", presence.min=20, k.splits=", k.splits,
            ", n.ensembles=", n.ensembles, ", KML.out=TRUE, models.save=TRUE",
            ", threshold.method='", var, "', threshold.sensitivity=", sensitivity, 
            ", threshold.PresenceAbsence=", PRESABS,
            ", factors=", factor.string, ", dummy.vars=", dummy.string,
            ", ENSEMBLE.best=", ENSEMBLE.best, ", ENSEMBLE.min=", ENSEMBLE.min, ", ENSEMBLE.exponent=", ENSEMBLE.exponent,
            ", MAXENT=", MAXENT, ", GBM=", GBM, ", GBMSTEP=", GBMSTEP, ", RF=", RF, ", GLM=", GLM, ", GLMSTEP=", GLMSTEP, ", GAM=", GAM,
            ", GAMSTEP=", GAMSTEP, ", MGCV=", MGCV, ", MGCVFIX=", MGCVFIX, ", EARTH=", EARTH, ", RPART=", RPART, ", NNET=", NNET,
            ", FDA=", FDA, ", SVM=", SVM, ", SVME=", SVME, ", BIOCLIM=", BIOCLIM, ", DOMAIN=", DOMAIN, ", MAHAL=", MAHAL,
            ", PROBIT=", PROBIT, ", AUC.weights=", AUC.weights, "), file=output.filename)", sep=""))
        if (is.null(presence.focal) == F) {
        if (n.ensembles > 1) {
            models.file <- paste(getwd(), "//models//", species.last, "_ENSEMBLE_", n.ensembles, "_models", sep="")
            models.file <- normalizePath(models.file, mustWork=F)
            if (file.exists(models.file) == T) {           
                assign("ensmodels.file", models.file, envir=.GlobalEnv)
                load(ensmodels.file)
                assign("ensmodels.focal", ensemble.models, envir=.GlobalEnv)          
                logger(paste("Focal ensemble models (object ensmodels.focal) loaded from: ", models.file, sep=""))
            }
        }else{
            models.file <- paste(getwd(), "//models//", species.last, "_models", sep="")
            models.file <- normalizePath(models.file, mustWork=F)
            if (file.exists(models.file) == T) {
                assign("ensmodels.file", models.file, envir=.GlobalEnv)
                load(ensmodels.file)
                assign("ensmodels.focal", ensemble.models, envir=.GlobalEnv)          
                logger(paste("Focal ensemble models (object ensmodels.focal) loaded from: ", models.file, sep=""))
            }
        }
        }
        putRcmdr("dialog.values", list())
        activateMenus()
	closeDialog()
	tkfocus(CommanderWindow())
        }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('ensemble.batch', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    capturebutton <- tkbutton(buttonsFrame, text="capture output in file", width="25", command=onCapture)
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(left1Frame, text="number of ensembles: ", width=20), n.ensembleEntry, sticky="w")
    tkgrid(tklabel(left1Frame, text="k.splits: ", width=20), k.splitsEntry, sticky="w")
    tkgrid(tklabel(left1Frame, text="ENSEMBLE.best: ", width=20), ENSEMBLE.bestEntry, sticky="w")
    tkgrid(tklabel(left1Frame, text="ENSEMBLE.min: ", width=20), ENSEMBLE.minEntry, sticky="w")
    tkgrid(tklabel(left1Frame, text="ENSEMBLE.exponent: ", width=20), ENSEMBLE.exponentEntry, sticky="w")
    tkgrid(PROBITCheckBox, tklabel(left2Frame, text="PROBIT transformations      "), sticky="w")
    tkgrid(AUC.weightsCheckBox, tklabel(left2Frame, text="Use mean AUC as weights"), sticky="w")
    tkgrid(tklabel(left2Frame, text="    ", width=5), tklabel(left2Frame, text="    ", width=30), sticky="w")
    tkgrid(left1Frame, sticky="w")
    tkgrid(left2Frame, sticky="w")
    tkgrid(tklabel(info1Frame, text="calibration: ", width=15), stack.focalEntry, sticky="w")
    tkgrid(tklabel(info1Frame, text="presence: ", width=15), presence.focalEntry, sticky="w")
    tkgrid(tklabel(info1Frame, text="absence: ", width=15), absence.focalEntry, sticky="w")
    tkgrid(tklabel(info2Frame, text="    ", width=35), absence.focalEntry, sticky="w")
    tkgrid(tklabel(info2Frame, text="(Better not to edit the parameters above", width=35), sticky="w")
    tkgrid(tklabel(info2Frame, text="but to use other menus to select the", width=35), sticky="w")
    tkgrid(tklabel(info2Frame, text="calibration stack, presence and absence)", width=35), sticky="w")
    tkgrid(info1Frame, sticky="w")
    tkgrid(info2Frame, sticky="w")
    tkgrid(leftFrame, infoFrame, sticky="w")
    tkgrid(firstFrame, sticky="w")
    tkgrid(tklabel(stackFrame, text="Select stacks to predict"), sticky="w")
    tkgrid(stackBox, stackScroll, sticky="w")
    tkgrid(tklabel(right1Frame, text="Select threshold method"), sticky="w")
    tkgrid(subsetBox, subsetScroll, sticky="w")
    tkgrid(varFrame, sticky="w")
    tkgrid(tklabel(right2Frame, text="threshold.sensitivity: ", width=20), sensitivityEntry, sticky="w")
    tkgrid(PRESABSCheckBox, tklabel(right3Frame, text="PresenceAbsence package"), sticky="w")
    tkgrid(right1Frame, sticky="w")
    tkgrid(right2Frame, sticky="w")
    tkgrid(right3Frame, sticky="w")
    tkgrid(stackFrame, rightFrame, sticky="w")
    tkgrid(secondFrame, sticky="w")
    tkgrid(tklabel(algoFrame, text="Select algorithms", width=15), sticky="w")
    tkgrid(MAXENTCheckBox, tklabel(algo1Frame, text="MAXENT"), sticky="w")
    tkgrid(GLMCheckBox, tklabel(algo1Frame, text="GLM"), sticky="w")
    tkgrid(MGCVCheckBox, tklabel(algo1Frame, text="MGCV"), sticky="w")
    tkgrid(NNETCheckBox, tklabel(algo1Frame, text="NNET"), sticky="w")
    tkgrid(BIOCLIMCheckBox, tklabel(algo1Frame, text="BIOCLIM"), sticky="w")
    tkgrid(tklabel(algo1Frame, text="    ", width=5), tklabel(algo1Frame, text="    ", width=12), sticky="w")
    tkgrid(GBMCheckBox, tklabel(algo2Frame, text="GBM"), sticky="w")
    tkgrid(GLMSTEPCheckBox, tklabel(algo2Frame, text="GLMSTEP"), sticky="w")
    tkgrid(MGCVFIXCheckBox, tklabel(algo2Frame, text="MGCVFIX"), sticky="w")
    tkgrid(FDACheckBox, tklabel(algo2Frame, text="FDA"), sticky="w")
    tkgrid(DOMAINCheckBox, tklabel(algo2Frame, text="DOMAIN"), sticky="w")
    tkgrid(tklabel(algo2Frame, text="    ", width=5), tklabel(algo2Frame, text="    ", width=12), sticky="w")
    tkgrid(GBMSTEPCheckBox, tklabel(algo3Frame, text="GBMSTEP"), sticky="w")
    tkgrid(GAMCheckBox, tklabel(algo3Frame, text="GAM"), sticky="w")
    tkgrid(EARTHCheckBox, tklabel(algo3Frame, text="EARTH"), sticky="w")
    tkgrid(SVMCheckBox, tklabel(algo3Frame, text="SVM"), sticky="w")
    tkgrid(MAHALCheckBox, tklabel(algo3Frame, text="MAHAL"), sticky="w")
    tkgrid(tklabel(algo3Frame, text="    ", width=5), tklabel(algo3Frame, text="    ", width=12), sticky="w")
    tkgrid(RFCheckBox, tklabel(algo4Frame, text="RF"), sticky="w")
    tkgrid(GAMSTEPCheckBox, tklabel(algo4Frame, text="GAMSTEP"), sticky="w")
    tkgrid(RPARTCheckBox, tklabel(algo4Frame, text="RPART"), sticky="w")
    tkgrid(SVMECheckBox, tklabel(algo4Frame, text="SVME"), sticky="w")
    tkgrid(tklabel(algo4Frame, text="     "), sticky="w")
    tkgrid(tklabel(algo4Frame, text="    ", width=5), tklabel(algo4Frame, text="    ", width=12), sticky="w")
    tkgrid(algo1Frame, algo2Frame, algo3Frame, algo4Frame, sticky="w")
    tkgrid(algoFrame, sticky="w")
    tkgrid(OKbutton, capturebutton, tklabel(buttonsFrame, text="     "), cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(subsetScroll, sticky="ns")
    tkgrid.configure(stackScroll, sticky="ns")
    tkselection.set(subsetBox, 0)
    tkselection.set(stackBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(subsetBox)
    tkwait.window(top)
}

ensemble.plot.GUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Plot suitability maps")
    firstFrame <- tkframe(top, relief="groove", borderwidth=2)
    left1Frame <- tkframe(firstFrame)
    right1Frame <- tkframe(firstFrame)
    species.names <- eval(parse(text=paste("levels(droplevels(factor(", presence.focal, "[, 1])))", sep="")), envir=.GlobalEnv)
    speciesFrame <- tkframe(left1Frame, relief="groove", borderwidth=2)
    speciesBox <- tklistbox(speciesFrame, width=40, height=6,
        selectmode="single", background="white", exportselection="FALSE") 
    speciesScroll <- tkscrollbar(speciesFrame, repeatinterval=5, command=function(...) tkyview(speciesBox, ...))
    tkconfigure(speciesBox, yscrollcommand=function(...) tkset(speciesScroll, ...))
    for (x in species.names) tkinsert(speciesBox, "end", x)
    m.variables <- c("suitability", "presence", "count")
    methodFrame <- tkframe(right1Frame, relief="groove", borderwidth=2)
    methodBox <- tklistbox(methodFrame, width=40, height=6,
        selectmode="single", background="white", exportselection="FALSE") 
    methodScroll <- tkscrollbar(methodFrame, repeatinterval=5, command=function(...) tkyview(methodBox, ...))
    tkconfigure(methodBox, yscrollcommand=function(...) tkset(methodScroll, ...))
    for (x in m.variables) tkinsert(methodBox, "end", x)
    secondFrame <- tkframe(top, relief="groove", borderwidth=2)
    left2Frame <- tkframe(secondFrame)
    right2Frame <- tkframe(secondFrame)
    update.stacklist()
    s.variables <- stack.list
    stackFrame <- tkframe(left2Frame, relief="groove", borderwidth=2)
    stackBox <- tklistbox(stackFrame, width=40, height=6,
        selectmode="single", background="white", exportselection="FALSE") 
    stackScroll <- tkscrollbar(stackFrame, repeatinterval=5, command=function(...) tkyview(stackBox, ...))
    tkconfigure(stackBox, yscrollcommand=function(...) tkset(stackScroll, ...))
    for (x in s.variables) tkinsert(stackBox, "end", x)
    variables <- c("spec_sens", "equal_sens_spec", "sensitivity", "no_omission", "prevalence", "kappa", "threshold.min", "threshold.mean", 
        "Sens=Spec", "MaxSens+Spec", "MaxKappa", "MaxPCC", "PredPrev=Obs", "ObsPrev", "MeanProb", "MinROCdist", "ReqSens")
    varFrame <- tkframe(right2Frame, relief="groove", borderwidth=2)
    subsetBox <- tklistbox(varFrame, width=40, height=6,
        selectmode="single", background="white", exportselection="FALSE")
    subsetScroll <- tkscrollbar(varFrame, repeatinterval=5, command=function(...) tkyview(subsetBox, ...))
    tkconfigure(subsetBox, yscrollcommand=function(...) tkset(subsetScroll, ...))
    for (x in variables) tkinsert(subsetBox, "end", x)
    thirdFrame <- tkframe(top, relief="groove", borderwidth=2)
    left3Frame <- tkframe(thirdFrame)
    right3Frame <- tkframe(thirdFrame)
    f.variables <- c("(none)", "_MEAN_", "_ENSEMBLE_", "_ENSEMBLE_1_", "_ENSEMBLE_2_", "_ENSEMBLE_3_", "_ENSEMBLE_4_", "_ENSEMBLE_5_",
        "_ENSEMBLE_6_", "_ENSEMBLE_7_", "_ENSEMBLE_8_", "_ENSEMBLE_9_", "_ENSEMBLE_10_")
    filterFrame <- tkframe(left3Frame, relief="groove", borderwidth=2)
    filterBox <- tklistbox(filterFrame, width=40, height=8,
        selectmode="single", background="white", exportselection="FALSE") 
    filterScroll <- tkscrollbar(filterFrame, repeatinterval=5, command=function(...) tkyview(filterBox, ...))
    tkconfigure(filterBox, yscrollcommand=function(...) tkset(filterScroll, ...))
    for (x in f.variables) tkinsert(filterBox, "end", x)
    right31Frame <- tkframe(right3Frame)
    right32Frame <- tkframe(right3Frame)
    right33Frame <- tkframe(right3Frame)
    right34Frame <- tkframe(right3Frame)
    fixedThresholdVariable <- tclVar("-1")
    fixedThresholdEntry <- tkentry(right31Frame, width=20, textvariable=fixedThresholdVariable)
    sensitivityVariable <- tclVar("0.9")
    sensitivityEntry <- tkentry(right31Frame, width=20, textvariable=sensitivityVariable)
    PRESABSVariable <- tclVar("0")
    PRESABSCheckBox <- tkcheckbutton(right32Frame, variable=PRESABSVariable)
    absVariable <- tclVar("6")
    absEntry <- tkentry(right33Frame, width=20, textvariable=absVariable)
    presVariable <- tclVar("6")
    presEntry <- tkentry(right33Frame, width=20, textvariable=presVariable)
    locVariable <- tclVar("0")
    locCheckBox <- tkcheckbutton(right34Frame, variable=locVariable)
    onOK <- function(){
        speciesValue <- species.names[as.numeric(tkcurselection(speciesBox))+1]
        stackValue <- s.variables[as.numeric(tkcurselection(stackBox))+1]
        methodValue <- m.variables[as.numeric(tkcurselection(methodBox))+1]
        filterValue <- f.variables[as.numeric(tkcurselection(filterBox))+1]
        var <- variables[as.numeric(tkcurselection(subsetBox))+1]
        sensitivity <- tclvalue(sensitivityVariable)
        fixedThreshold <- tclvalue(fixedThresholdVariable)
        PRESABS <- tclvalue(PRESABSVariable) == "1"
        if (var == "Sens=Spec"){PRESABS <- TRUE}
        if (var == "MaxSens+Spec"){PRESABS <- TRUE}
        if (var == "MaxKappa"){PRESABS <- TRUE}
        if (var == "MaxPCC"){PRESABS <- TRUE}
        if (var == "PredPrev=Obs"){PRESABS <- TRUE}
        if (var == "ObsPrev"){PRESABS <- TRUE}
        if (var == "MeanProb"){PRESABS <- TRUE}
        if (var == "MinROCdist"){PRESABS <- TRUE}
        if (var == "ReqSens"){PRESABS <- TRUE}
        presValue <- tclvalue(presVariable)
        absValue <- tclvalue(absVariable)
        if (filterValue == "(none)") {
            doItAndPrint(paste("ensemble.plot(RASTER.species.name='", speciesValue, "', RASTER.stack.name='", stackValue,
                "', plot.method='", methodValue, "', positive.filters=c('grd'), threshold=", fixedThreshold, ", p=", presence.focal, ", a=", absence.focal, 
                ", threshold.method='", var, "', threshold.sensitivity=", sensitivity, 
                ", threshold.PresenceAbsence=", PRESABS, 
                ", abs.breaks=", absValue, ", pres.breaks=", presValue, ")", sep=""))
        }else{
            doItAndPrint(paste("ensemble.plot(RASTER.species.name='", speciesValue, "', RASTER.stack.name='", stackValue,
                "', plot.method='", methodValue, "', positive.filters=c('grd', '", filterValue, 
                "'), threshold=", fixedThreshold, ", p=", presence.focal, ", a=", absence.focal, 
                ", threshold.method='", var, "', threshold.sensitivity=", sensitivity, 
                ", threshold.PresenceAbsence=", PRESABS, 
                ", abs.breaks=", absValue, ", pres.breaks=", presValue, ")", sep=""))
        }
        LOC <- tclvalue(locVariable) == "1"
        if (LOC == T) {doItAndPrint(paste("points(", presence.focal, "[which(", presence.focal, "[, 1] == '", 
            speciesValue, "'), c(2:3)], pch=4, col='black', cex=0.9, lwd=1.5)", sep=""))}
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('ensemble.plot', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(speciesFrame, text="Species variable"), sticky="w")
    tkgrid(speciesBox, speciesScroll, sticky="w")
    tkgrid(speciesFrame, sticky="w")
    tkgrid(tklabel(methodFrame, text="Suitability method"), sticky="w")
    tkgrid(methodBox, methodScroll, sticky="w")
    tkgrid(methodFrame, sticky="w")
    tkgrid(left1Frame, right1Frame, sticky="w")
    tkgrid(firstFrame, sticky="w")
    tkgrid(tklabel(stackFrame, text="Stack"), sticky="w")
    tkgrid(stackBox, stackScroll, sticky="w")
    tkgrid(stackFrame, sticky="w")
    tkgrid(tklabel(right2Frame, text="Select threshold method"), sticky="w")
    tkgrid(subsetBox, subsetScroll, sticky="w")
    tkgrid(varFrame, sticky="w")
    tkgrid(left2Frame, right2Frame, sticky="w")
    tkgrid(secondFrame, sticky="w")
    tkgrid(tklabel(filterFrame, text="Filter"), sticky="w")
    tkgrid(filterBox, filterScroll, sticky="w")
    tkgrid(filterFrame, sticky="w")
    tkgrid(tklabel(right31Frame, text="threshold value: ", width=20), fixedThresholdEntry, sticky="w")
    tkgrid(tklabel(right31Frame, text="threshold.sensitivity: ", width=20), sensitivityEntry, sticky="w")
    tkgrid(PRESABSCheckBox, tklabel(right32Frame, text="PresenceAbsence package"), sticky="w")
    tkgrid(tklabel(right33Frame, text="breaks for absence: ", width=20), absEntry, sticky="w")
    tkgrid(tklabel(right33Frame, text="breaks for presence: ", width=20), presEntry, sticky="w")
    tkgrid(locCheckBox, tklabel(right34Frame, text="Add presence locations"), sticky="w")
    tkgrid(right31Frame, sticky="w")
    tkgrid(right32Frame, sticky="w")
    tkgrid(right33Frame, sticky="w")
    tkgrid(right34Frame, sticky="w")
    tkgrid(left3Frame, right3Frame, sticky="w")
    tkgrid(thirdFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(speciesScroll, sticky="ns")
    tkgrid.configure(stackScroll, sticky="ns")
    tkgrid.configure(methodScroll, sticky="ns")
    tkgrid.configure(filterScroll, sticky="ns")
    tkgrid.configure(subsetScroll, sticky="ns")
    tkselection.set(speciesBox, 0)
    tkselection.set(stackBox, 0)
    tkselection.set(methodBox, 0)
    tkselection.set(filterBox, 0)
    tkselection.set(subsetBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(speciesBox)
    tkwait.window(top)
}


model.select.menu <- function(){
    logger(paste("Select models", sep=""))
    models.dir <- normalizePath(paste(getwd(), "//models//*models", sep=""), mustWork=F)
    assign("models.default", models.dir, envir=.GlobalEnv)
    doItAndPrint(paste("models.file <- choose.files(default=models.default, multi=F)", sep=""))
    models.file <- normalizePath(models.file, mustWork=F)
    if (file.exists(models.file) == T) {
        assign("ensmodels.file", models.file, envir=.GlobalEnv)
        load(ensmodels.file)
        assign("ensmodels.focal", ensemble.models, envir=.GlobalEnv)          
        logger(paste("Focal ensemble models (object ensmodels.focal) loaded from: ", models.file, sep=""))
    }
    putRcmdr("dialog.values", list())
    activateMenus()
}


eval.strip.GUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Evaluation strip")
    update.stacklist()
    s.variables <- stack.list
    stackFrame <- tkframe(top, relief="groove", borderwidth=2)
    stackBox <- tklistbox(stackFrame, width=65, height=8,
        selectmode="single", background="white", exportselection="FALSE") 
    stackScroll <- tkscrollbar(stackFrame, repeatinterval=5, command=function(...) tkyview(stackBox, ...))
    tkconfigure(stackBox, yscrollcommand=function(...) tkset(stackScroll, ...))
    for (x in s.variables) tkinsert(stackBox, "end", x)
    mvariables <- eval(parse(text=paste("names(ensmodels.focal$output.weights)[ensmodels.focal$output.weights>0]", sep="")), envir=.GlobalEnv)
    mvariables <- c("ENSEMBLE", mvariables)
    layerFrame <- tkframe(top, relief="groove", borderwidth=2)
    layerBox <- tklistbox(layerFrame, width=65, height=8,
        selectmode="single", background="white", exportselection="FALSE") 
    layerScroll <- tkscrollbar(layerFrame, repeatinterval=5, command=function(...) tkyview(layerBox, ...))
    tkconfigure(layerBox, yscrollcommand=function(...) tkset(layerScroll, ...))
    for (x in mvariables) tkinsert(layerBox, "end", x)
    onOK <- function(){
        stackValue <- s.variables[as.numeric(tkcurselection(stackBox))+1]
        command <- paste("evaluation.strip.data(xn=", stackValue, ", models.list=ensmodels.focal, steps=50)", sep="")
        logger(paste("strip.data <- ", command, sep=""))
        assign("strip.data", justDoIt(command), envir=.GlobalEnv)
    }
    onplot <- function(){
        layers <- mvariables[as.numeric(tkcurselection(layerBox))+1]
        doItAndPrint(paste("evaluation.strip.plot(data=strip.data, model='", layers, "')", sep=""))
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('evaluation.strip.data', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="prepare data", width="15", command=onOK, default="active")
    plotbutton <- tkbutton(buttonsFrame, text="plot", width="12", command=onplot, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(stackFrame, text="Select stack to prepare data"), sticky="w")
    tkgrid(stackBox, stackScroll, sticky="w")
    tkgrid(stackFrame, sticky="w")
    tkgrid(tklabel(layerFrame, text="Select model to plot"), sticky="w")
    tkgrid(layerBox, layerScroll, sticky="w")
    tkgrid(layerFrame, sticky="w")
    tkgrid(OKbutton, plotbutton, tklabel(buttonsFrame, text="    "), cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(layerScroll, sticky="ns")
    tkgrid.configure(stackScroll, sticky="ns")
    tkselection.set(layerBox, 0)
    tkselection.set(stackBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(layerBox)
    tkwait.window(top)
}

ens.predict.GUI <- function(){
    top <- tktoplevel()
    tkwm.title(top, "Predict ensemble models")
    update.stacklist()
    s.variables <- stack.list
    stackFrame <- tkframe(top, relief="groove", borderwidth=2)
    stackBox <- tklistbox(stackFrame, width=50, height=8,
        selectmode="single", background="white", exportselection="FALSE") 
    stackScroll <- tkscrollbar(stackFrame, repeatinterval=5, command=function(...) tkyview(stackBox, ...))
    tkconfigure(stackBox, yscrollcommand=function(...) tkset(stackScroll, ...))
    for (x in s.variables) tkinsert(stackBox, "end", x)
    onOK <- function(){
        stackValue <- s.variables[as.numeric(tkcurselection(stackBox))+1]
        doItAndPrint(paste("ensemble.raster(xn=", stackValue, ", models.list=ensmodels.focal, RASTER.species.name = ensmodels.focal$species.name)", sep=""))
    }
    onCancel <- function() {
        tkgrab.release(top)
        tkfocus(CommanderWindow())
        tkdestroy(top)  
        }
    buttonsFrame <- tkframe(top)
    onHelp <- function() {
        if (.Platform$OS.type != "windows") tkgrab.release(top)
        doItAndPrint(paste("help('ensemble.raster', help_type='html')", sep=""))
        }
    helpButton <- tkbutton(buttonsFrame, text="Help", width="12", command=onHelp)
    OKbutton <- tkbutton(buttonsFrame, text="OK", width="15", command=onOK, default="active")
    cancelButton <- tkbutton(buttonsFrame, text="Cancel", width="12", command=onCancel)
    tkgrid(tklabel(stackFrame, text="Select stack to predict suitabilities"), sticky="w")
    tkgrid(stackBox, stackScroll, sticky="w")
    tkgrid(stackFrame, sticky="w")
    tkgrid(OKbutton, cancelButton, helpButton)
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(stackScroll, sticky="ns")
    tkselection.set(stackBox, 0)
    for (row in 0:6) tkgrid.rowconfigure(top, row, weight=0)
    for (col in 0:0) tkgrid.columnconfigure(top, col, weight=0)
    .Tcl("update idletasks")
    tkwm.resizable(top, 0, 0)
    tkwm.deiconify(top)
    tkgrab.set(top)
    tkfocus(stackBox)
    tkwait.window(top)
}
