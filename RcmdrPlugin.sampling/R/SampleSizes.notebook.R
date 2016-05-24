################################################################ 
### Sample size windows 
### Notebooks removed so it is compatible with Linux/Unix systems and easier to control
### Author: Susie Jentoft
### Last edited: 12.11.2012
### Notes: Sample size calculations includes 3 windows with 3 functions: Data sets, PSU data set options, Stratification options
###	
################################################################ 



### Window to select the two main data sets ###
#################################################################
Sample.datasets<-function(){
	initializeDialog(title="Sample Sizes - Data sets")
	variables<-listDataSets()

# PSU Box
	PSUBox<-variableListBox2(top, variableList=variables, title="Select the PRIMARY SAMPLING UNIT data set:")

# StrataBox with check button
	strataBox<-variableListBox2(top, variableList=variables, title="Select the STRATIFICATION data set:")
	stratayesno<-tclVar("0")
	strataCB<-tkcheckbutton(top, text="No stratification dataset",command=function()disable.fn2(stratayesno, strataBox))
	tkconfigure(strataCB, variable=stratayesno)

#forward function
	onOK<-function(){
		df1<-getSelection(PSUBox)
		stratayesno <-tclvalue(stratayesno) #0 is yes stratadataset, 1 is create dataset
		if (stratayesno==0) {df2<-getSelection(strataBox)} else {df2<-getSelection(PSUBox)}

		if (length(df1)==0) {check.fn2("PSU DATASET") 
			return()}

		if (stratayesno==0 & length(df2)==0) {check.fn2("STRATIFICATION DATASET") 
			return()}

		tkdestroy(top)
		PSU.options(df1=df1, df2=df2, stratayesno=stratayesno)
	}

#Grid
	tkgrid(tklabel(top, text="To calculate sample sizes, two data sets are required. \nOne should contain details on the Primary sampling units (PSUs) \nand the other on the stratification of these PSUs.", justify="left"), sticky="w")
	tkblank(top)
	tkgrid(getFrame(PSUBox), sticky="w")
	tkblank(top)
	tkgrid(getFrame(strataBox), sticky="w")
	tkgrid(strataCB, sticky="w")
	tkblank(top)

#help and buttons
	OKCancelHelp2(window=top, onHelp=dataonHelp)
	tkgrid(buttonsFrame, columnspan="2", sticky="w")
	dialogSuffix(rows=3, window=top, focus=buttonsFrame)
}





### Window to select PSU and SSU variable DETAILS ###
###################################################################

PSU.options <- function(df1, df2, stratayesno) {
	initializeDialog(title="Sample Sizes - PSU Data set")
	variables2<-names(get(df1))

	tkgrid(tklabel(top, text="Please select the variables to use in the PSU data set.                 ", justify="left"), sticky="w")
	tkblank(top)

#PSU List Box - PSU dataframe
	PSUBox<-variableListBox2(top, variableList=variables2, title="The PSU variable is:")
	tkgrid(getFrame(PSUBox), sticky="w")
	tkblank(top)

#Size list box - PSU dataframe
	sizeBox<-variableListBox2(top, variableList=variables2, title="The PSU SIZE variable is:")
	tkgrid(getFrame(sizeBox), sticky="w")
	tkblank(top)

#Strata Strata List Box
	strataBox<-variableListBox2(top, variableList=variables2, title="The STAGE 1 STRATIFICATION variable is:")
	tkgrid(getFrame(strataBox), sticky="w")
	tkblank(top)

#Domain options
	demoframe <- tkframe(top, relief="groove", borderwidth=2)
	radioButtons2(demoframe, name="demoyesno", buttons=c("no", "yes"), values=c(1,0),
		click.command=function()demodisable(),labels=c("no", "yes - the STAGE 2 stratification variables are: [select all]"), 
		title="Is there STAGE 2 STRATIFICATION?", right.buttons=FALSE)
	tkgrid(demoyesnoFrame, sticky="w")
	tkgrid(demoframe, sticky="w")
	demoBox<-variableListBox2(demoframe, variableList=variables2, selectmode="multiple", title="", listHeight=8)
	tkgrid(getFrame(demoBox), sticky="w")

	radioButtons2(demoframe, name="demoValue", buttons=c("prop2", "fixed2"), values=c(1,2),
		labels=c("Proportional to size","Equal sample size                                                        "), 
		title="How is the sample distributed among the STAGE 2 STRATA?        ", right.buttons=FALSE)
	tkblank(demoframe)
	tkgrid(demoValueFrame, sticky="w",columnspan=1)
	tkgrid(demoframe)

	demodisable<-function(){
		box.vars <- list(prop2Button, fixed2Button, prop2Label, fixed2Label, demoValueTitle)
		if (tclvalue(demoyesnoVariable)==1) {
			setstate2(box.vars, "disabled")
			tkconfigure(demoValueTitle, foreground="grey")
			setstate(demoBox, "disabled")
			}
		if (tclvalue(demoyesnoVariable)==0) {
			setstate2(box.vars, "normal")
			tkconfigure(demoValueTitle, foreground="blue")
			setstate(demoBox, "normal")
			}
	}
	demodisable()

# On ok
	onOK<-function(){
		PSU1 <- match(getSelection(PSUBox), variables2)
	 	strata1 <- match(getSelection(strataBox), variables2)
	 	size1 <- match(getSelection(sizeBox), variables2)
		demoyesno <- tclvalue(demoyesnoVariable) #0 for stage 2 strata, 1 for no stratification
		if (demoyesno==0) {demo1 <- match(getSelection(demoBox), variables2)} else {demo1 <- ""}
#		demo1 <- ifelse(demoyesno==0, match(getSelection(demoBox), variables2), "") #Doesn't work (with new box?)
		domainalloc <- tclvalue(demoValueVariable)

		if (length(strata1)==0) {check.fn2("STAGE 1 STRATIFICATION") 
			return()}
		if (length(PSU1)==0) {check.fn2("PSU")
			return()}
		if (length(size1)==0) {check.fn2("SIZE")
			return()}
		if (any(is.na(demo1))) {check.fn2("STAGE 2 STRATIFICATION") 
			return()}
		if (demoyesno==0) {
			if (!all(rowSums(get(df1)[,demo1])==get(df1)[,size1])) {
			calcsizes<-tkmessageBox(message = "Stage 2 strata sizes do not add up to the PSU size. \nWOuld you like to calcuate the PSU size based on the Stage 2 sizes?", icon = "warning", type = "yesnocancel", default = "yes")
			if (tclvalue(calcsizes)== "yes") {
				doItAndPrint(paste(df1,"$PSU_SIZE <- rowSums(", df1, "[, c(", list(demo1),")])", sep=""))
				size1 <- match("PSU_SIZE", names(get(df1)))
				}
			}}
		tkdestroy(top)
		strata.options(df1, df2, stratayesno, PSU1, strata1, size1, demoyesno, demo1, domainalloc)
		return()
	}

#help and buttons
	OKCancelHelp2(window=top, onHelp=PSUonHelp)
	tkgrid(buttonsFrame, columnspan="2", sticky="w")
	dialogSuffix(rows=3, window=top, focus=buttonsFrame)
	}


### Window for Stage 1 Stratification Options ###
###########################################################################

strata.options<-function(df1, df2, stratayesno, PSU1, strata1, size1, demoyesno, demo1, domainalloc ) {
	initializeDialog(title="Sample Sizes - Stage 1 Stratification")
	if (stratayesno==0) {
		variables<-names(get(df2))
		variables2 <- names(get(df1)) } else { 
		variables<-names(get(df1))
		variables2 <- names(get(df1))}

# autoSelection between PSU and strata Datasets for strata variable
	strata.selection<-NULL
	if (stratayesno==0) {
		strata.variable<-names(get(df1))[strata1]
		strata.selection<-ifelse(is.na(match(strata.variable, variables)), NA, match(strata.variable,variables)-1)
		}

#Strata List Box
	strataBox<-variableListBox2(top, variableList=variables, title="The Stage 1 STRATIFICATION variable is:",
		initialSelection = strata.selection)
	tkgrid(getFrame(strataBox), sticky="w")
	tkblank(top)

#mk listbox and CB
	mkframe <- tkframe(top, relief="groove", borderwidth=4)
	mkBoxframe <- tkframe(mkframe)
	mkBox<-variableListBox2(mkBoxframe, variableList=variables, title="The Stage 1 SAMPLE SIZE variable is:")
	mkValue<-tclVar("0")
	mkCB<-tkcheckbutton(mkBoxframe, text="Calculate STAGE 1 SAMPLE SIZES                  ",
		command=function()mk.disable())
	tkconfigure(mkCB, variable=mkValue)
	tkgrid(getFrame(mkBox), sticky="nw")
	tkgrid(mkCB, sticky ="w")

	costframe <- tkframe(mkframe, relief="groove", borderwidth=1)
	PSUcost <- tclVar("")
	Parcost <- tclVar("")
	PSUmean <- tclVar("")
	PSUV <- tclVar("")
	PSUcostBox <- ttkcombobox(costframe, values=variables, textvariable=PSUcost) 
	ParcostBox<- ttkcombobox(costframe, values=variables, textvariable=Parcost)
	PSUmeanBox <- ttkcombobox(costframe, values=variables2, textvariable=PSUmean)
	PSUVBox <- ttkcombobox(costframe, values=variables2, textvariable=PSUV)
	PSUcostLabel <- tklabel(costframe, text="The PSU COST variable is: ", fg="blue")
	ParcostLabel <- tklabel(costframe, text="The PARTICIPANT COST variable is: ", fg="blue")
	PSUmeanLabel <- tklabel(costframe, text="The PSU MEAN variable is:", fg="blue")
	PSUVLabel <- tklabel(costframe, text="The PSU VARIANCE variable is:", fg="blue")
	tkblank(costframe)
	tkgrid(PSUcostLabel, PSUcostBox, sticky="w")
	tkgrid(ParcostLabel, ParcostBox, sticky="w")
	tkgrid(PSUmeanLabel, PSUmeanBox, sticky="w")
	tkgrid(PSUVLabel, PSUVBox, sticky="w")
	tkblank(costframe)
	tkgrid(mkBoxframe, costframe, sticky="nw")
	tkgrid(mkframe, sticky="w")
	

#Strata Sizes
	ssframe <- tkframe(top, relief="groove", borderwidth=4)
	ssBoxframe <- tkframe(ssframe)
	ssBox<-variableListBox2(ssBoxframe, variableList=variables, title="The Stage 1 STRATA SAMPLE SIZE variable is:      ")
	ssValue<-tclVar("0")
	ssCB<-tkcheckbutton(ssBoxframe, text="Calculate STRATA SAMPLE SIZES",
		command=function()ss.disable())
	tkconfigure(ssCB, variable=ssValue)
	tkgrid(getFrame(ssBox), sticky="nw")
	tkgrid(ssCB, sticky ="w")
	
#Total Sample Size
	allocframe<-tkframe(ssframe, relief="groove", borderwidth=1)
	sstotalframe<-tkframe(allocframe)
	sstot<-tclVar(2000)
	ssField<-tkentry(sstotalframe, width="10", textvariable=sstot)
	sstotallabel<-tklabel(sstotalframe, text="The total sample size is approximately:", fg="blue")
	
#Allocation options
	radioButtons2(allocframe, name="ppsValue", buttons=c("rbpps1", "rbpps2", "rbpps3", "rbpps4"), values=c(1,2,3,4),
		labels=c("Proportional to size", "Proportional to size with a minimum in each stratum of: ", 
		"Fixed sample size per stratum", "Neymans allocation using variable:"), title="How is the sample distributed among the STRATA?",
		right.buttons=FALSE, click.command=function()psdisable())
	minsize<-tclVar(20)
	minField<-tkentry(allocframe, width="5", textvariable=50)
	neyBox<-variableListBox2(allocframe, variableList=variables, title="")

#Gridding
	tkgrid(sstotallabel, ssField, columnspan=1, sticky="w")
	tkgrid(sstotalframe, sticky="w")
	tkgrid(tklabel(allocframe, text=""))
	tkgrid(ppsValueFrame, minField, sticky="w",columnspan=1)
	tkgrid(getFrame(neyBox), sticky="w")
	tkgrid(ssBoxframe, allocframe, sticky="nw", columnspan=2)
	tkblank(top)
	tkgrid(ssframe, sticky="w", columnspan=2)

# Disabling
	mk.disable <- function() {
		boxes <- list(PSUcostBox, ParcostBox, PSUmeanBox, PSUVBox)
		labelnames <- list(PSUcostLabel, ParcostLabel, PSUmeanLabel, PSUVLabel)
		if (tclvalue(mkValue)==0) {
			setstate2(boxes, "disabled")
			setstate2(labelnames, "disabled")
			setstate(mkBox, "normal")
			}
		if (tclvalue(mkValue)==1) {
			setstate2(boxes, "readonly")
			setstate2(labelnames, "normal")
			setstate(mkBox, "disabled")
			}
	}

	psdisable<-function(){
		if (tclvalue(ppsValueVariable)==2) tkconfigure(minField,state="normal")
		if (tclvalue(ppsValueVariable)!=2) tkconfigure(minField,state="disabled")
		if (tclvalue(ppsValueVariable)==4) setstate(neyBox,state="normal")
		if (tclvalue(ppsValueVariable)!=4) setstate(neyBox,state="disabled")
	}
	
	ss.disable <- function() {
		boxes <- list(minField, rbpps1Button, rbpps2Button, rbpps3Button, rbpps4Button, 
		rbpps1Label, rbpps2Label, rbpps3Label, rbpps4Label, sstotallabel, ssField)
		if (tclvalue(ssValue)==0) {
			setstate2(boxes, "disabled")
			setstate(neyBox, "disabled")
			tkconfigure(ppsValueTitle, foreground="grey")
			setstate(ssBox, "normal")
			}
		if (tclvalue(ssValue)==1) {
			setstate2(boxes, "normal")
			setstate(neyBox, "normal")
			tkconfigure(ppsValueTitle, foreground="blue")
			setstate(ssBox, "disabled")
			}
		psdisable()
	}

	mk.disable()
	ss.disable()


#forward function #add in Ney condition check?
	onOK<-function(){
		strata2<-match(getSelection(strataBox), variables)
		mkvar<-match(getSelection(mkBox), variables)
		mkvalue<-tclvalue(mkValue) #0 is strataset contains mk value, 1 is optim
		ssyesno<-tclvalue(ssValue) # 1 for no sss variable, 0 for sss variable
		ssvar<-match(getSelection(ssBox), variables)
		sstotal<-as.numeric(tclvalue(sstot)) #total sample size
		ppstype<-tclvalue(ppsValueVariable) #1 for PPS, 2 for PPS with min, 3 for fixed
		nminimum<-as.numeric(tclvalue(minsize))
		neyvar <- match(getSelection(neyBox), variables)

		Parcost <- match(tclvalue(Parcost), variables)
		PSUcost <- match(tclvalue(PSUcost), variables)
		PSUV <- match(tclvalue(PSUV), variables2)
		PSUmean <- match(tclvalue(PSUmean), variables2)

		if (length(strata2)==0) {check.fn2("STRATIFICATION") 
			return()}
		if (length(mkvar)==0 & mkvalue==0) {check.fn2("Stage 1 sample size")
			return()}
		if (length(ssvar)==0 & ssyesno==0) {check.fn2("Strata sample size")
			return()}
		if (mkvalue==1 & is.na(Parcost)) {check.fn2("participant cost")
			return()}
		if (mkvalue==1 & is.na(PSUcost)) {check.fn2("PSU cost")
			return()}
		if (mkvalue==1 & is.na(PSUV)) {check.fn2("PSU variance")
			return()}
		if (mkvalue==1 & is.na(PSUmean)) {check.fn2("PSU mean")
			return()}
		if (any(is.na(match(levels(get(df1)[,strata1]), levels(get(df2)[,strata2]))))) {
			tkmessageBox(message = "Strata levels do not match between the data sets. Please check them.")
			return()}

### Add in report
		datetime <- paste(format(Sys.time(), "%d%m%Y"), "_", format(Sys.time(), "%H%M"), sep="")
		reportname <- paste("SampleSizes_", datetime, sep="")
		if (mkvalue==1) costtotal<-round(sum(as.double(SampleSizeDatasets$StrataSS$ST1_COST)))
		df1names<-names(get(df1))
		if (stratayesno==0) df2names<-names(get(df2))
		if (stratayesno==1) df2names<-names(get(df1))
		
		cat("Documentation for EHES sampling - Calculating Sample Sizes\n\nDate:",date(),"\nR Version:", paste(getRversion()),"\nRcmdrPlugin.sampling Version: ", packageDescription("RcmdrPlugin.sampling")$Version, "\nWorking Directory:", getwd(),
		"\n\nDatasets...\nThe dataset containing Primary Sampling Units was called:", df1,
		ifelse(stratayesno==0, "\nThe dataset containing stratification information was called:", "\nNo separate stratification dataset was used"), 
		if (stratayesno==0) df2,
		"\n\nPSU Dataset Details...\nThe Primary Sampling Unit variable was:", df1names[PSU1], "\nThe Size variable was: ", df1names[size1],
		if (demoyesno==0) c("\nThe second stage stratification (domain) variables were: ", df1names[demo1]),
		if (demoyesno==1) "\nNo second stage stratification (domain) variables were used",
		if (stratayesno==0) c("\n\nStratification Dataset Details...\nThe stratitfication variable was: ", df2names[strata2]), 
		if (mkvalue==0) "\nThe stage 1 sample size variable was called:",
		if (mkvalue==0) ifelse(stratayesno==0, df2names[mkvar], df1names[mkvar]),
		if (mkvalue==1) "\nOptimisation of stage 1 sample sizes was selected",
		if (ssyesno==0) "\nThe datset contained the strata sample sizes variable: ",
		if (ssyesno==0) ifelse(stratayesno==0, df2names[ssvar], df1names[ssvar]),
		if (ssyesno==1) c("\nThe dataset did not contain a strata sample size variable"),
		if (ssyesno==1 & ppstype==1) "\nThe strata sample sizes were selected to be distributed proportional to size",
		if (ssyesno==1 & ppstype==2) c("\nThe strata sample sizes were selected to be distributed proportional to size with a minumum of: ", nminimum),
		if (ssyesno==1 & ppstype==3) c("\nThe strata sample sizes were selected to be a fixed size."),
		if (ssyesno==1 & ppstype==4) c("\nThe strata sample sizes were selected to distributed using Neymans Allocation"),
		if (mkvalue==1) c("\nOptimal sample sizes were calculated using the PSU dataset mean and variance variables: ", df1names[PSUmean], df1names[PSUV]),
		if (mkvalue==1) c("\nThe cost was calculated using the variables: ", df2names[PSUcost], df2names[Parcost]),
		"\n\nAdditional Details...\nThe total sample size was: ", sstotal,
		if (mkvalue==1) c("\nThe total cost came to: ", costtotal),
		if (mkvalue==1 & demoyesno==1) c("\nThe weighted mean for the key variable is estimated as: ", signif(T.Var,digits=4), "\nThe overall variance for the key variable is estimated as: ", signif(TotalVar, digits=4)),
		if (mkvalue==1 & demoyesno==1) c("\nThe Coefficient of Variation (as a percentage) for the key variable is estimated as: ", signif(COEF, digits=4), "\nThe Design effect is estimated as: ", signif(Deff, digits=4)),
		"\n", file=reportname)

		command <- paste("SampleSizes(df1=", df1, ", df2=", df2, ", stratayesno=", stratayesno,  
			", PSU1=", PSU1, ", strata1=", strata1,", size1=", size1,", demo1=", c(demo1),", demoyesno=", demoyesno, ", domainalloc=", domainalloc, 
			", strata2=",strata2, ", mkvar=", mkvar,", mkvalue=", mkvalue,", ssyesno=", ssyesno, ", ssvar=", ssvar, ", sstotal=", sstotal, ", ppstype=", ppstype,
			", nminimum=", nminimum, ", Parcost=", Parcost, ", PSUcost=",PSUcost,
			", PSUV=",PSUV,", PSUmean=", PSUmean, ")\n", sep="")
		doItAndPrint(paste("SampleSizeDatasets <- ", command, sep=""))
		Strataname <- paste("StrataSizes_", datetime, sep="")
		PSUname <- paste("PSUSizes_", datetime, sep="")
		doItAndPrint(paste(Strataname, " <- SampleSizeDatasets$StrataSS", sep=""))
		doItAndPrint(paste(PSUname, " <- SampleSizeDatasets$PSUSS\n", sep=""))

#finish up
		tkdestroy(top)
		tkfocus(CommanderWindow())
		}
	

#help and buttons
	tkblank(top)
	OKCancelHelp2(window=top, onHelp=strataonHelp)
	tkgrid(buttonsFrame, columnspan=2, sticky="w")
	dialogSuffix(rows=3, window=top, focus=buttonsFrame)
	}


###############################################################################