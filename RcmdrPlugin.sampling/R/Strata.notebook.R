#### Notebook for drawing a stratified sample
#### Author: Susie Jentoft, Statistics Norway
#### Last edited: 16 October 2012


Strata.notebook<-function() {
	initializeDialog(title=gettextRcmdr("One-Stage Stratified Sample"))
	.activeDataSet <- ActiveDataSet()
	vars<-names(get(.activeDataSet))

#Introduction text
	tkgrid(tklabel(top, text="Please select the stratification variable, sample size and sample allocation.", justify="left"), columnspan=2, sticky="w")

#Strata Box
	stratabox <- variableListBox(top, vars , title = "Select the stratification variable: ")
	tkgrid(tklabel(top, text=""))
	tkgrid(getFrame(stratabox), columnspan=2, sticky="w")

# Total Sample Size
	sstotalframe<-tkframe(top)
	sstotal<-tclVar(50)
	sstotallabel<-tklabel(sstotalframe, text="The total sample size is:", fg="blue")
	ssField<-tkentry(sstotalframe, width="10", textvariable=sstotal)

	tkgrid(sstotallabel, ssField)
	tkgrid(tklabel(top, text=""))
	tkgrid(sstotalframe, columnspan=2, sticky="w")

#Strata distribution radio buttons
	buttonframe <-tkframe(top)
	radioButtons(window=buttonframe, name="dist", buttons=c("pro","min","fixed","opt"), 
		labels=c("Proportional to size","Proportional to size with a minimum of:","Fixed, equal size in each stratum","Neymans allocation using the variable:") , title="", right.buttons=FALSE)
	minframe <- tkframe(buttonframe)
	min<-tclVar(100)
	minField<-tkentry(minframe, width="5", textvariable=min)
	tkgrid(minField)
	tkblank(minframe)
	tkgrid(distFrame, minframe, sticky="w")

#Neymans variable box
	neybox<-variableListBox(top, vars, selectmode="single", title="")	
	tkgrid(tklabel(top, text=""))
	tkgrid(tklabel(top, text="How will the sample be distributed among strata?", fg="blue"), sticky="w", columnspan=2)
	tkgrid(buttonframe, sticky="w", columnspan=2)
	tkgrid(getFrame(neybox), sticky="w", columnspan=2)
	
#Sample variable # Take out for time being
#	nameframe <- tkframe(top)
#	today<-Sys.Date()
#	namedefault<- paste("Sample_",format(today, "%d%m%Y"), sep="")
#	nameVar <- tclVar(namedefault)
#   	nameEntry <- tkentry(nameframe, width="20", textvariable=nameVar)
#	tkgrid(tklabel(nameframe, text=gettextRcmdr("Indicator variable name: "), fg="blue"), nameEntry, sticky="w", columnspan=3)
#	tkgrid(tklabel(top, text="")) 
#	tkgrid(nameframe, sticky="w", columnspan=2) 

#Disable radio button function
	rbdisable<-function(){
		if (tclvalue(distVariable)=="pro"|tclvalue(distVariable)=="fixed") {
			tkconfigure(minField, state="disabled")
			tkconfigure(neybox$listbox, state="disabled")
			}
		if (tclvalue(distVariable)=="min") {
			tkconfigure(minField, state="normal")
			tkconfigure(neybox$listbox, state="disabled")
			}
		if (tclvalue(distVariable)=="opt") {
			tkconfigure(minField, state="disabled")
			tkconfigure(neybox$listbox, state="normal")
			}
		}
	tkbind(proButton, "<ButtonPress-1>", function()tclvalue(distVariable)<-"pro")
	tkbind(minButton, "<ButtonPress-1>", function()tclvalue(distVariable)<-"min")
	tkbind(fixedButton, "<ButtonPress-1>", function()tclvalue(distVariable)<-"fixed")
	tkbind(optButton, "<ButtonPress-1>", function()tclvalue(distVariable)<-"opt")
	tkbind(proButton, "<ButtonRelease-1>", function()rbdisable())
	tkbind(minButton, "<ButtonRelease-1>", function()rbdisable())
	tkbind(fixedButton, "<ButtonRelease-1>", function()rbdisable())
	tkbind(optButton, "<ButtonRelease-1>", function()rbdisable())
	rbdisable()

# on OK
	onOK<-function(){

		#From Main Variables
	 	strata1 <- as.numeric(match(getSelection(stratabox), vars))

		#from strata allocation
		sstotal<-as.numeric(tclvalue(sstotal)) #total sample size
		ppstype<-tclvalue(distVariable) #pro, min, fixed , opt
		nminimum<-as.numeric(tclvalue(min))
		neyvar<-as.numeric(match(getSelection(neybox), vars))
		nvar <- paste("Sample_", format(Sys.time(), "%d%m%Y"), "_", format(Sys.time(), "%H%M"), sep="")

	### CHECKS ###

		if (length(strata1)==0) {
 			tkmessageBox(message="A stratification variable is required") 
			return()}
		if (length(unique(get(.activeDataSet)[,strata1]))==1) {
 			tkmessageBox(message="Strata variable must have more than one stratum.") 
			return()}
		if (length(unique(get(.activeDataSet)[,strata1]))>=sstotal) {
 			tkmessageBox(message="There are more strata than the sample size allows. \nPlease change the stratification variable or increase the sample size.") 
			return()}
		if (length(neyvar)==0 & ppstype=="opt") {
			tkmessageBox(message="A variable is required") 
			return()}
		if (ppstype=="ney") {
			tkmessageBox(message="sorry this option is not available yet") 
			return()}

	### Print out function ###	wrong
		closeDialog()
		command <- paste(nvar, " <- StrataSample(data = ", .activeDataSet, ", strata1 = ", strata1, ", sstotal= ", sstotal,
			", ppstype = ", dQuote(ppstype), ", nminimum = ", nminimum, ", neyvar = ", neyvar, ")", sep="")
		doItAndPrint(command)
		activeDataSet(nvar)
		doItAndPrint(paste("View(", nvar, ")", sep=""))
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="strata", reset="Srata.notebook")
	tkgrid(tklabel(top, text="")) 
   	tkgrid(buttonsFrame, sticky="w", columnspan=2)
	dialogSuffix(rows=4, columns=2, focus=ssField)
}



### Stratification sample function ###
StrataSample <- function(data, strata1, sstotal, ppstype, nminimum, neyvar) {
	
	# Standardise
	data<-data[order(data[,strata1]),]
	data$strata.id<-as.numeric(data[,strata1])
	temp <- table(data$strata.id)[]

	# n strata calculations
		if (ppstype=="pro"){		
			n<-(temp/sum(temp))*sstotal
			}
		if (ppstype=="min"){
			n<-(temp/sum(temp))*sstotal
			n<-ifelse (n < nminimum, nminimum, n)
			}
		if (ppstype=="fixed"){
			nfixed <- sstotal/ length(temp)
			n<-rep(nfixed, length(temp))
			}
		if (ppstype=="opt"){
			tkmessageBox(message="neymans allocation has not been written yet") #need to add this in
			n <- 1
			}
		n <- special.round(n) #not sure if right to put round in here? Doesn't work without

	if (any(n==0)) {tkmessageBox(message="The sample size is too small to sample every strata") 
			Strata.notebook()}

### take sample - not sure about special rounding here...
	sampleunits <- sampling::strata(data=data, stratanames="strata.id", size=n, method="srswor") 
	sampleunits <-getdata(data, sampleunits)

### clean up
	sampleunits$strata.id <- sampleunits$Stratum <- sampleunits$ID_unit <-NULL
	return(sampleunits)
}