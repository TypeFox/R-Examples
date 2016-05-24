### Stage 1 sampling menu and function for a 2 stage design 
### Author: Susie Jentoft
### Last edited 22.10.2012
### Taken out notebook function, check digits and report


Stage1.notebook<-function() {
	one<-inclus<-samplingtypeFrame<-run1<-NULL
	initializeDialog(title="Stage 1 Sampling")
	variables2<-names(get(activeDataSet()))

# PSU Box
	PSUselection<-ifelse(is.na(match("PSU_SN", variables2)), NA, match("PSU_SN", variables2)-1)
	PSUBox<-variableListBox2(top, variableList=variables2, title="Select the PRIMARY SAMPLING UNIT variable:",
		initialSelection=PSUselection)
	tkgrid(getFrame(PSUBox), sticky="w")
	
# Strata Box
	strataselection<-ifelse(is.na(match("STRATUM_ID", variables2)), NA, match("STRATUM_ID", variables2)-1)
	strataBox<-variableListBox2(top, variableList=variables2, title="Select the STAGE 1 STRATIFICATION variable:",
		initialSelection=strataselection)
	tkblank(top)
	tkgrid(getFrame(strataBox), sticky="w")

# inclusion probability Box
	inclusselection<-ifelse(is.na(match("ST1_PROB", variables2)), NA, match("ST1_PROB", variables2)-1)
	inclusBox<-variableListBox2(top, variableList=variables2, title="Select the STAGE 1 INCLUSION PROBABILITY variable:",
		initialSelection=inclusselection)
	tkblank(top)
	tkgrid(getFrame(inclusBox), sticky="w")

# Sampling type variable
	radioButtons(top, name="samplingtype", buttons=c("Maxentropy","system"), values=c(1,2), 
		labels=c("Maximum entropy sampling (recommended) ","Systematic Sampling"), 
		title="What type of sampling do you want to do?",right.buttons=FALSE)
	tkblank(top)
	tkgrid(samplingtypeFrame, sticky="w")

# Randomise PSUs
	RandomValue<-tclVar("1")
	RandomCB<-tkcheckbutton(top, text="Randomize the PSU order within each stratum")
	tkconfigure(RandomCB, variable=RandomValue)
	tkblank(top)
	tkgrid(RandomCB, sticky="w")	

	onOK<-function(){
		stage1df<-activeDataSet()
		PSU <- as.double(tkcurselection(PSUBox$listbox))+1
		strata <- as.double(tkcurselection(strataBox$listbox))+1
		inclus<- as.double(tkcurselection(inclusBox$listbox))+1
		RandomPSU<-tclvalue(RandomValue)
		if (tclvalue(samplingtypeVariable)=="1") {samplingtype<-"UPmaxentropy"} else {samplingtype<-"UPsystematic"}

		if (length(PSU)==0) {check.fn2("Primary Sampling Unit")
			return()}
		if (length(strata)==0) {check.fn2("STRATIFICATION")
			return()}
		if (length(inclus)==0) {check.fn2("First Stage Inclusion Probability")
			return()}

		closeDialog()
		
		doItAndPrint(paste("Stage1Sample <- stage1Sample.fn(stage1df=",stage1df,", PSU=",PSU,", strata=",strata,", inclus=",inclus, ", samplingtype=",samplingtype,", randomPSU=", RandomPSU,")\n",sep=""))
		activeDataSet("Stage1Sample")
	}

	OKCancelHelp2(window=top, onHelp=Stage1onHelp)
	tkgrid(buttonsFrame, columnspan="2", sticky="w")	
	dialogSuffix(rows=3, columns=2, window=top, focus=buttonsFrame)
}



### Sample function for stage 1 sampling ###


stage1Sample.fn <-function(stage1df, PSU, strata, inclus, samplingtype, randomPSU){

#take out demographic variables
x<-stage1df
x<-x[order(x[,strata],x[,PSU]),]

y<-unique(x[,c(PSU,strata,inclus)])

PSU2<-match(names(x)[PSU], names(y))
strata2<-match(names(x)[strata],names(y))
inclus2<-match(names(x)[inclus],names(y))

y$strata.id<-as.numeric(y[,strata2])

#Set indicator I to 1 for those with inclusion probability =1
y$I<-ifelse(y[,inclus2]==1,1,0)

#select those with inclusion probability <1
sample<-y[which(y[,inclus2]<1),]
if (nrow(y)!=nrow(sample)) {
	message("Sampling has proceeded, however this dataset contains PSU(s) with inclusion probabilities = 1")
	}

#sample those with inclusion prob<1
total<-max(sample$strata.id)
pb<-tkProgressBar(title = "progress bar", min = 0,
                    max = total, width = 300)
selectionI<-matrix(NA, nrow=total, ncol=1)
for (i in 1: max(sample$strata.id)){
	selectionI[i]<-list(samplingtype(sample[,inclus2][which(sample$strata.id==i)]))
	setTkProgressBar(pb, i, label=paste( round(i/total*100, 0),"% done"))
	}
close(pb)
sample$I<-c(selectionI,recursive=TRUE)
sample<-rbind(subset(sample,sample$I==1),y[which(y$I==1),])

#Combine 
Selected.PSU<-x[which(x[,PSU] %in% sample[,PSU2]),]

#Randomize PSU order if check button selected
if (randomPSU==1){
	Selected.PSU$rand<-runif(nrow(Selected.PSU))
	rrr<-tapply(Selected.PSU$rand, FUN=min, INDEX=factor(Selected.PSU[,PSU]))
	Selected.PSU$rrr<-rrr[match(Selected.PSU[,PSU], rownames(rrr))]
	Selected.PSU<-Selected.PSU[order(Selected.PSU[,strata], Selected.PSU$rrr),]

}

#Tidy up 
row.names(Selected.PSU)<-Selected.PSU$I<-Selected.PSU$strata.id<-Selected.PSU$PSUorder<-Selected.PSU$rand<-Selected.PSU$rrr<-NULL
return(Selected.PSU)
}