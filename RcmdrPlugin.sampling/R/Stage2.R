### Notebook and function for stage 2 sampling
### Author: Susie Jentoft
### Last edited: 22.10.2012
### Taken out check digits and report

Stage2.notebook<-function() {
	top<-NULL
	initializeDialog(title="Sampling Stage 2", window=top)
	variables2<-names(get(activeDataSet()))

# PSU List Box
	PSUselection<-ifelse(is.na(match("PSU_SN", variables2)), NA, match("PSU_SN", variables2)-1)
	PSUBox<-variableListBox2(top, variableList=variables2, title="Select the PRIMARY SAMPLING UNIT variable:", 
		initialSelection=PSUselection)
	tkgrid(getFrame(PSUBox), sticky="w")

# Strata Variable
	strataselection<-ifelse(is.na(match("STRATUM_ID", variables2)), NA, match("STRATUM_ID", variables2)-1)
	strataBox<-variableListBox2(top, variableList=variables2, title="Select the STRATIFICATION variable:",
		initialSelection=strataselection)
	tkblank(top)
	tkgrid(getFrame(strataBox), sticky="w")

# Domain variable function
	demselection<-ifelse(is.na(match("DOMAIN_ID", variables2)), NA, match("DOMAIN_ID", variables2)-1)
	demBox<-variableListBox2(top, variableList=variables2, title="Select the STAGE 2 STRATIFICATION (domain) variable:",
		initialSelection=demselection)
	tkblank(top)
	tkgrid(getFrame(demBox), sticky="w")
	demoValue<-tclVar("0")
	demoCB<-tkcheckbutton(top, text="No second stage stratification Variables", command=function()disable.fn(demoValue, demBox$listbox, demBox$label))
	tkconfigure(demoCB, variable=demoValue)
	tkgrid(demoCB, sticky="w")

# Inclusion probability 2 variable
	inclus2selection<-ifelse(is.na(match("ST2_PROB", variables2)), NA, match("ST2_PROB", variables2)-1)
	inclus2Box<-variableListBox2(top, variableList=variables2, title="Select the STAGE 2 INCLUSION PROBABILITY variable:",
		initialSelection=inclus2selection)
	tkblank(top)
	tkgrid(getFrame(inclus2Box), sticky="w")

# on ok function
	onOK<-function(){
		df1 <- activeDataSet()
		PSU <- as.double(tkcurselection(PSUBox$listbox))+1
		strata <- as.double(tkcurselection(strataBox$listbox))+1
		dem <- as.double(tkcurselection(demBox$listbox))+1
		demyesno <- tclvalue(demoValue)
		inclus2 <- as.double(tkcurselection(inclus2Box$listbox))+1
		
	#checks
		if (length(PSU)==0) {check.fn2("PSU")
				return()}
		if (length(strata)==0) {check.fn2("STRATIFICATION")
				return()}
		if (length(inclus2)==0) {check.fn2("SECOND STAGE INCLUSION PROBABILITY")
				return()}
	closeDialog(top)
	command<-paste("SampleUnits<-stage2Sample.fn(stage2df=",df1,", PSU=",PSU,", strata=",strata,", dem=", dem, ", demyesno=", demyesno, ", inclus2=", inclus2, ")", sep="")
	doItAndPrint(command)	
	}

	OKCancelHelp2(window=top, onHelp=Stage2onHelp)
	tkgrid(buttonsFrame, columnspan="2", sticky="w")
	dialogSuffix(rows=2, columns=3, window=top, focus=buttonsFrame)
}




### SAMPLING FUNCTION ###

stage2Sample.fn<-function(stage2df, PSU, strata, dem, inclus1, inclus2) {
	stage2df[,PSU]<-factor(stage2df[,PSU])
	x<-stage2df

#Create random numbers
	x$rand<-runif(nrow(x))

#order by PSU, demo, random no.
	if (!is.na(dem)) {
		x<-x[order(x[,strata],x[,PSU],x[,dem],x$rand),]
		} else	{
			x<-x[order(x[,strata],x[,PSU],x$rand),]
			}
	
#calculate cumulative inclus probs.
	x$b<-cumsum(x[,inclus2])
	x$a[1]<-0
	x$a[2:nrow(x)]<-x$b[1:(nrow(x)-1)]

#calculate a random number between 0 and 1 to start
	r<-runif(1)

#calculate rounded down a value
	x$base<-floor(x$a)

#Add in modified r value
	x$r<-r+x$base

#select sample 
	x$I<-ifelse(x$r>=x$a & x$r<x$b, 1, 0)
	
#create frame with selected individs and tidy up

	data<-x[which(x$I==1),]
	
	data$inclusprob<-data[,inclus1]*data[,inclus2]
	data$sample.weights<-1/data$inclusprob
	data$b<-data$a<-data$I<-data$rand<-data$base<-data$r<-NULL

	if (!is.na(dem)) {data<-data[order(data[,PSU],data[,dem]),]} else {
		data<-data[order(data[,PSU]),]
		}

#Tidy up and save working file
	row.names(data)<-NULL
	data$strat.id<-NULL
	return(data)
	}