optimal.thresholds<-function(	DATA=NULL,
					threshold=101,
					which.model=1:(ncol(DATA)-2),
					model.names=NULL,
					na.rm=FALSE,
					opt.methods=NULL,
					req.sens,
					req.spec,
					obs.prev=NULL,
					smoothing=1,
					FPC,FNC){
### Finds the optimal thresholds by given methods

# opt.methods	optimization methods:
#	1  "Default"	threshold=0.5
#	2  "Sens=Spec"	threshold where sensitivity equals specificity
#	3  "MaxSens+Spec"	threshold that maximizes the sum of sensitivity and specificity
#	4  "MaxKappa"	threshold that maximizes Kappa
#	5  "MaxPCC"		threshold that maximizes PCC (percent correctly classified)
#	6  "PredPrev=Obs"	threshold where predicted prevalence equals observed prevalence
#	7  "ObsPrev"	threshold set to Observed prevalence
#	8  "MeanProb"	threshold set to mean predicted probability
#	9  "MinROCdist"	threshold where ROC plot makes closest approach to (0,1)
#	10 "ReqSens"	highest threshold where sensitivity meets user defined requirement
#	11 "ReqSpec"	lowest threshold where specificity meets user defined requirement
#	11 "Cost"		point where slope of cost function meets ROC plot

#
### DATA is a matrix (or dataframe) of observed and predicted values where: 
### 	the first column is the plot id, 
###   the second column is the observed values (either 0/1 or actual values),
###   the remaining columns are the predicted probabilities for the model.
###
###	DATA		matrix nrow=number of plots,
###				col1=PLOTID
###				col2=observed   (0 / 1)
###			  	col3=prediction probabilities from first model
###			  	col4=prediction probabilities from second model, etc...
###
###	threshold	cutoff values for translating predicted probabilities into
###			0 /1 values.
###	It can be specified as either:
###		a single threshold (a number between 0 and 1)
###		a vector of thresholds (all between 0 and 1)
###		an interger representing the number of evenly spaced thresholds to calculate
###				
###	Note: to produce useful plots requires a large number of thresholds
###
###	which.model a number indicating which model in DATA should be used for
###			calculating the confusion matrix
###	na.rm		should rows containing NA's be removed from the dataset
###			NOTE:	if ra.rm=FALSE, and NA's are present, 
###			function will return NA
###	model.names	names of each model for the legend
###	opt.methods	methods for optimizing thresholds
###	req.sens	User defind required sensitivity
###	req.spec	User defind required specificity
###	obs.prev	observed prevalence to be used in 'opt.meth="PredPrev=Obs" and "ObsPrev"'
###			defaults to observed prevalence from 'DATA'
###	smoothing	smoothing factor for finding optimized thresholds


### Check optimization methods ###

POSSIBLE.meth<-c(	"Default",
			"Sens=Spec",
			"MaxSens+Spec",
			"MaxKappa",
			"MaxPCC",
			"PredPrev=Obs",
			"ObsPrev",	
			"MeanProb",
			"MinROCdist",
			"ReqSens",
			"ReqSpec",
			"Cost")

### Special Case: 'DATA'=NULL , returns list of possible thresholds

if(is.null(DATA)==TRUE){
	return(POSSIBLE.meth)
}else{

### 'DATA!=NULL


### check opt.methods ###

if(is.null(opt.methods)==TRUE){
	opt.methods<-POSSIBLE.meth}

N.meth<-length(opt.methods)

if(is.numeric(opt.methods)==TRUE){
	if(sum(opt.methods%in%(1:length(POSSIBLE.meth)))!=N.meth){
		stop("invalid optimization method")
	}else{
		opt.methods<-POSSIBLE.meth[opt.methods]}}

	if(sum(opt.methods%in%POSSIBLE.meth)!=N.meth){
		stop("invalid optimization method")}

### check Req stuff ###
if("ReqSens"%in%opt.methods){
	if(missing(req.sens)){
		warning("req.sens defaults to 0.85")
		req.sens<-0.85}}
if("ReqSpec"%in%opt.methods){
	if(missing(req.spec)){
		warning("req.spec defaults to 0.85")
		req.spec<-0.85}}

### check cost benafit stuff ###
if("Cost"%in%opt.methods){
	if(missing(FPC) || missing(FNC)){
		warning("costs assumed to be equal")
		FPC<-1
		FNC<-1}
	if(FPC<=0 || FNC<=0){stop("costs must be positive")}}



### checking 'smoothing' ###

if(length(smoothing)!=1){
	stop("'smoothing' must be a single number greater than or equal to 1")
}else{
	if(is.numeric(smoothing)==FALSE){
		stop("'smoothing' must be a single number greater than or equal to 1")
	}else{
		if(smoothing<1){	
			stop("'smoothing' must be a single number greater than or equal to 1")}}}

### check for and deal with NA values: ###

	if(sum(is.na(DATA))>0){
     		if(na.rm==TRUE){
            	NA.rows<-apply(is.na(DATA),1,sum)
            	warning(	length(NA.rows[NA.rows>0]),
					" rows ignored due to NA values")
           		DATA<-DATA[NA.rows==0,]
      	}else{return(NA)}}

### translate actual observations from values to presence/absence ###

DATA[DATA[,2]>0,2]<-1
N.models<-ncol(DATA)-2

### check 'obs.prev' ###

if(is.null(obs.prev)==TRUE){
	obs.prev<-sum(DATA[,2])/nrow(DATA)}

if(obs.prev<0 || obs.prev>1){
	stop("'obs.prev' must be a number between zero and one")}
if(obs.prev==0){
	warning("because your observed prevalence was zero, results may be strange")}
if(obs.prev==1){
	warning("because your observed prevalence was one, results may be strange")}

### check that modelling makes sense (i.e. sensitivity and specificity exist) ###

OBS<-DATA[,2]
if(length(OBS[OBS==0])==0){
	stop(	"no observed absences in dataset, therefore specificity does not",
		"exist, and modeling, much less threshold optimization, is not very",
		"meaningful")}
if(length(OBS[OBS==1])==0){
	stop(	"no observed presences in dataset, therefore sensitivity does not",
		"exist, and modeling, much less threshold optimization, is not very",
		"meaningful")}

### Check that if 'which.model' is specified, it is an integer and not greater than number of models in DATA ###

if(min(which.model)<1 || sum(round(which.model)!=which.model)!=0){
	stop("values in 'which.model' must be positive integers!")}
if(max(which.model) > N.models){
	stop("values in 'which.model' must not be greater than number of models in 'DATA'")}

### check length of model.names matches number of models, if needed, generate model names ### 

if(is.null(model.names)==TRUE){
	model.names<-if(is.null(names(DATA))==FALSE){names(DATA)[-c(1,2)]}else{paste("Model",1:N.models)}
}

if(N.models!=length(model.names) && (length(which.model) != 1 || length(model.names) != 1)){
	stop(	"If 'model.names' is specified it must either be a single name, or a vector",
		"of the same length as the number of model predictions in 'DATA'")}	

### Pull out data from which model ###

DATA<-DATA[,c(1,2,which.model+2)]
if(length(model.names)!=1){model.names<-model.names[which.model]}

### find the number of models in DATA ###

N.dat<-ncol(DATA)-2

### count number of thresholds ###
	
N.thr<-length(threshold)

### check that the 'threshold' is valid ###

if(min(threshold)<0){
	stop("'threshold' can not be negative")}

if(max(threshold)>1){
	if(N.thr==1 && round(threshold)==threshold){
		threshold<-seq(length=threshold,from=0,to=1)
		N.thr<-length(threshold)
	}else{
		stop("non-interger 'threshold' greater than 1")}}

### set up return matrix ###

OPT.THRESH<-data.frame(matrix(0,N.meth,N.dat))
names(OPT.THRESH)<-model.names

### Calculate optimal thresholds: ###

for(dat in 1:N.dat){
	ACC<-presence.absence.accuracy(	DATA,
							which.model=dat,
							threshold=threshold,
							find.auc=FALSE,
							st.dev=FALSE)
	for(meth in 1:N.meth){

		if(opt.methods[meth]=="Default"){
			OPT.THRESH[meth,dat]<-0.5}

		if(opt.methods[meth]=="Sens=Spec"){
			SENS.SPEC<-abs(ACC$sensitivity-ACC$specificity)
			OPT.THRESH[meth,dat]<-mean(ACC$threshold[rank(SENS.SPEC,ties.method="min")<=smoothing])}

		if(opt.methods[meth]=="MaxSens+Spec"){
			SENS.SPEC<-ACC$sensitivity+ACC$specificity
			OPT.THRESH[meth,dat]<-mean(ACC$threshold[rank(-SENS.SPEC,ties.method="min")<=smoothing])}

		if(opt.methods[meth]=="MaxKappa"){
			OPT.THRESH[meth,dat]<-mean(ACC$threshold[rank(-ACC$Kappa,ties.method="min")<=smoothing])}

		if(opt.methods[meth]=="MaxPCC"){
			OPT.THRESH[meth,dat]<-mean(ACC$threshold[rank(-ACC$PCC,ties.method="min")<=smoothing])}
	
		if(opt.methods[meth]=="PredPrev=Obs"){
			PREV<-predicted.prevalence(	DATA,
								which.model=dat,
								threshold=threshold)
			PREV.diff<-abs(obs.prev-PREV[,3])
			OPT.THRESH[meth,dat]<-mean(PREV$threshold[rank(PREV.diff,ties.method="min")<=smoothing])}

		if(opt.methods[meth]=="ObsPrev"){
			OPT.THRESH[meth,dat]<-obs.prev}
		
		if(opt.methods[meth]=="MeanProb"){
			OPT.THRESH[meth,dat]<-mean(DATA[,2+dat])}

		if(opt.methods[meth]=="MinROCdist"){
			ROC<-(ACC$specificity-1)^2+(1-ACC$sensitivity)^2
			OPT.THRESH[meth,dat]<-mean(ACC$threshold[rank(ROC,ties.method="min")<=smoothing])}

		if(opt.methods[meth]=="ReqSens"){
			OPT.THRESH[meth,dat]<-max(ACC$threshold[(ACC$sensitivity)>=req.sens])}

		if(opt.methods[meth]=="ReqSpec"){
			OPT.THRESH[meth,dat]<-min(ACC$threshold[(ACC$specificity)>=req.spec])}

		if(opt.methods[meth]=="Cost"){
			if(obs.prev==0){obs.prev<-0.000001}
			sl<-(FPC/FNC)*(1-obs.prev)/obs.prev
			x<-(1-ACC$specificity)
			y<-ACC$sensitivity
			rad<-(x^2+y^2)^.5
			theta<-atan2(y,x)
			theta.sl<-atan(sl)
			theta.new<-theta-theta.sl
			x.new<-rad*cos(theta.new)
			y.new<-rad*sin(theta.new)
			OPT.THRESH[meth,dat]<-mean(ACC$threshold[rank(-y.new,ties.method="min")<=smoothing])}
	}
}
OPT.THRESH<-cbind(Method=opt.methods,OPT.THRESH)
return(OPT.THRESH)
}}
		