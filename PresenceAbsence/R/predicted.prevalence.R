
predicted.prevalence<-function(DATA,threshold=.5,which.model=(1:N.models),na.rm=FALSE){
### Calculates the observed and predicted prevalence for the species
### Function will work for one model and multiple thresholds, or one threshold
### and multiple models, or multiple models and multiple thresholds.
###
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
###		It can be specified as either:
###			a single threshold (a number between 0 and 1)
###			a vector of thresholds (all between 0 and 1)
###			an interger representing the number of evenly spaced thresholds to calculate
###
###	which.model a number or vector indicating which models in DATA should be used 
###	na.rm		should rows containing NA's be removed from the dataset
###			NOTE:	if ra.rm=FALSE, and NA's are present, 
###				function will return NA

### check logicals

if(is.logical(na.rm)==FALSE){
	stop("'na.rm' must be of logical type")}

### check for and deal with NA values:

	if(sum(is.na(DATA))>0){
     		if(na.rm==TRUE){
            	NA.rows<-apply(is.na(DATA),1,sum)
            	warning(	length(NA.rows[NA.rows>0]),
					" rows ignored due to NA values")
           		DATA<-DATA[NA.rows==0,]
      	}else{return(NA)}}

###translate actual observations from values to presence/absence###

	DATA[DATA[,2]>0,2]<-1

### Check that if 'which.model' is specified, it is an integer and not greater than number of models in DATA

	N.models<-ncol(DATA)-2

	if(min(which.model)<1 || sum(round(which.model)!=which.model)!=0){
		stop("values in 'which.model' must be positive integers")}
	if(max(which.model) > N.models){
		stop("values in 'which.model' must not be greater than number of models in 'DATA'")}

### Pull out data from 'which.model' model

	DATA<-DATA[,c(1,2,which.model+2)]

###check that length(threshold) matches number of models###

	N.thresh<-length(threshold)
	N.dat<-ncol(DATA)-2

	if(min(threshold)<0){
		stop("'threshold' can not be negative")}

	if(max(threshold)>1){
		if(N.thresh==1 && round(threshold)==threshold){
			threshold<-seq(length=threshold,from=0,to=1)
			N.thresh<-length(threshold)
		}else{
			stop("'threshold is a non-integer greater than 1")}
	}

###Calculate Observed Prevalence#####

	N.plots<-nrow(DATA)
	N.observed<-sum(DATA[,2])
	Prev.observed<-N.observed/N.plots

###Calculate Predicted Prevalence###

	PREVALENCE<-data.frame(matrix(0,N.thresh,N.dat+2))
	names(PREVALENCE)<-c(	"threshold",
				"Obs.Prevalence",
				if(is.null(names(DATA))==FALSE){names(DATA)[-c(1,2)]}else{paste("Model",1:N.models)[which.model]})

	PREVALENCE[,1]<-threshold
	PREVALENCE[,2]<-rep(Prev.observed,N.thresh)
	for(dat in 1:N.dat){
		PRED<-matrix(0,N.plots,N.thresh)
		for(thresh in 1:N.thresh){
			if(thresh==1){
				PRED[DATA[,dat+2]>=threshold[thresh],thresh]<-1
			}else{
				PRED[DATA[,dat+2]>threshold[thresh],thresh]<-1}}
		N.pred<-apply(PRED,2,sum)
		PREVALENCE[,dat+2]<-N.pred/N.plots
	}
	return(PREVALENCE)
}


