
roc.plot.calculate<-function(DATA,threshold=101,which.model=1,na.rm=FALSE){
### Calculates PCC, sensitivity, specificity, and Kappa for a single presence absence model 
### at a series of thresholds in preparation for creating a ROC plot.
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
###	It can be specified as either:
###		a single threshold (a number between 0 and 1)
###		a vector of thresholds (all between 0 and 1)
###		an interger representing the number of evenly spaced thresholds to calculate
###				
###	Note: to produce useful plots requires a large number of thresholds
###
###	which.model a number indicating which model in DATA should be used for
###			calculating the confusion matrix
###
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

### Check that 'which.model' is a single integer and not greater than number of models in DATA

	if(length(which.model)!=1){
		stop("this function will only work for a single model, 'which.model' must be of length one")}
	if(which.model<1 || round(which.model)!=which.model){
		stop("'which.model' must be a positive integer")}
	if(which.model+2 > ncol(DATA)){
		stop("'which.model' must not be greater than number of models in 'DATA'")}

### Pull out data from single model ###

	DATA<-DATA[,c(1,2,which.model+2)]

### count number of thresholds ###
	
	N.thr<-length(threshold)

###check that the 'threshold' is valid

	if(min(threshold)<0){
		stop("'threshold' can not be negative!")}

	if(max(threshold)>1){
		if(N.thr==1 && round(threshold)==threshold){
			threshold<-seq(length=threshold,from=0,to=1)
			N.thr<-length(threshold)
		}else{
			stop("non-interger 'threshold' greater than 1")}
	}
	
###set up error vectors###

	PCC      <-rep(0,N.thr)
	SENSITIVITY<-rep(0,N.thr)
	SPECIFICITY<-rep(0,N.thr)
	KAPPA      <-rep(0,N.thr)

###do calculations for each 'threshold'###

	for(thresh in 1:N.thr){
		CMX<-cmx(DATA=DATA,threshold=threshold[thresh])
		PCC[thresh]      <-pcc(CMX=CMX,st.dev=FALSE)
		SENSITIVITY[thresh]<-sensitivity(CMX=CMX,st.dev=FALSE)
		SPECIFICITY[thresh]<-specificity(CMX=CMX,st.dev=FALSE)
		KAPPA[thresh]      <-Kappa(CMX=CMX,st.dev=FALSE)
	}
	return(data.frame(threshold=threshold,PCC=PCC,sensitivity=SENSITIVITY,specificity=SPECIFICITY,Kappa=KAPPA))
}
