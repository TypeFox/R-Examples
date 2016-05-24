
presence.absence.accuracy<-function(DATA,threshold=.5,find.auc=TRUE,st.dev=TRUE,which.model=(1:(ncol(DATA)-2)),na.rm=FALSE){
### Calculates five accuracy measures for presence absence data and their standard
### deviations.
### Function will work for one model and multiple thresholds, or one threshold
### and multiple models, or multiple models each with their own threshold.
###
### Note: the standard errors are for comparing an individual model to random
###       assignment (i.e. AUC=.5), if you want to compare two models to each other  
###       it is necessary to account for correlation due to the fact that they 
###       use the same test set, which this function will not do.
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
###	find.auc	should auc be calculated
###	st.dev	should the standard deviation be calculated
###	which.model a number or vector indicating which models in DATA should be used 
###	na.rm		should rows containing NA's be removed from the dataset
###			NOTE:	if ra.rm=FALSE, and NA's are present, 
###				function will return NA

### check logicals

if(is.logical(find.auc)==FALSE){
	stop("'find.auc' must be of logical type")}

if(is.logical(st.dev)==FALSE){
	stop("'st.dev' must be of logical type")}

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
		if(min(which.model)<1 || sum(round(which.model)!=which.model)!=0){
			stop("values in 'which.model' must be positive integers")}
		if(max(which.model)+2 > ncol(DATA)){
			stop("values in 'which.model' must not be greater than number of models in 'DATA'")}

### Pull out data from which.model model

	DATA<-DATA[,c(1,2,which.model+2)]

### check for names ###

	N.models<-ncol(DATA)-2
	model.names<-if(is.null(names(DATA))){paste("Model",1:N.models)}else{names(DATA)[-c(1,2)]}

###check that length(threshold) matches number of models###
	
	N.thr<-length(threshold)
	N.dat<-ncol(DATA)-2
	REP.dat<-FALSE

	if(min(threshold)<0){
		stop("'threshold' can not be negative")}

	if(max(threshold)>1){
		if(N.thr==1 && round(threshold)==threshold && (N.dat==1 || N.dat==threshold)){
			threshold<-seq(length=threshold,from=0,to=1)
			N.thr<-length(threshold)
		}else{
			stop(	"either length of 'threshold' doesn't equal number of",
 				"models, or 'threshold is a non-integer greater than 1")}
	}

	if(N.thr==1 && N.dat>1){
		threshold<-rep(threshold,N.dat)
		N.thr<-length(threshold)}

	if(N.thr!=N.dat){
		if(N.dat==1){
			REP.dat<-TRUE
		}else{
			stop("length of 'threshold' doesn't equal number of models!")} 
	}
		
######Calculate Accuracy######
###Rep.dat=FALSE###

	if(REP.dat==FALSE){
		if(st.dev==FALSE){
			ERROR<-data.frame(matrix(0,N.dat,7))
			names(ERROR)<-c(	"model","threshold","PCC","sensitivity","specificity","Kappa","AUC")
			ERROR[,1]<-model.names
			ERROR[,2]<-threshold

			for(dat in 1:N.dat){
				CMX<-cmx(DATA=DATA,threshold=threshold[dat],which.model=dat)
				ERROR[dat,3]<-pcc(CMX=CMX,st.dev=FALSE)
				ERROR[dat,4]<-sensitivity(CMX=CMX,st.dev=FALSE)
				ERROR[dat,5]<-specificity(CMX=CMX,st.dev=FALSE)
				ERROR[dat,6]<-Kappa(CMX=CMX,st.dev=FALSE)
				if(find.auc==TRUE){
					ERROR[dat,7]<-auc(DATA=DATA,st.dev=FALSE,which.model=dat)}
			}

   		}else{

			ERROR<-data.frame(matrix(0,N.dat,12))
			names(ERROR)<-c(	"model","threshold",
						"PCC","sensitivity","specificity","Kappa","AUC",
						"PCC.sd","sensitivity.sd","specificity.sd","Kappa.sd","AUC.sd")
			ERROR[,1]<-model.names
			ERROR[,2]<-threshold

			for(dat in 1:N.dat){
				CMX<-cmx(DATA=DATA,threshold=threshold[dat],which.model=dat)
				ERROR[dat,c(3,8)]<-pcc(CMX)
				ERROR[dat,c(4,9)]<-sensitivity(CMX)
				ERROR[dat,c(5,10)]<-specificity(CMX)
				ERROR[dat,c(6,11)]<-Kappa(CMX)
				if(find.auc==TRUE){
					ERROR[dat,c(7,12)]<-auc(DATA=DATA,which.model=dat)}
			}
		}
###Rep.dat=TRUE###
	}else{
		if(st.dev==FALSE){
			ERROR<-data.frame(matrix(0,N.thr,7))
			names(ERROR)<-c(	"model","threshold","PCC","sensitivity","specificity","Kappa","AUC")
			ERROR[,1]<-model.names
			ERROR[,2]<-threshold

			for(thresh in 1:N.thr){
				CMX<-cmx(DATA=DATA,threshold=threshold[thresh])
				ERROR[thresh,3]<-pcc(CMX=CMX,st.dev=FALSE)
				ERROR[thresh,4]<-sensitivity(CMX=CMX,st.dev=FALSE)
				ERROR[thresh,5]<-specificity(CMX=CMX,st.dev=FALSE)
				ERROR[thresh,6]<-Kappa(CMX=CMX,st.dev=FALSE)
				}
			if(find.auc==TRUE){
				ERROR[,7]<-auc(DATA=DATA,st.dev=FALSE)}
		}else{

			ERROR<-data.frame(matrix(0,N.thr,12))
			names(ERROR)<-c(	"model","threshold",
						"PCC","sensitivity","specificity","Kappa","AUC",
						"PCC.sd","sensitivity.sd","specificity.sd","Kappa.sd","AUC.sd")
			ERROR[,1]<-model.names
			ERROR[,2]<-threshold

			for(thresh in 1:N.thr){
				CMX<-cmx(DATA=DATA,threshold=threshold[thresh])
				ERROR[thresh,c(3,8)]<-pcc(CMX)
				ERROR[thresh,c(4,9)]<-sensitivity(CMX)
				ERROR[thresh,c(5,10)]<-specificity(CMX)
				ERROR[thresh,c(6,11)]<-Kappa(CMX)
			}
			if(find.auc==TRUE){
				area<-auc(DATA)
				ERROR[,7]<-area$AUC
				ERROR[,12]<-area$AUC.sd}
		}
	}
	if(find.auc==TRUE){
		return(ERROR)
	}else{
		if(st.dev==FALSE){
			return(ERROR[,1:6])
		}else{
			return(ERROR[,c(1:6,8:11)])
		}
	}
}

