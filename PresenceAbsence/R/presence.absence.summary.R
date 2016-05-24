presence.absence.summary<-function(	DATA,
						threshold=101,
						find.auc=TRUE,
						which.model=1,
						na.rm=FALSE,
						main=NULL,
						model.names=NULL,
						alpha=0.05,
						N.bins=5,
						N.bars=10,
						truncate.tallest=FALSE,
						opt.thresholds=NULL,
						opt.methods=NULL,
						req.sens,
						req.spec,
						obs.prev=NULL,
						smoothing=1,
						vert.lines=FALSE,
						add.legend=TRUE,
						add.opt.legend=TRUE,
						legend.cex=0.6,
						opt.legend.cex=0.6,
						pch=NULL,
						FPC,FNC,
						cost.line=FALSE){
### Takes a single model and creates three types of plots. 
### It is not quite as flexible as the individual plot functions, some arguments
### are preset so that the plots will be comparable, but the remaining arguments
### have the same meaning.
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
###			It can be specified as either:
###				a single threshold (a number between 0 and 1)
###				a vector of thresholds (all between 0 and 1)
###				an interger representing the number of evenly spaced
###				thresholds to calculate
###
###	Note: to produce useful plots requires a large number of thresholds
###			
###   find.auc		should auc be calculated
###	which.model 	a number indicating which model in DATA should be used for
###				calculating the confusion matrix
###	req.sens		User defind required sensitivity
###				Note: 'req.sens' only used if opt.methods contain 'ReqSens'
###	req.spec		User defind required specificity
###				Note: 'req.spec' only used if opt.methods contain 'ReqSpec'
###	obs.prev		observed prevalence to be used in 'opt.meth="PredPrev=Obs" and "ObsPrev"'
###				defaults to observed prevalence from 'DATA'
###	smoothing		smoothing factor for finding optimized thresholds	
###	na.rm			should rows containing NA's be removed from the dataset
###				NOTE:	if ra.rm=FALSE, and NA's are present, 
###					function will return NA
###	main			main title
###	opt.thresholds	should optimal thresholds be calculated and plotted
###	opt.methods		optimization methods for thresholds
###				Note: also controls which error
###				statistics added to error.threshold.plot
###	model.names		names of each model for the legend 
###	alpha			two sided alpha value for the confidence intervals
###	N.bins		integer giving number of bins for predicted probabilities for calibration plot
###	truncate.tallest	should highest bar in histogram be truncated
###	add.legend		TRUE = add legends to the plot 
###	vert.lines 		should vertical lines be drawn at each optimized threshold
###	pch			point type for marking specific thresholds

### check logicals ###

if(is.logical(find.auc)==FALSE){
	stop("'find.auc' must be of logical type")}

if(is.logical(na.rm)==FALSE){
	stop("'na.rm' must be of logical type")}

if(is.logical(add.legend)==FALSE){
	stop("'add.legend' must be of logical type")}

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

if(!is.null(opt.methods) && is.null(opt.thresholds)){opt.thresholds<-TRUE}
if(is.null(opt.methods) && is.null(opt.thresholds)){opt.thresholds<-FALSE}
if(is.null(opt.methods)){opt.methods<-c(1,2,4)}

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
	if(FPC<=0 || FNC<=0){stop("costs must be positive")}
	if(is.logical(cost.line)==FALSE){
		stop("'cost.line' must be of logical type")}
	if(!"Cost"%in%opt.methods){cost.line<-FALSE}}

### checking 'smoothing'

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

### Check that 'which.model' is a single integer and not greater than number of models in DATA ###

	if(length(which.model)!=1){
		stop("this function will only work for a single model, 'which.model' must be of length one")}
	if(which.model<1 || round(which.model)!=which.model){
		stop("'which.model' must be a positive integer")}
	if(which.model+2 > ncol(DATA)){
		stop("'which.model' must not be greater than number of models in 'DATA'")}

###check length of model.names matches number of models, if needed, generate model names ###

	if(is.null(model.names)==TRUE){
		model.names<-if(is.null(names(DATA))==FALSE){names(DATA)[-c(1,2)]}else{paste("Model",1:(ncol(DATA)-2))}
	}

	if(ncol(DATA)-2!=length(model.names) && (length(which.model) != 1 || length(model.names) != 1)){
		stop(	"if 'model.names' is specified it must either be a single name, or a vector",
			"of the same length as the number of model predictions in 'DATA'")}	

### Pull out data from which.model model ###

	DATA<-DATA[,c(1,2,which.model+2)]
	if(length(model.names)!=1){model.names<-model.names[which.model]}

### translate observations from values to presence/absence ###

	DATA[DATA[,2]>0,2]<-1

### check 'obs.prev' ###

if(is.null(obs.prev)==TRUE){
	obs.prev<-sum(DATA[,2])/nrow(DATA)}

if(obs.prev<0 || obs.prev>1){
	stop("'obs.prev' must be a number between zero and one")}
if(obs.prev==0){
	warning("because your observed prevalence was zero, results may be strange")}
if(obs.prev==1){
	warning("because your observed prevalence was one, results may be strange")}

### Generate main title ###

	if(is.null(main)==TRUE){main<-paste("Accuracy Plots for",model.names)}

### set up plot region ###
	op1<-par(mfrow=c(2,2),cex=.7,oma=c(0,0,5,0))

### presence absence histogram ###
	presence.absence.hist(	DATA=DATA,
					N.bars=N.bars,
					truncate.tallest=truncate.tallest,
					main="",
					opt.thresholds=opt.thresholds,
					opt.methods=opt.methods,
					req.sens=req.sens,
					req.spec=req.spec,
					obs.prev=obs.prev,
					add.opt.legend=FALSE,
					legend.cex=legend.cex,
					FPC=FPC,FNC=FNC)

### calibration plot ###
	calibration.plot(	DATA=DATA,N.bins=N.bins,
			ylab="observed as proportion of bin",
			main="Observed vs. Predicted",
			alpha=alpha)
	
### roc plot ###
	auc.roc.plot(	DATA=DATA,
				threshold=threshold,
				find.auc=find.auc,
				main="ROC Plot",
				model.names=model.names,
				opt.thresholds=opt.thresholds,
				opt.methods=opt.methods, 
				req.sens=req.sens,
				req.spec=req.spec,
				obs.prev=obs.prev,
				smoothing=smoothing,
				line.type=FALSE,
				add.legend=add.legend,
				add.opt.legend=add.opt.legend,
				legend.cex=legend.cex,
				opt.legend.cex=opt.legend.cex,
				FPC=FPC,FNC=FNC,
				cost.line=cost.line)

### error threshold plot ###
	line.type<-c(2,3,1)
	er.thresh<-error.threshold.plot(	DATA=DATA,
							threshold=threshold,
							main="Error Rate verses Threshold",
							opt.thresholds=opt.thresholds,
							opt.methods=opt.methods,
							req.sens=req.sens,
							req.spec=req.spec,
							obs.prev=obs.prev,
							smoothing=smoothing,
							add.legend=add.legend,
							vert.lines=vert.lines,
							add.opt.legend=FALSE,
							legend.cex=legend.cex,
							FPC=FPC,FNC=FNC)$threshold

	
### main title ###
	mtext(main,side=3,cex=1.2,line=2,outer=TRUE)

### restore original parameters ###
	par(op1)
}
