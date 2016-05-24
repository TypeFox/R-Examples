presence.absence.hist<-function(	DATA,
						which.model=1,
						na.rm=FALSE,
						xlab="predicted probability",
						ylab="number of plots",
						main=NULL,
						model.names=NULL,
						color=NULL,
						N.bars=20,
						truncate.tallest=FALSE,
						ylim=1.25*range(0,apply(counts,2,sum)),
						opt.thresholds=NULL,
						threshold=101,
						opt.methods=NULL,
						req.sens,
						req.spec,
						obs.prev=NULL,
						smoothing=1,
						add.legend=TRUE,
						legend.text=c("present","absent"),
						legend.cex=0.8,
						add.opt.legend=TRUE,
						opt.legend.text=NULL,
						opt.legend.cex=0.7,
						pch=NULL,
						FPC,FNC){
### Takes a single model and creates a stacked histogram of the predicted 
### probabilities with true presences and true absences in different colors
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
###	which.model 	a number indicating which model in DATA should be used for
###				calculating the confusion matrix
###	na.rm			should rows containing NA's be removed from the dataset
###				NOTE:	if ra.rm=FALSE, and NA's are present, 
###					function will return NA
###	main			main title
###	model.names		names of each model for the legend 
###	color			allows specification of colors for Presence/Absence
###				defaults to Presence = dark gray, Absence = light gray
###	N.bars		number of bin breaks for histograms
###	truncate.tallest	If one bar is much taller than others, should it be truncated to 
###				fit on plot
###	ylim			limit for y axis
###	opt.thresholds TRUE = calculate the optimal threshold by 4 meathods
###				if 'plot.it' is true, they are marked on the plot
###				if 'plot.it' is false, the values are returned along with their error
###				statistics.
###	threshold	cutoff values for translating predicted probabilities into
###			0 /1 values. (Only used if 'opt.thresholds'=TRUE.
###	It can be specified as either:
###		a single threshold (a number between 0 and 1)
###		a vector of thresholds (all between 0 and 1)
###		an interger representing the number of evenly spaced thresholds to calculate
###	
###	opt.methods		optimization methods for thresholds
###				Note: also controls which error statistics added to error.threshold.plot
###	req.sens		User defind required sensitivity
###				Note: 'req.sens' only used if opt.methods contain 'ReqSens'
###	req.spec		User defind required specificity
###				Note: 'req.spec' only used if opt.methods contain 'ReqSpec'
###	obs.prev	observed prevalence to be used in 'opt.meth="PredPrev=Obs" and "ObsPrev"'
###			defaults to observed prevalence from 'DATA'
###	smoothing		smoothing factor for finding optimized thresholds
###	add.legend		Should legend be added
###	legend.text		defaults to c("Present","Absent")
###	legend.cex		legend cex
###	add.opt.legend	TRUE = add opt legends to the plot
###	opt.legend.text	defaults to 'opt.methods'
###	opt.legend.cex	cex for opt legend
###	pch			point type for marking specific thresholds
	
### check logicals ###

if(is.logical(na.rm)==FALSE){
	stop("'na.rm' must be of logical type")}

if(is.logical(truncate.tallest)==FALSE){
	stop("'truncate.tallest' must be of logical type")}

if(is.logical(add.legend)==FALSE){
	stop("'add.legend' must be of logical type")}

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

### Check that 'which.model' is a single integer and not greater than number of models in DATA ###

	N.models<-ncol(DATA)-2

	if(length(which.model)!=1){
		stop("this function will only work for a single model, 'which.model' must be of length one")}
	if(which.model<1 || round(which.model)!=which.model){
		stop("'which.model' must be a positive integer")}
	if(which.model > N.models){
		stop("'which.model' must not be greater than number of models in 'DATA'")}

###check length of model.names matches number of models, if needed, generate model names ###

	if(is.null(model.names)==TRUE){
		model.names<-if(is.null(names(DATA))==FALSE){names(DATA)[-c(1,2)]}else{paste("Model",1:N.models)}
	}

	if(N.models!=length(model.names) && (length(which.model) != 1 || length(model.names) != 1)){
		stop(	"If 'model.names' is specified it must either be a single name, or a vector",
			"of the same length as the number of model predictions in 'DATA'")}	

### Pull out data from single model ###

	DATA<-DATA[,c(1,2,which.model+2)]
	model.names<-model.names[which.model]

### count number of thresholds ###
	
N.thr<-length(threshold)

###check that the 'threshold' is valid ###

if(min(threshold)<0){
	stop("'threshold' can not be negative")}

if(max(threshold)>1){
	if(N.thr==1 && round(threshold)==threshold){
		threshold<-seq(length=threshold,from=0,to=1)
		N.thr<-length(threshold)
	}else{
		stop("non-interger 'threshold' greater than 1")}
}

### Generate main title ###

	if(is.null(main)==TRUE){main<-model.names}

### Check legend text ###

	if(length(legend.text)!=2){
		stop("'legend.text' must be vector of length 2")}

	if(add.legend==FALSE){
		legend.text<-FALSE}

### Check colors ###

	if(is.null(color)==TRUE){
		color=c("gray20","gray80")}

	if(length(color)!=2){
		stop("color must be vector of length 2")}

### translate observations from values to presence/absence ###

	DATA[DATA[,2]>0,2]<-1

### check 'obs.prev'### 

if(is.null(obs.prev)==TRUE){
	obs.prev<-sum(DATA[,2])/nrow(DATA)}

if(obs.prev<0 || obs.prev>1){
	stop("'obs.prev' must be a number between zero and one")}
if(obs.prev==0){
	warning("because your observed prevalence was zero, results may be strange")}
if(obs.prev==1){
	warning("because your observed prevalence was one, results may be strange")}

### Generate Breaks ###

	breaks<-(0:N.bars)/N.bars

### Calculate bin counts for histogram ###

	counts<-rbind(	hist(DATA[DATA[,2]==1,3],breaks=breaks,plot=FALSE)$counts,
				hist(DATA[DATA[,2]==0,3],breaks=breaks,plot=FALSE)$counts)

	counts.stacked<-apply(counts,2,sum)
	MAX.BAR<-(rank(-counts.stacked, ties.method= "first")==1)
	SECOND.BAR<-(rank(-counts.stacked, ties.method= "first")==2)

### make plot ###
	spacer<-0.025
	xlim<-c(-spacer*N.bars,(1.2+spacer)*N.bars)

	if(truncate.tallest==TRUE && counts.stacked[MAX.BAR]>2*counts.stacked[SECOND.BAR] && counts.stacked[SECOND.BAR]!=0){
		print("height of tallest bar truncated to fit on plot")
		mids<-hist(DATA[DATA[,2]==1,3],breaks=breaks,plot=FALSE)$mids
		MAX.BAR.count<-counts[,MAX.BAR]
		MAX.BAR.stacked<-max(counts.stacked)
		MAX.BAR.trunc<-1.2*counts.stacked[SECOND.BAR]
		MAX.BAR.mid<-mids[MAX.BAR]

		new.counts<-counts[,MAX.BAR]*MAX.BAR.trunc/max(counts.stacked)
		counts[,MAX.BAR]<-new.counts

		hatching<-rep(0,N.bars)
		hatching[MAX.BAR]<-MAX.BAR.trunc
					
		xmids<-barplot(	counts,
					col=color,
					xlab=xlab, ylab=ylab, main=main,
					xlim=xlim,ylim=ylim)
		par(new=TRUE)
		barplot(	hatching,
				density=20, col="black",
				xlab="", ylab="", main="",
				xlim=xlim,ylim=ylim)
		text(xmids[MAX.BAR]+spacer*N.bars,1.05*MAX.BAR.trunc,(MAX.BAR.count[1]+MAX.BAR.count[2]))

		par(new=TRUE)
		xlim<-c(min(xlim)/max(xlim),1-min(xlim)/max(xlim))
		plot(	0:1,0:1,xlim=xlim,ylim=0:1,type="n",bty="n",xlab="",ylab="",yaxt="n")
		# add legend #
		inset<-c(.02,.02)
		if(add.legend==TRUE){
			leg.loc<-legend(	x="topright",inset=c(.02,.02),
						legend=legend.text,fill=color,cex=legend.cex)
			inset<-c(0.02,(1-(leg.loc$rect$top-leg.loc$rect$h))+0.05)}	

	}else{
		barplot(	counts,
				col=color,
				xlab=xlab, ylab=ylab,
				main=main,
				ylim=ylim)
		par(new=TRUE)
		xlim<-c(min(xlim)/max(xlim),1-min(xlim)/max(xlim))
		plot(	0:1,0:1,xlim=xlim,ylim=0:1,type="n",bty="n",xlab="",ylab="",yaxt="n")

		# add legend #
		inset<-c(.02,.02)
		if(add.legend==TRUE){
			leg.loc<-legend(	x="topright",inset=c(.02,.02),
						legend=legend.text,fill=color,cex=legend.cex,
						bg="white")
			inset<-c(0.02,(1-(leg.loc$rect$top-leg.loc$rect$h))+0.05)}
	}

### OPTIMAL THRESHOLDS ###

if(!is.null(opt.methods) && is.null(opt.thresholds)){opt.thresholds<-TRUE}
if(is.null(opt.methods) && is.null(opt.thresholds)){opt.thresholds<-FALSE}
if(is.null(opt.methods)){opt.methods<-c(1,2,4)}

if(is.logical(opt.thresholds)==TRUE){
if(opt.thresholds==TRUE){

	### Check optimization methods

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

	N.meth<-length(opt.methods)

	if(is.numeric(opt.methods)==TRUE){
		if(sum(opt.methods%in%(1:length(POSSIBLE.meth)))!=N.meth){
			stop("invalid optimization method")
		}else{
			opt.methods<-POSSIBLE.meth[opt.methods]}}

	if(sum(opt.methods%in%POSSIBLE.meth)!=N.meth){
		stop("invalid optimization method")}

	if(is.null(opt.legend.text)==TRUE){opt.legend.text<-opt.methods}
	if(length(opt.legend.text)!=N.meth){
		stop("'opt.legend.text' must be of same length as 'opt.methods'")}

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

	###calculate optimal thresholds###
	
	OBS<-DATA[,2]
	if(length(OBS[OBS==0])==0){
		stop(	"no observed absences in dataset, therefore specificity does not",
			"exist, and modeling, much less threshold optimization, is not very",
			"meaningful")}
	if(length(OBS[OBS==1])==0){
		stop(	"no observed presences in dataset, therefore sensitivity does not",
			"exist, and modeling, much less threshold optimization, is not very",
			"meaningful!")}
			
	###calculate optimal thresholds

	OPT.THRESH<-optimal.thresholds(	DATA,
							threshold=threshold,
							model.names=model.names,
							opt.methods=opt.methods,
							req.sens=req.sens,
							req.spec=req.spec,
							obs.prev=obs.prev,
							smoothing=smoothing,
							na.rm=na.rm,
							FPC=FPC,FNC=FNC)

	###set pch###

	if(is.null(pch)==TRUE){
		pch<-c(1,5,2,16,15,17,8,6,9,12,4,7)[match(opt.methods,POSSIBLE.meth)]
	}else{
		pch<-rep(pch,length(opt.methods))[1:length(opt.methods)]}
}}

### special case: specified thresholds ###

if(is.logical(opt.thresholds)==FALSE){

	if(!is.numeric(opt.thresholds)){
		stop("'opt.thresholds' must be 'TRUE', 'FALSE', or numeric")}

	if(min(opt.thresholds)<0){
		stop("'opt.thresholds' can not be negative")}

	if(max(opt.thresholds)>1){
		if(N.thr==1 && round(opt.thresholds)==opt.thresholds){
			opt.thresholds<-seq(length=opt.thresholds,from=0,to=1)
			N.thr<-length(opt.thresholds)
		}else{
			stop("non-interger, non-logical 'opt.thresholds' greater than 1")}
	}
	
	N.opt.thresh<-length(opt.thresholds)
	if(is.null(opt.legend.text)){opt.legend.text<-rep("threshold",N.opt.thresh)}
	if(length(opt.legend.text)!=N.opt.thresh){
		stop("length of 'opt.legend.text' does not match number of specified thresholds")}

	if(is.null(pch)){
		pch<-1:N.opt.thresh}
	if(length(pch)!=N.opt.thresh){
		stop("length of 'pch' does not match number of specified thresholds")}

	OPT.THRESH<-data.frame(opt.legend.text,opt.thresholds)
	opt.thresholds<-TRUE
}

###plot the optimal tresholds###	

if(opt.thresholds==TRUE){

	abline(v=OPT.THRESH[,2],lty=1,lwd=1,col="lightgray")
	points(OPT.THRESH[,2],rep(0.95,length(pch)),pch=pch,cex=2)

###add opt legend###

	if(add.opt.legend==TRUE){
		Mark.pretty<-round(OPT.THRESH[,2],2)
		Mark.pretty.char<-as.character(Mark.pretty)
		Mark.pretty.char[Mark.pretty==0]<-"0.00"
		Mark.pretty.char[Mark.pretty==1]<-"1.00"
		Mark.pretty.char[nchar(Mark.pretty.char)==3]<-paste(Mark.pretty.char[nchar(Mark.pretty.char)==3],"0",sep="")

		legend(	x="topright",
				inset=inset,
				legend=paste(Mark.pretty.char,"  ",opt.legend.text,sep=""),
				pch=pch,pt.cex=1,
				cex=opt.legend.cex,
				bg="white")}
}
}
