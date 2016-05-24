
error.threshold.plot<-function(	DATA,
						threshold=101,
						which.model=1,
						na.rm=FALSE,
						xlab="Threshold",
						ylab="Accuracy Measures",
						main=NULL,
						model.names=NULL,
						color=NULL,
						line.type=NULL,
						lwd=1,
						plot.it=TRUE,
						opt.thresholds=NULL,
						opt.methods=NULL,
						req.sens,
						req.spec,
						obs.prev=NULL,
						smoothing=1,
						vert.lines=FALSE,
						add.legend=TRUE,
						legend.text=legend.names,
						legend.cex=0.8,
						add.opt.legend=TRUE,
						opt.legend.text=NULL,
						opt.legend.cex=0.7,
						pch=NULL,
						FPC,FNC
						){

### Takes a single model and calculates sensitivity, specificity, and Kappa
### as a function of threshold. 'max.error' is the maximum allowable error for 
### points with positive observations.
### Note: this function is not elegant. The default for 'threshold' results in just calculating 
### a large number of evenly spaced thresholds, and picking the best from that list. If more
### than one threshold are tied, it takes the mean of the ties. This is good enough for graphs. 
### An elegant function would need to calculate the list of unique
### thresholds, the same length as the total number of plots, but not evenly spaced.
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
###	na.rm		should rows containing NA's be removed from the dataset
###			NOTE:	if ra.rm=FALSE, and NA's are present, 
###			function will return NA
###	xlab		label for x axis
###	ylab		label for y axis
###	main		main title
###	model.names	names of each model for the legend
###	color		each error statistic in a different color
###			color can be a vector of color codes, or alogical argument where:
###			TRUE  = each model in a different color
###			FALSE = each model in a different line style
###	line.type	line types for the graph
###			this can be a vector of line types, or a logical argument where:
###			TRUE  = each AUC in a different line type (default for black and white)
###			FALSE = each AUC as a solid line (default for color)	###	lwd		line width for the accuracy measures
###	plot.it	TRUE = A plot is generated
###			FALSE= calulations are caried out but no plot is generated
###	opt.thresholds TRUE = calculate the optimal threshold by 4 meathods
###			if 'plot.it' is true, they are marked on the plot
###			if 'plot.it' is false, the values are returned along with their error
###			statistics.
###	opt.methods	optimization methods for thresholds
###			See 'optimal.thresholds' for details
###			Note: also controls which error statistics added to error.threshold.plot
###	req.sens	User defind required sensitivity
###			Note: 'req.sens' only used if opt.methods contain 'ReqSens'
###	req.spec	User defind required specificity
###			Note: 'req.spec' only used if opt.methods contain 'ReqSpec'
###	obs.prev	observed prevalence to be used in 'opt.meth="PredPrev=Obs" and "ObsPrev"'
###			defaults to observed prevalence from 'DATA'
###	smoothing	smoothing factor for finding optimized thresholds
### 	vert.lines	should vertical lines be added to plot at each optimized threshold
###	add.legend		Should legend be added
###	legend.text		defaults to 'legend.names' for each accuracy statistic
###	legend.cex		legend cex
###	add.opt.legend	TRUE = add opt legends to the plot
###	opt.legend.text	defaults to 'opt.methods'
###	opt.legend.cex	cex for opt legend
###	pch		point type for marking specific thresholds

### check logicals ###

if(is.logical(na.rm)==FALSE){
	stop("'na.rm' must be of logical type")}

if(is.logical(plot.it)==FALSE){
	stop("'plot.it' must be of logical type")}

if(is.logical(add.legend)==FALSE){
	stop("'add.legend' must be of logical type")}

if(is.logical(add.opt.legend)==FALSE){
	stop("'add.opt.legend' must be of logical type")}

if(is.logical(vert.lines)==FALSE){
	stop("'vert.lines' must be of logical type")}

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

### Check that if 'which.model' is specified, it is an integer and not greater than number of models in DATA ###
if(min(which.model)<1 || sum(round(which.model)!=which.model)!=0){
	stop("values in 'which.model' must be positive integers")}
if(max(which.model) > N.models){
	stop("values in 'which.model' must not be greater than number of models in 'DATA'")}

###check length of model.names matches number of models, if needed, generate model names ### 

if(is.null(model.names)==TRUE){
	model.names<-if(is.null(names(DATA))==FALSE){names(DATA)[-c(1,2)]}else{paste("Model",1:N.models)}
}

if(N.models!=length(model.names) && (length(which.model) != 1 || length(model.names) != 1)){
	stop(	"If 'model.names' is specified it must either be a single name, or a vector",
		"of the same length as the number of model predictions in 'DATA'")}	

### Pull out data from single model ###

DATA<-DATA[,c(1,2,which.model+2)]
if(length(model.names)!=1){model.names<-model.names[which.model]}

### Generate main title ###

if(is.null(main)==TRUE){main<-model.names}

### checking 'smoothing' ###

if(length(smoothing)!=1){
	stop("'smoothing' must be a single number greater than or equal to 1")
}else{
	if(is.numeric(smoothing)==FALSE){
		stop("'smoothing' must be a single number greater than or equal to 1")
	}else{
		if(smoothing<1){	
			stop("'smoothing' must be a single number greater than or equal to 1")}}}

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

### set up colors and line types ###

N.lines<-2+sum(opt.methods %in% c("MaxKappa","MaxPCC","MinROCdist","MaxSens+Spec"))

if(is.null(color)==TRUE){ 
	colors<-rep(1,N.lines)
}else{
	if(is.logical(color)==TRUE){
		if(color==FALSE){
			colors<-rep(1,N.lines)
		}else{
			colors<-(1:N.lines)+1
			lwd<-2*lwd}
	}else{
		colors<-rep(color,N.lines)[1:N.lines]
		lwd<-2*lwd}}

if(is.null(line.type)==TRUE){ 
	line.type<-N.lines:1
}else{
	if(is.logical(line.type)==TRUE){
		if(line.type==FALSE){
			line.type<-rep(1,N.lines)
		}else{
			line.type<-N.lines:1}
	}else{
		line.type<-rep(line.type,N.lines)[1:N.lines]}}
			

### calculate error statistics ###

Model.dat<-roc.plot.calculate(DATA=DATA,threshold=threshold)

### plot lines for each error statistic ###

if(plot.it==TRUE){
	legend.names<-c("sensitivity","specificity")
	op<-par(pty="s")
	plot(c(0,1),c(0,1),type="n",xlab=xlab,ylab=ylab,main=main)
	lines(Model.dat$threshold,Model.dat$sensitivity,lty=line.type[1],lwd=lwd,col=colors[1])
	lines(Model.dat$threshold,Model.dat$specificity,lty=line.type[2],lwd=lwd,col=colors[2])
	flag<-2
	if("MaxSens+Spec" %in% opt.methods){
		flag<-flag+1
		lines(Model.dat$threshold,(Model.dat$sensitivity+Model.dat$specificity)-1,lty=line.type[flag],lwd=lwd,col=colors[flag])
		legend.names<-c(legend.names,"(Sens+Spec)-1")}
	if("MaxKappa" %in% opt.methods){
		flag<-flag+1
		lines(Model.dat$threshold,Model.dat$Kappa,lty=line.type[flag],lwd=lwd,col=colors[flag])
		legend.names<-c(legend.names,"Kappa")}
	if("MaxPCC" %in% opt.methods){
		flag<-flag+1
		lines(Model.dat$threshold,Model.dat$PCC,lty=line.type[flag],lwd=lwd,col=colors[flag])
		legend.names<-c(legend.names,"PCC")}
	if("MinROCdist" %in% opt.methods){
		flag<-flag+1
		lines(Model.dat$threshold,((Model.dat$specificity-1)^2+(1-Model.dat$sensitivity)^2)^.5,lty=line.type[flag],lwd=lwd,col=colors[flag])
		legend.names<-c(legend.names,"ROC Distance")}

	###restore original parameters###
	par(op)
}

### if 'opt.thresholds' equals true ###

if(is.logical(opt.thresholds)==TRUE){
if(opt.thresholds==TRUE){
	
	if(is.null(opt.legend.text)==TRUE){opt.legend.text<-opt.methods}
	if(length(opt.legend.text)!=N.meth){
		stop("'opt.legend.text' must be of same length as 'opt.methods'")}

### calculate optimal thresholds ###
	
	OBS<-DATA[,2]
	if(length(OBS[OBS==0])==0){
		stop(	"no observed absences in dataset, therefore specificity does not",
			"exist, and modeling, much less threshold optimization, is not very",
			"meaningful")}
	if(length(OBS[OBS==1])==0){
		stop("no observed presences in dataset, therefore sensitivity does not",
			"exist, and modeling, much less threshold optimization, is not very",
			"meaningful")}
			
	### check cost benafit stuff ###

	if(missing(FPC) || missing(FNC)){
		if("Cost"%in%opt.methods){
			warning("costs assumed to be equal")}
		FPC<-1
		FNC<-1}

	### calculate error statisics for optimal thresholds ###

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

	OPT.ERROR<-roc.plot.calculate(DATA=DATA,threshold=OPT.THRESH[,2])
	OPT.ERROR<-cbind(opt.methods,OPT.ERROR)

	### calculate y values for points ###
		
	OPT.Y<-rep(0,length(opt.methods))		
	OPT.Y[opt.methods %in% c("Default","Sens=Spec","PredPrev=Obs","ObsPrev","MeanProb","ReqSens","Cost")]<-
		OPT.ERROR[opt.methods %in% c("Default","Sens=Spec","PredPrev=Obs","ObsPrev","MeanProb","ReqSens","Cost"),]$sensitivity
	OPT.Y[opt.methods == "MaxSens+Spec"]<-
		(OPT.ERROR[opt.methods == "MaxSens+Spec",]$sensitivity+OPT.ERROR[opt.methods == "MaxSens+Spec",]$specificity)-1
	OPT.Y[opt.methods == "MaxKappa"]<-OPT.ERROR[opt.methods == "MaxKappa",]$Kappa
	OPT.Y[opt.methods == "ReqSpec"]<-OPT.ERROR[opt.methods == "ReqSpec",]$specificity
	OPT.Y[opt.methods == "MaxPCC"]<-OPT.ERROR[opt.methods == "MaxPCC",]$PCC
	OPT.Y[opt.methods == "MinROCdist"]<-
		((OPT.ERROR[opt.methods == "MinROCdist",]$specificity-1)^2+(1-OPT.ERROR[opt.methods == "MinROCdist",]$sensitivity)^2)^.5

	### set pch for points ###

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
	OPT.ERROR<-roc.plot.calculate(DATA=DATA,threshold=OPT.THRESH[,2])
	OPT.ERROR<-cbind(opt.legend.text,OPT.ERROR)
	opt.thresholds<-TRUE
	vert.lines<-TRUE
}

### plot the optimal tresholds ###
if(opt.thresholds==TRUE){	
	if(plot.it==TRUE){

		if("ReqSens" %in% opt.methods){
			abline(h=req.sens,lty=1,lwd=1,col="lightgray")}
		if("ReqSpec" %in% opt.methods){
			abline(h=req.spec,lty=1,lwd=1,col="lightgray")}

		if(vert.lines==TRUE){
			abline(v=OPT.THRESH[,2],lty=1,lwd=1,col="lightgray")
			OPT.Y<-rep(1,length(OPT.THRESH[,2]))}

		points(OPT.THRESH[,2],OPT.Y,pch=pch,cex=2)

	### add legends ###
		inset<-c(.02,.02)
		if(add.opt.legend==TRUE){

			Mark.pretty<-round(OPT.THRESH[,2],2)
			Mark.pretty.char<-as.character(Mark.pretty)
			Mark.pretty.char[Mark.pretty==0]<-"0.00"
			Mark.pretty.char[Mark.pretty==1]<-"1.00"
			Mark.pretty.char[nchar(Mark.pretty.char)==3]<-paste(Mark.pretty.char[nchar(Mark.pretty.char)==3],"0",sep="")

			leg.loc<-legend(	x="bottomright",inset=inset,
						legend=paste(Mark.pretty.char,"  ",opt.legend.text,sep=""),
						pch=pch,pt.cex=1,
						bg="white", cex=opt.legend.cex)
			inset<-c(0.02,leg.loc$rect$top+0.05)}
		if(add.legend==TRUE){
			if(length(legend.text)!=N.lines){
				stop("length of 'legend.text' must match number of lines drawn on plot")}
			legend(	x="bottomright",inset=inset,
					legend=legend.text,
					lty=line.type,cex=legend.cex,
					col=colors,lwd=lwd,bg="white" )}

	}
	return(OPT.ERROR)
}else{
	inset<-c(.02,.02)
	if(add.legend==TRUE & plot.it==TRUE){
		if(length(legend.text)!=N.lines){
			stop("length of 'legend.text' must match number of lines drawn on plot")}
		legend(	x="bottomright",inset=inset,
				legend=legend.text,
				lty=line.type,cex=legend.cex,
				col=colors,lwd=lwd,bg="white" )}
}
}

