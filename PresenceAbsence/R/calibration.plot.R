

calibration.plot<-function(	DATA,
				which.model=1,
				na.rm=FALSE,
				alpha=0.05,
				N.bins=5,
				xlab="Predicted Probability of Occurrence",
				ylab="Observed Occurrence as Proportion of Sites Surveyed",
				main=NULL,
				color=NULL,
				model.names=NULL){

### Takes a single model and creates a goodness of fit plot of Observed verses 
### predicted values. The predicted values are grouped into 10 bins, then the ratio
### of plots with observed values of present verses total plots is calculated for each bin. 
### The confidence interval for each bin is also ploted, along with the total number of 
### plots for the bin.
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
###	which.model a number indicating which model in DATA should be used for
###			calculating the confusion matrix
###	na.rm		should rows containing NA's be removed from the dataset
###			NOTE:	if ra.rm=FALSE, and NA's are present, 
###				function will return NA
###	alpha		two sided alpha value for the confidence intervals
###	N.bins		integer giving number of bins to split predicted probabilities into
###	color		logical value or color codes
###	model.names	names of each model for the legend

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

	N.models<-ncol(DATA)-2

	if(length(which.model)!=1){
		stop("this function will only work for a single model, 'which.model' must be of length one")}
	if(which.model<1 || round(which.model)!=which.model){
		stop("'which.model' must be a positive integer")}
	if(which.model > N.models){
		stop("'which.model' must not be greater than number of models in 'DATA'")}


###check length of model.names matches number of models, if needed, generate model names

	if(is.null(model.names)==TRUE){
		model.names<-if(is.null(names(DATA))==FALSE){names(DATA)[-c(1,2)]}else{paste("Model",1:N.models)}
	}

	if(ncol(DATA)-2!=length(model.names) && (length(which.model) != 1 || length(model.names) != 1)){
		stop(	"If 'model.names' is specified it must either be a single name, or a vector",
			"of the same length as the number of model predictions in 'DATA'")}	

###set colors, line widths###

	N.dat<-ncol(DATA)-2

	if(is.null(color)==TRUE){ 
		col<-1
		lwd=1
	}else{
		if(is.logical(color)==TRUE){
			if(color==FALSE){
				col<-1
				lwd=1
			}else{
				col<-((1:N.dat)+1)[which.model]
				lwd=2}
		}else{
			col<-rep(color,N.dat)[1:N.dat]
			col<-col[which.model]
			lwd=2}}

### Pull out data from single model

	DATA<-DATA[,c(1,2,which.model+2)]
	N.dat<-1
	if(length(model.names)!=1){model.names<-model.names[which.model]}

### Generate main title

	if(is.null(main)==TRUE){
		main<-paste("Observed vs. Predicted (",model.names,")",sep="")}

### translate observations from values to presence/absence

	DATA[DATA[,2]>0,2]<-1

###define the bins###
	bin.cuts<-(0:N.bins)/N.bins
	bin.centers<-((1:N.bins-.5))/N.bins

###split data into bins###
	DATA.breaks<-cut(DATA[,3],breaks=bin.cuts,labels=1:N.bins)

###find bin counts###
	N.total<-tapply(DATA[,2],DATA.breaks,length)
	N.presence<-tapply(DATA[,2],DATA.breaks,sum)

###find empty bins###
	Empty<-is.na(N.total)==TRUE

###remove NA's from the bin counts###
	N.total[is.na(N.total)==TRUE]<-0
	N.presence[is.na(N.presence)==TRUE]<-0

###find observed proportions for each bin###
	OBS.proportion<-N.presence/N.total
	OBS.proportion[Empty]<-NA

###find lower and upper bounds for the confidence intervals###
	df1.low<-2*(N.total-N.presence+1)
	df2.low<-2*N.presence

	df1.up<-2*(N.presence+1)
	df2.up<-2*(N.total-N.presence)

	Lower<-rep(0,N.bins)
	Upper<-rep(1,N.bins)

	TF<-N.presence!=0 #a true/false vector of bins that should have a lower bound greater than zero.
	Lower[TF]<-N.presence[TF]/
				(N.presence[TF]+((N.total[TF]-N.presence[TF]+1)*qf(1-alpha/2,df1.low[TF],df2.low[TF])))
	Lower[Empty]<-NA

	TF<-N.presence<N.total #a true/false vector of bins that should have an upper bound less than one.
	Upper[TF]<-((N.presence[TF]+1)*qf(1-alpha/2,df1.up[TF],df2.up[TF]))/
				(N.total[TF]-N.presence[TF]+((N.presence[TF]+1)*qf(1-alpha/2,df1.up[TF],df2.up[TF])))
	Upper[Empty]<-NA

###make plot###
	op<-par(pty="s")
	plot(c(-.05,1.05),c(-.05,1.05),type="n",xlab=xlab,ylab=ylab,main=main)
	abline(a=0,b=1,lty=3)
	for(i in 1:N.bins){
		lines(rep(bin.centers[i],2),c(Lower[i],Upper[i]),col=col,lwd=lwd)}
	points(bin.centers,OBS.proportion,pch=16,col=col,cex=1.5)
	text(bin.centers,Upper+.07,labels=N.total)

###restore original parameters###
	par(op)

###find average prediction per bin###

	BinPred<-tapply(DATA[,3],DATA.breaks,mean)


###create output dataframe###

	output<-data.frame(	BinCenter=bin.centers,
					NBin=N.total,
					BinObs=OBS.proportion,
					BinPred=BinPred,
					BinObsCIlower=Lower,
					BinObsCIupper=Upper)

	return(output)

}


