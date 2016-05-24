####################################################################################################################

model.interaction.plot <- 
function(	model.obj=NULL,
		x = NULL,					# the first variable to be plotted
		y = NULL,					# the second variable to be plotted
		response.category=NULL,			# for categorical response, which category to show in plot 
		quantiles=NULL,				# if QRF model, which quantile to show in plot
		all=FALSE,					# if QRF model,
		obs=1,					# if QRF model
		qdata.trainfn=NULL,			# if RF model was built outside ModelMap, need to supply training data 
		folder=NULL,				# No ending slash, to output to working dir = getwd()
		MODELfn=NULL,
		PLOTfn=NULL,
		pred.means = NULL,			# allows specification of values for other variables
		xlab = NULL,				# allows manual specification of the x label
		ylab = NULL,				# and y label
		x.range = NULL,				# manual range specification for the x variable
		y.range = NULL,				# and the y
		z.range = NULL,				# allows control of the vertical axis
		ticktype = "detailed",			# specifiy detailed types - otherwise "simple"
		theta = 55,           			# rotation 
		phi=40,               			# and elevation
		smooth = "none",      			# controls smoothing of the predicted surface
#		mask = FALSE,         			# controls masking using a sample intensity model
		plot.type = NULL,   			# controls whether a "persp" or "image" plot is drawn
	# Graphics Arguments
		device.type=NULL,	
		res=NULL,
		jpeg.res=72,
		device.width=7,
		device.height=7,
		units="in",
		pointsize=12,
		cex=par()$cex,
		col=NULL,
		xlim = NULL,
		ylim = NULL,
		zlim = NULL,
		...)                  			# allows the passing of additional arguments to plotting routine
                          	 			# useful options include shade, ltheta, lphi for controlling illumination
                          	 			# and cex for controlling text size - cex.axis and cex.lab have no effect
{
# this function is adapted from
# gbm.perspec version 2.9 April 2007
# J Leathwick/J Elith
#
# takes a rf or gbm regression tree object produced by model.build and
# plots a perspective plot showing predicted values for two predictors
# as specified by number using x and y
# values for all other variables are set at their mean by default
# but values can be specified by giving a list consisting of the variable name
# and its desired value, e.g., c(name1 = 12.2, name2 = 57.6)





########################################################################################
###################### GUI and checking for needed info ################################
########################################################################################


### Check Platform ###

Rplatform<-.Platform$OS.type


### Add to Filters Table ###


## Adds to file filters to Cran R Filters table.
if(.Platform$OS.type=="windows"){
	Filters<-rbind(Filters,csv=c("Comma-delimited files (*.csv)", "*.csv"))}

### check plot.type ###

## If the plot.type variable is NULL, then the user selects variable from pop-up list.
if (is.null(plot.type)){
	plot.type <- select.list(c("persp","image"), title="Select plot type.")
	if(plot.type=="" || is.null(plot.type)){
		plot.type=="persp"}	
}

plot.type<-switch(plot.type,	"persp"="persp",
					"p"="persp",
					"perspective"="persp",
					"image"="image",
					"i"="image",
					"image.plot"="image",
					"persp")

### Check Graphics Device Type ###

device.type<-check.device.type(device.type)

if(is.null(res)){res<-jpeg.res}

### Select Output Folder ###

if(is.null(folder)){
	if(any(device.type%in%c("jpeg","pdf","postscript"))){
		if(.Platform$OS.type=="windows"){
			folder<-choose.dir(default=getwd(), caption="Select directory")
		}else{
			folder<-getwd()}
	}
}

### check col.ramp ###

if(plot.type=="image"){
	 if( is.null(col)){
		l <- seq(100,0,length.out=101)
		c <- seq(0,100,length.out=101)
		col.ramp <- hcl(h = 120, c = c, l = l)
	}else{
		col.ramp <- col
	}
}

### check model.obj ###

if(is.null(model.obj)){
	if(is.null(MODELfn)){
		if(.Platform$OS.type=="windows"){
			MODELfn <- choose.files(caption="Select model object", filters = Filters["All",], multi = FALSE)
			if(is.null(MODELfn)){stop("must provide a model object")}
		}else{stop("must provide a model object")}
	}
	modelname<-load(MODELfn)
	if(length(modelname)!= 1){
		stop("file must contain single model object")}
	assign("model.obj",get(modelname))
}

### Extract Model Type from model.obj ###

model.type<-check.model.type(model.obj)

if(model.type=="QRF"){if("QRF"%in%names(model.obj)){model.obj<-model.obj$QRF}}

### Load modelling packages ###

if(model.type=="CF"){REQUIRE.party()}
if(model.type=="QRF"){REQUIRE.quantregForest()}
if(model.type=="SGB" || model.type=="QSGB"){REQUIRE.gbm()}

if(model.type=="CF"){
	WARN.IM<-TRUE
}else{
	WARN.IM<-FALSE
}

### Extract predList from model.obj ###

predList<-check.predList(model.obj=model.obj,model.type=model.type)
n.preds <- length(predList)
  
### Extract Response Type from model.obj ###

response.type<-check.response.type(model.obj=model.obj,model.type=model.type,ONEorTWO="model.obj")

### Extract Response Name from model.obj ###

if(model.type=="CF"){
	if(length(model.obj@responses@variables)==1){
		response.name <- names(model.obj@responses@variables)
	}else{
		stop("ModelMap does not support multivariate response")
	}
	#response.name <- as.character(model.obj@data@formula$response[2])
}else{
	if(!is.null(model.obj$response)){
#*#
		response.name<-model.obj$response
#*#
	}else{
		response.name<-""
	}
}

### check response.category ###

if(response.type=="categorical"){
	if(model.type=="RF"){ response.levels<-colnames(model.obj$votes)}
	if(model.type=="CF"){ response.levels<-model.obj@responses@levels[[1]]}
	if(model.type=="SGB"){response.levels<-model.obj$classes}

	if(is.null(response.category)){
		response.category <- select.list(response.levels, title="Select response category.")
		if(response.category=="" || is.null(response.category)){
			stop("must choose response category")
		}
	}
	if(!response.category%in%response.levels){
		stop("supplied 'response.category' of ",response.category," is not one of the response categories of 'model.obj'")}
	response.category<-as.character(response.category)
}	

### check quantiles ###

if(model.type=="QRF"){
	if(is.null(quantiles)){
		quantiles <- select.list(c("0.01","0.10","0.50","0.90","0.99"), title="Quantile", multiple = FALSE)	
		if(quantiles=="" || is.null(quantiles)){quantiles<-"0.50"}
		quantiles<-as.numeric(quantiles)
	}
	
	if(length(quantiles)>1){
		warning("'quantiles' is of length greater than 1, only first element will be used",immediate. = WARN.IM )
		quantiles<-quantiles[1]
	}
	
	if(!is.numeric(quantiles)){
		stop("'quantiles' must be a number between 0 and 1")}
	
	if(quantiles<=0 || quantiles>=1){
		stop("'quantiles' must be a number between 0 and 1")}
}

### define x labels ###

if (is.null(x)){
	x <- select.list(predList, title="Select predictor for X axis.")
	if(x=="" || is.null(x)){
		stop("must give predictor to use for X axis")}
}	

if(length(x)!=1){
	stop("'x' must be a single name or number")}

if(is.numeric(x)){
	if(x>=1 && x<=n.preds){
		x.name<-predList[x]
	}else{stop("if 'x' is numeric 'x' must be number between 1 and number of predictor variables in 'model.obj'") }
}else{
	if(x%in%predList){
		x.name<-x
		x<-match(x.name,predList)
	}else{stop("if 'x' is character string, 'x' must be name of a predictor variable in 'model.obj'") }
}

if(is.null(xlab)){xlab <- x.name} 

### define y labels ###

if (is.null(y)){
	y <- select.list(predList, title="Select predictor for Y axis.")
	if(y=="" || is.null(y)){
		stop("must give predictor to use for Y axis")}
}	

if(length(y)!=1){
	stop("'y' must be a single name or number")}
	
if(is.numeric(y)){
	if(y>=1 && y<=n.preds){
		y.name<-predList[y]
	}else{stop("if 'y' is numeric 'y' must be number between 1 and number of predictor variables in 'model.obj'") }
}else{
	if(y%in%predList){
		y.name<-y
		y<-match(y.name,predList)
	}else{stop("if 'y' is character string 'y' must be name of a predictor variable in 'model.obj'") }
}

if(is.null(ylab)){ylab <- y.name}

## extract var.type and var.levels

var.type<-NULL
var.levels<-NULL

#RF + QRF
if(model.type=="RF" || model.type=="QRF"){
	var.factors<-!sapply(model.obj$forest$xlevels,identical,0)
	var.type<-rep(0,n.preds)
	var.levels<-as.list(rep(0,n.preds))
	if( any(var.factors)){
		var.type[var.factors]<-sapply(model.obj$forest$xlevels[var.factors],length)
		var.levels[var.factors]<-model.obj$forest$xlevels[var.factors]
	}	
}

##CF
if(model.type=="CF"){
	inputs<-model.obj@data@get("input")
	var.factors<-sapply(inputs,is.factor)
	var.type<-rep(0,n.preds)
	var.levels<-as.list(rep(0,n.preds))
	if( any(var.factors)){
		var.type[var.factors]<-sapply(lapply(inputs,levels),length)[var.factors]
		var.levels[var.factors]<-lapply(inputs,levels)[var.factors]
	}	
}

##SGB
if(model.type=="SGB"){
	if(!is.null(model.obj$var.type)){var.type<-model.obj$var.type}
	if(!is.null(model.obj$var.levels)){var.levels<-model.obj$var.type}
}

### define predictor data ###

qdata<-NULL

if(model.type!="CF"){
	if(!is.null(model.obj$predictor.data)){
		qdata <- model.obj$predictor.data
	}else{
		if(model.type=="SGB" && !is.null(model.obj$data$x)){
			qdata<-model.obj$data$x
			Nrow<-length(qdata)/n.preds
			qdata<-data.frame(matrix(qdata,Nrow,n.preds))
			names(qdata)<-predList
	
			for(p in (1:n.preds)[var.type!=0]   ){
				qdata[,p]<-factor(qdata[,p],labels=var.levels[[p]])
			}
		}
	}
}

if(model.type=="CF"){
	qdata<-inputs
}
	
if(is.null(qdata)){

	## If training data is NULL, then the user selects file from pop-up browser.
	if (is.null(qdata.trainfn)){
		if(.Platform$OS.type=="windows"){
			qdata.trainfn <- choose.files(caption="Select data file", filters = Filters["csv",], multi = FALSE)
			if(is.null(qdata.trainfn)){stop("")}
		}else{stop("if RF model built outside of ModelMap package you must provide qdata.trainfn")}
	}
	
	## Check if file name is full path or basename
	if(is.matrix(qdata.trainfn)!=TRUE && is.data.frame(qdata.trainfn)!=TRUE){
		if(identical(basename(qdata.trainfn),qdata.trainfn)){
			qdata.trainfn<-file.path(folder,qdata.trainfn)}
	}
	
	## Read in training data
	if(is.matrix(qdata.trainfn)==TRUE || is.data.frame(qdata.trainfn)==TRUE){
		qdata<-qdata.trainfn
		qdata<-data.frame(qdata)
	}else{
		qdata<-read.table(file=qdata.trainfn,sep=",",header=TRUE,check.names=FALSE,as.is=TRUE)
	}

	## Extract predictor data
	qdata<-qdata[,predList]

	## apply levels to qdata
	for(p in (1:n.preds)[var.type!=0]   ){
		qdata[,p]<-factor(qdata[,p],labels=var.levels[[p]])
	}
}	


########################################################################################
############################ Actual Plot Setup Stuff ###################################
########################################################################################

### deal with xlim vs x.range ###

if(!is.null(xlim)){
	if(!is.null(x.range)){
		if(!identical(x.range,xlim)){stop("'x.range' and 'xlim' control the same thing, do not supply both arguments")}
	}else{
		x.range<-xlim
	}
}

if(!is.null(ylim)){
	if(!is.null(y.range)){
		if(!identical(y.range,ylim)){stop("'y.range' and 'ylim' control the same thing, do not supply both arguments")}
	}else{
		y.range<-ylim
	}
}

if(!is.null(zlim)){
	if(!is.null(z.range)){
		if(!identical(z.range,zlim)){stop("'z.range' and 'zlim' control the same thing, do not supply both arguments")}
	}else{
		z.range<-zlim
	}
}

### build grid of x and y values ###

if(var.type[x]==0){
	xaxt<-"s"
	N.x<-50
  	if (is.null(x.range)) {
    		x.var <- seq(min(qdata[,x],na.rm=TRUE),max(qdata[,x],na.rm=TRUE),length = N.x)
  	}else{x.var <- seq(x.range[1],x.range[2],length = N.x)}
}else{
	xaxt<-"n"
	N.x<-var.type[x]
	x.var <- as.factor(levels(qdata[,x]))
}

if(var.type[y]==0){
	yaxt<-"s"
	N.y<-50
  	if (is.null(y.range)) {
    		y.var <- seq(min(qdata[,y],na.rm=TRUE),max(qdata[,y],na.rm=TRUE),length = N.y)
  	}else{y.var <- seq(y.range[1],y.range[2],length = N.y)}
}else{
	yaxt<-"n"
	N.y<-var.type[y]
	y.var <- as.factor(levels(qdata[,y]))
}

pred.frame <- expand.grid(list(x.var,y.var))
names(pred.frame) <- c(x.name,y.name)

j <- 3
for (i in 1:n.preds) {
	if (i != x && i != y) {
		# find mean of continuous predictors
		if (var.type[i]==0) {
			m <- match(predList[i],names(pred.means))
			if (is.na(m)) {
				pred.frame[,j] <- mean(qdata[,i],na.rm=TRUE)
			}else{
				pred.frame[,j] <- pred.means[m]
			}
		}
		# find most common value of factored predictors #
		if (var.type[i]!=0) { 			
			m <- match(predList[i],names(pred.means))
			temp.table <- table(qdata[,i])
			temp.table <- sort(temp.table,decreasing = TRUE)
			if (is.na(m)) {
				pred.frame[,j] <- rep(names(temp.table)[1],N.x*N.y)
			}else{
				if(pred.means[m]%in%names(temp.table)){
					pred.frame[,j] <- pred.means[m]
				}else{
					stop("'pred.means' contains a value for factored predictor ",predList[i]," not found in training data")
				}
			}
			pred.frame[,j] <- factor(pred.frame[,j],levels=var.levels[[i]])
		}
		names(pred.frame)[j] <- predList[i]
		j <- j + 1
	}
}

pred.frame<-pred.frame[,predList]


### form the prediction ###

if(model.type=="RF"){
	if(response.type=="binary"){
		prediction<-predict(model.obj, pred.frame,type="vote")[,"1"]
	}

	if(response.type=="categorical"){
		prediction<-predict(model.obj, pred.frame,type="vote")[,response.category]
	}

	if(response.type=="continuous"){
		prediction<-predict(model.obj, pred.frame)
	}
}

if(model.type=="QRF"){

	prediction<-predict(model.obj, pred.frame, quantiles=quantiles, all=all, obs=obs)

}

if(model.type=="CF"){
	if(response.type=="binary"){
		prediction<-CF.list2df(predict(model.obj, newdata=pred.frame, OOB=FALSE, type="prob"))
		colnames(prediction)<-sapply(strsplit(colnames(prediction),paste(response.name,".",sep="")),'[',2)
		prediction<-prediction[,2]
	}

	if(response.type=="categorical"){
		prediction<-CF.list2df(predict(model.obj, newdata=pred.frame, OOB=FALSE, type="prob"))
		colnames(prediction)<-sapply(strsplit(colnames(prediction),paste(response.name,".",sep="")),'[',2)
		prediction<-prediction[,response.category]
	}

	if(response.type=="continuous"){
		prediction<-predict(model.obj, newdata=pred.frame, OOB=FALSE)
	}
}

if(model.type=="SGB"){
	if(is.null(model.obj$best.iter)){
		n.trees<-model.obj$n.trees
	}else{
		n.trees<-model.obj$best.iter
	}
	prediction <- gbm::predict.gbm(model.obj,pred.frame,n.trees = n.trees, type="response")
}


### model smooth if specified ###

if (smooth == "model") {
	
	if(any(is.factor(x.var),is.factor(y.var))){
		warning("smoothing is not appropriate for factored predictors",immediate. = WARN.IM)}

	pred.glm <- glm(prediction ~ ns(pred.frame[,1], df = 8) * ns(pred.frame[,2], df = 8), data=pred.frame,family=poisson)
	prediction <- fitted(pred.glm)
}

### report the maximum value and set up realistic ranges for z ###
 
max.pred <- max(prediction)
cat("maximum value = ",round(max.pred,2),"\n")

#print("STARTING Z.RANGE")
if (is.null(z.range)) {
	if (response.type == "binary" || response.type == "categorical") {
		z.range <- c(0,1)
	} 
	if(response.type == "poisson"){
      	z.range <- c(0,max.pred * 1.1)
    	}
	if(response.type == "continuous"){
		z.min <- min(prediction,na.rm=TRUE)
		z.max <- max(prediction,na.rm=TRUE)
		z.delta <- z.max - z.min
		if(plot.type=="persp"){
			z.range <- c(z.min - (1.1 * z.delta), z.max + (1.1 * z.delta))
		}else{
			z.range <- range(prediction,na.rm=TRUE)
		}
	}
}

if(z.range[2]==z.range[1]){z.range[2]<-z.range[2]+0.1}

### form the matrix ###
#print("starting pred.matrix")

	pred.matrix <- matrix(prediction,ncol=N.y,nrow=N.x)

### kernel smooth if specified ###

if (smooth == "average") {  #apply a 3 x 3 smoothing average
	if(any(is.factor(x.var),is.factor(y.var))){
		warning("smoothing is not appropriate for factored predictors",immediate. = WARN.IM)}

	pred.matrix.smooth <- pred.matrix
		for (i in 2:49) {
			for (j in 2:49) {
				pred.matrix.smooth[i,j] <- mean(pred.matrix[c((i-1):(i+1)),c((j-1):(j+1))])
			}
		}
	pred.matrix <- pred.matrix.smooth
}

### mask out values inside hyper-rectangle but outside of sample space ###

#if (mask) {
#	mask.trees <- mask.object$gbm.call$best.trees
#	point.prob <- gbm::predict.gbm(mask.object[[1]],pred.frame, n.trees = mask.trees, type="response")
#	point.prob <- matrix(point.prob,ncol=50,nrow=50)
#	pred.matrix[point.prob < 0.5] <- 0.0
#}

########################################################################################
###################################  Plot ##############################################
########################################################################################

### check output filename ###
#print("check output filename")
if(is.null(PLOTfn)){
	if(is.null(MODELfn)){
		if(response.type=="categorical"){
			PLOTfn<- paste(model.type,"_",response.type,"_",response.name,"_",response.category,"_",plot.type,"_",x.name,"_",y.name,sep="")
		}else{
			if(model.type=="QRF"){
				PLOTfn<- paste(model.type,"_",response.type,"_",response.name,"_",plot.type,"_",x.name,"_",y.name,"_",quantiles,"_quantile",sep="")
			}else{
				PLOTfn<- paste(model.type,"_",response.type,"_",response.name,"_",plot.type,"_",x.name,"_",y.name,sep="")
			}
		}
	}else{
		if(response.type=="categorical"){
			PLOTfn<- paste(MODELfn,"_",response.category,"_",plot.type,"_",x.name,"_",y.name,sep="")
		}else{
			if(model.type=="QRF"){
				PLOTfn<- paste(MODELfn,"_",plot.type,"_",x.name,"_",y.name,"_",quantiles,"_quantile",sep="")
			}else{
				PLOTfn<- paste(MODELfn,"_",plot.type,"_",x.name,"_",y.name,sep="")
			}
		}
	}
}

#print("PLOTfn:")
#print(PLOTfn)

if(identical(basename(PLOTfn),PLOTfn)){
	PLOTfn<-file.path(folder,PLOTfn)}


### loop thru devices ###
#print("starting loops")
#if(!"none"%in%device.type){
for(i in 1:length(device.type)){

#print(paste("Loop =" i))

### initialize graphics device ###

initialize.device(	PLOTfn=PLOTfn,DEVICE.TYPE=device.type[i],
				res=res,device.width=device.width,device.height=device.height,
				units=units,pointsize=pointsize,cex=cex)

###################################################################################

if (plot.type=="image") {

	zlab<-""
	if(response.type=="categorical"){
		zlab<-paste("probability of", response.category)}
	if(model.type=="QRF"){
		zlab<-paste("prediction for ", quantiles*100, "% quantile", sep="")}


	x.min<-min(unclass(x.var))
	x.max<-max(unclass(x.var))
	x.delta<-(x.max-x.min)
	x.trans<-x.delta/(2*(N.x-1))
	xlim=c((x.min-x.trans),(x.max+x.trans))

	y.min<-min(unclass(y.var))
	y.max<-max(unclass(y.var))
	y.delta<-(y.max-y.min)
	y.trans<-y.delta/(2*(N.y-1))
	ylim=c((y.min-y.trans),(y.max+y.trans))

	image.plot(	x = unclass(x.var), y = unclass(y.var), 
		z = pred.matrix, zlim = z.range, 
		xlab = xlab, ylab = ylab,
		legend.lab=zlab, 
		legend.line=2.5,
		xaxs="i",yaxs="i",
		xlim=xlim,ylim=ylim,
		xaxt=xaxt,yaxt=yaxt,
		col=col.ramp,
		...)
	
	
	if(xaxt=="n"){
		par(xpd=TRUE)
		x.loc<-unclass(x.var)
		offset1<-y.delta*.02       		
		offset2<-y.delta*.04			
		text(x=x.loc,y=ylim[1]-offset2,labels=x.var,pos=1, offset=0,xpd=TRUE)
		mtext("(categorical)",side=1,line=1.8)
		for(j in 1:length(x.loc)){lines(c(x.loc[j],x.loc[j]),c(ylim[1],ylim[1]-offset1    ) )  }
		par(xpd=FALSE)
	}

 	if(yaxt=="n"){
		par(xpd=TRUE)
		y.loc<-unclass(y.var)
		offset1<-x.delta*.02
		offset2<-x.delta*.04
		text(x=xlim[1]-offset2,y=unclass(y.var),labels=y.var,pos=2, offset=0,xpd=TRUE)
		mtext("(categorical)",side=2,line=2)
		for(j in 1:length(y.loc)){lines(c(xlim[1],xlim[1]-offset1),c(y.loc[j],y.loc[j])) }
		par(xpd=FALSE)
	}
	#par(xpd=TRUE)
	#mtext(zlab,side=4,line=4.5,outer=TRUE)
	#par(xpd=FALSE)
}
if(plot.type=="persp"){

	zlab<-"fitted value"
	if(response.type=="categorical"){
		zlab<-paste("probability of", response.category)}
	if(model.type=="QRF"){
		zlab<-paste("prediction for ", quantiles*100, "% quantile", sep="")}

    	if(!any(is.factor(x.var),is.factor(y.var))){
		xlab<-paste("\n\n",xlab,sep="")
		ylab<-paste("\n\n",ylab,sep="")
		zlab<-paste("\n\n",zlab,sep="")
    		persp(x=unclass(x.var), y=unclass(y.var), z=pred.matrix, zlim= z.range,      	# input vars
      		xlab = xlab, ylab = ylab, zlab = zlab,   						# labels
      		theta=theta, phi=phi, r = sqrt(10), d = 3,               			# viewing pars
      		ticktype = ticktype, mgp = c(4,1,0), ...)   
    
	}else{
	
		### calculate axis values ###
		
		a=1
		if(theta>=0 && theta<=90){a=1}
		if(theta>90 && theta<=180){a=2}
		if(theta>180 && theta<=270){a=3}
		if(theta>270 && theta<=360){a=4}

		PN<-data.frame(	X=c(1,1,-1,-1),Y=c(-1,1,1,-1),  #PositiveNegative
					ZX=c(-1,1,1,-1),ZY=c(-1,-1,1,1)) 
		MM<-(PN/2)+1.5 #Max/Min

		if(is.factor(x.var)){
			x.loc<-unclass(x.var)
			xlim<-range(x.loc)
			x.axis<-levels(x.var)
			x.delta<-(xlim[2]-xlim[1])/1000
			x.dloc<-sort(c(x.loc,(x.loc[-1]+x.delta),max(x.loc)+1)) 	#data loc (with doubled values)
			pred.matrix<-pred.matrix[rep(1:N.x,each=2),]
			x.tloc<-c(x.loc,max(x.loc)+1)						#tick loc
			x.loc<-x.loc+.5								#axis values loc
			xlim<-range(x.dloc)
		}else{
			x.loc<-pretty(x.var)
			xlim<-range(x.var)
			x.loc<-x.loc[x.loc>xlim[1]&x.loc<xlim[2]]
			x.axis<-format(x.loc)
			x.dloc<-unclass(x.var)
			x.tloc<-x.loc
		}
		
		#distances for Y axis
		X1<-xlim[MM$X[a]]							#inner tick
		X2<-xlim[MM$X[a]]+0.08*PN$X[a]*(xlim[2]-xlim[1])	#outer tick
		X3<-xlim[MM$X[a]]+0.15*PN$X[a]*(xlim[2]-xlim[1]) 	#axis values
		X4<-xlim[MM$X[a]]+0.30*PN$X[a]*(xlim[2]-xlim[1]) 	#axis label
		X5<-xlim[MM$X[a]]+0.40*PN$X[a]*(xlim[2]-xlim[1])  	#catagorical tag
	
		if(is.factor(y.var)){
			y.loc<-unclass(y.var)
			ylim<-range(y.loc)
			y.axis<-levels(y.var)
			y.delta<-(ylim[2]-ylim[1])/1000
			y.dloc<-sort(c(y.loc,(y.loc[-1]+y.delta),max(y.loc)+1))
			pred.matrix<-pred.matrix[,rep(1:N.y,each=2)]
			y.tloc<-c(y.loc,max(y.loc)+1)	
			y.loc<-y.loc+.5
			ylim<-range(y.dloc)
		}else{
			y.loc<-pretty(y.var)
			ylim<-range(y.var)
			y.loc<-y.loc[y.loc>ylim[1]&y.loc<ylim[2]]
			y.axis<-format(y.loc)
			y.dloc<-unclass(y.var)
			y.tloc<-y.loc
		}

		#distances for X axis
		Y1<-ylim[MM$Y[a]]
		Y2<-ylim[MM$Y[a]]+0.08*PN$Y[a]*(ylim[2]-ylim[1])
		Y3<-ylim[MM$Y[a]]+0.15*PN$Y[a]*(ylim[2]-ylim[1])
		Y4<-ylim[MM$Y[a]]+0.30*PN$Y[a]*(ylim[2]-ylim[1])
		Y5<-ylim[MM$Y[a]]+0.40*PN$Y[a]*(ylim[2]-ylim[1])


		ZX1<-xlim[MM$ZX[a]]
		ZX2<-xlim[MM$ZX[a]]+0.04*PN$ZX[a]*(xlim[2]-xlim[1])
		ZX3<-xlim[MM$ZX[a]]+0.06*PN$ZX[a]*(xlim[2]-xlim[1])
		ZX4<-xlim[MM$ZX[a]]+0.21*PN$ZX[a]*(xlim[2]-xlim[1])
		ZY1<-ylim[MM$ZY[a]]
		ZY2<-ylim[MM$ZY[a]]+0.04*PN$ZY[a]*(ylim[2]-ylim[1])
		ZY3<-ylim[MM$ZY[a]]+0.06*PN$ZY[a]*(ylim[2]-ylim[1])
		ZY4<-ylim[MM$ZY[a]]+0.21*PN$ZY[a]*(ylim[2]-ylim[1])


		z.loc<-pretty(z.range,shrink.sml=.5)
		z.loc<-z.loc[-c(1,length(z.loc))]


		### make plot ###

		VT<-persp(x=x.dloc, y=y.dloc, z=pred.matrix, zlim=z.range,      	
		      xlab = xlab, ylab = ylab, zlab = "fitted value",     
			theta=theta, phi=phi, r = sqrt(10), d = 3,
		      ticktype = ticktype, mgp = c(5,1,0), xaxt=xaxt, yaxt=yaxt,axes=FALSE,...)

		par(xpd=NA)

		### x axis labels ###
		#ticks#
		XY1<-trans3d(x.tloc,Y1,z.range[1],VT)
		XY2<-trans3d(x.tloc,Y2,z.range[1],VT)
		for(j in 1:length(x.tloc)){lines(c(XY1$x[j],XY2$x[j]),c(XY1$y[j],XY2$y[j]))}
		#axis values#
		N.xval<-length(x.loc) 	#number of x axis values
		XY<-trans3d(x.loc,Y3,z.range[1],VT)
		text(XY$x,XY$y,x.axis)
		#label#
		srt<-atan((XY$y[N.xval]-XY$y[1])/(XY$x[N.xval]-XY$x[1]))
		srt<-(srt*180)/pi
		XY<-trans3d(mean(x.loc),Y4,z.range[1],VT)
		text(XY$x,XY$y,xlab,adj=.5,srt=srt,font=2)
		if(is.factor(x.var)){
			XY<-trans3d(mean(x.loc),Y5,z.range[1],VT)
			text(XY$x,XY$y,"(categorical)",adj=.5,srt=srt,font=2)}



		### y axis labels ###
		#ticks#
		XY1<-trans3d(X1,y.tloc,z.range[1],VT)
		XY2<-trans3d(X2,y.tloc,z.range[1],VT)
		for(j in 1:length(y.tloc)){lines(c(XY1$x[j],XY2$x[j]),c(XY1$y[j],XY2$y[j]))}
		#axis values#
		N.yval<-length(y.loc) 	#number of y axis values
		XY<-trans3d(X3,y.loc,z.range[1],VT)
		text(XY$x,XY$y,y.axis)
		#label#
		srt<-atan((XY$y[N.yval]-XY$y[1])/(XY$x[N.yval]-XY$x[1]))
		srt<-(srt*180)/pi
		XY<-trans3d(X4,mean(y.loc),z.range[1],VT)
		text(XY$x,XY$y,ylab,adj=.5,srt=srt,font=2)
		if(is.factor(y.var)){
			XY<-trans3d(X5,mean(y.loc),z.range[1],VT)
			text(XY$x,XY$y,"(categorical)",adj=.5,srt=srt,font=2)}

		### z axis labels ###
		#ticks#
		XY1<-trans3d(ZX1,ZY1,z.loc,VT)
		XY2<-trans3d(ZX2,ZY2,z.loc,VT)
		for(j in 1:length(z.loc)){lines(c(XY1$x[j],XY2$x[j]),c(XY1$y[j],XY2$y[j]))}
		#values#
		XY<-trans3d(ZX3,ZY3,z.loc,VT)
		text(XY$x,XY$y,format(z.loc),adj=1)
		#label#
		
		srt<-atan((max(XY$y)-min(XY$y))/(max(XY$x)-min(XY$x)))
		srt<-(srt*180)/pi
		XY<-trans3d(ZX4,ZY4,mean(z.loc),VT)
		text(XY$x,XY$y,zlab,adj=.5,srt=180-srt,font=2)

		par(xpd=TRUE)

	}
}
#################################################################################

if(!device.type[i]%in%c("default","none")){dev.off()}
}
#}
#################################################################################

}
