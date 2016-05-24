
#############################################################################################
#############################################################################################
################################# Model Map Making ##########################################
#############################################################################################
#############################################################################################



model.mapmake<-function(model.obj=NULL,
				folder=NULL,		# No ending slash, to output to working dir = getwd()
				MODELfn=NULL,
				rastLUTfn=NULL,
			# Model Evaluation Arguments
				na.action=NULL,	# also used for mapping
				na.value=-9999,	# NA value for rasters, note, internal code changes this to NAval
			# Mapping arguments
				keep.predictor.brick=FALSE,	
				map.sd=FALSE,
				OUTPUTfn=NULL,
			# QRF arguments
				quantiles=NULL,
			# SGB arguments
				n.trees=NULL
){


## Note: Must have R version 2.5.1 or greater
## Note: Must have all rasters in same size, same resolution, and with -9999 for nodata values.
## Note: The names of the predictors must be the same as the predictor names in the training/test data.
##		If the predictors are part of a layer stack, the names must include the number of the
##		layer within the stack (Example: stacknameb1, stacknameb2, stacknameb3, ...)
## Note: If subsampled test data is desired, run training selection.R first.

## Note: Must have the following packages installed: rgdal, gbm, randomForest

# Note:	'rgdal' is only used to make map
#	'gbm'	is only used for SGB models
#	'randomForest' is used for RF models, and also for na.action="na.roughfix" for all model types.


#############################################################################################
################################### Check Platform ##########################################
#############################################################################################


Rplatform<-.Platform$OS.type


#############################################################################################
################################### Add to Filters Table ####################################
#############################################################################################


## Adds to file filters to Cran R Filters table.
if(.Platform$OS.type=="windows"){
	Filters<-rbind(Filters,img=c("Imagine files (*.img)", "*.img"))
	Filters<-rbind(Filters,csv=c("Comma-delimited files (*.csv)", "*.csv"))}

#####################################################################################
############################# check model.obj #######################################
#####################################################################################

if(is.null(model.obj)){
	if(is.null(MODELfn)){
		if(.Platform$OS.type=="windows"){
			MODELfn <- choose.files(caption="Select model object", filters = Filters["All",], multi = FALSE)
			if(is.null(MODELfn)){stop("must provide a model object")}
		}else{stop("must provide a model object")}
	}
	modelname<-load(MODELfn)
	if(length(modelname)!=1){
		stop("file must contain single model object")}
	assign("model.obj",get(modelname))
}



#####################################################################################
##################### Extract Model Type from model.obj #############################
#####################################################################################

model.type<-check.model.type(model.obj)

#############################################################################################
########################### Load modelling packages #########################################
#############################################################################################

if(model.type=="CF"){REQUIRE.party()}
if(model.type=="QRF"){REQUIRE.quantregForest()}
if(model.type=="SGB" || model.type=="QSGB"){REQUIRE.gbm()}


#####################################################################################
############################# Check if map.sd legal #################################
#####################################################################################

if(map.sd){
	if(!model.type%in%c("RF")){
		if(!(model.type=="QRF" && "RF"%in%names(model.obj))){
			warning("Standard deviation maps only available for \"RF\" models")
			map.sd<-FALSE
		}
	}
	if(response.type!="continuous"){
		warning("Standard deviation maps only available for continuous response models")
		map.sd<-FALSE}
}

#####################################################################################
############################# Check Quantiles #######################################
#####################################################################################

if(model.type=="QRF"){
	if(is.null(quantiles)){
		model.quantiles<-if("QRF"%in%names(model.obj)){model.obj$QRF$quantiles}else{model.obj$quantiles}
		quantiles<-if(!is.null(model.quantiles)){model.quantiles}else{c(0.1, 0.5, 0.9)}
	}
	if(!is.numeric(quantiles)){
		stop("'quantiles' must be numbers between 0 and 1")}
	if(any(quantiles<=0) || any(quantiles>=1)){
		stop("'quantiles' must be numbers between 0 and 1")}
}


#############################################################################################
########################### Extract predList from model.obj #################################
#############################################################################################

predList<-check.predList(model.obj=model.obj,model.type=model.type)

#############################################################################################
####################### Extract Factored Predictors from model.obj ##########################
#############################################################################################

model.levels<-check.model.levels(model.obj=model.obj,model.type=model.type)

#predFactor<-NULL
#if(!is.null(model.levels)){predFactor<-names(model.levels)}

#############################################################################################
###################### Extract response.name from model.obj #################################
#############################################################################################

response.name<-check.response.name(model.obj=model.obj,model.type=model.type,response.name=NULL,qdata=NA)

#############################################################################################
###################### Extract response.type from model.obj #################################
#############################################################################################

response.type<-check.response.type(model.obj=model.obj,model.type=model.type,ONEorTWO="model.obj")

#############################################################################################
################################ Load Libraries #############################################
#############################################################################################
## Loads necessary libraries.

##require(rgdal)
#require(raster)

#if(model.type=="RF" ){
#	require(randomForest)
#}
#if(model.type=="SGB"){
#	#require(gbm)
#}
##if(response.type=="categorical"){library(gdata)}

# Note:	'rgdal' is only used to make map
#	'gbm'	is only used for SGB models
#	'randomForest' is used for RF models, and also for na.action="na.roughfix" for all model types.
#	'gdata' is (no longer) used for mapping levels of categorical reponse variables



#############################################################################################
################################ Select Output Folder #######################################
#############################################################################################

if(is.null(folder)){
	if(.Platform$OS.type=="windows"){
		folder<-choose.dir(default=getwd(), caption="Select directory")
	}else{
		folder<-getwd()}
}

#############################################################################################
############################# Generate Output File Names ####################################
#############################################################################################

if(is.null(OUTPUTfn)){
	if(is.null(MODELfn)){
		OUTPUTfn<-paste(model.type,"_",response.type,"_",response.name,"_map",sep="")
	}else{
		OUTPUTfn<-paste(extension(MODELfn,""),"_map",sep="")
	}
}

OUTPUTfn <- FNcheck(	OUTPUTfn=OUTPUTfn,
				folder=folder,
				ERROR.NAME="OUTPUTfn")

### After FNcheck, output filename will have path, base and extension in every case

OUTPUTbase  <- basename(OUTPUTfn)				#name and extension, no path

OUTPUTsplit <- strsplit(OUTPUTbase,split="\\.")[[1]]
OUTPUTname <- OUTPUTsplit[1]  				#name, no extension or path
OUTPUText   <- paste(".",OUTPUTsplit[2],sep="")		#just extension

OUTPUTpath  <- dirname(OUTPUTfn)				#just path

OUTPUTfn.noext<-file.path(OUTPUTpath,OUTPUTname)	#path and name, no extension



#############################################################################################
############################### Determine NA.ACTION #########################################
#############################################################################################

#NA.ACTION is character string. Choices: "omit", "rough.fix".
#na.action is the actual function, no quotes. Function argument will staart with quotes then have them removed.
#model.NA.ACTION and model.na.action are similar but extracted from model object,
#   and need to be checked against 'model.diagnostics(na.action)'.

NAvalid<-c("na.omit")
model.NAvalid<-c("omit")

NA.ACTION<-NULL

if(!is.null(na.action)){
###if 'na.action' specified in 'model.mapmake()' check if valid, and if so use it.###
	if(is.function(na.action)){
		 stop("ModelMap requires the use of quotes when specifying the argument 'na.action'")
	}else{
		if(!na.action%in%NAvalid){
			stop("ModelMap currently supports only \"na.omit\" for map making 'na.action'")
		}else{
			#make NA.ACTION characters without "na."
			NA.ACTION<-switch(na.action,na.omit="omit",na.roughfix="roughfix","invalid")
			#turn na.action into function
			#na.action<-switch(na.action,na.omit=na.omit,na.roughfix=na.roughfix)
		}
	}
}else{
###if 'na.action not specified in 'model.map make()' then check if model was built with a valid na.action###
	if(model.type!="CF"){
		model.NA.ACTION<-if("QRF"%in%names(model.obj)){class(model.obj$QRF$na.action)}else{class(model.obj$na.action)}
		print(paste("model.NA.ACTION:",model.NA.ACTION))
		if(!is.null(model.NA.ACTION)){
			if(model.NA.ACTION!="NULL"){	
				if(!model.NA.ACTION%in%model.NAvalid){
					warning("Model Object built with 'na.action' ",model.NA.ACTION," not supported for map making therefore default 'na.action' is \"na.omit\"")
				}else{
					NA.ACTION<-model.NA.ACTION	
					#print("Using 'na.action' from 'model.obj'")	
					#na.action<-model.na.action
				}
			}
		}
	}
}
if(is.null(NA.ACTION)){
###if neither function call nor model object gives valid 'na.action' then defaults to "na.omit"###
	NA.ACTION<-"omit"	
	print("Using default 'na.action' of \"na.omit\"")
	#na.action<-na.omit
}

#############################################################################################
############################# SGB + CV: check for n.trees ###################################
#############################################################################################

if(model.type=="SGB" && is.null(n.trees)){
	if(!is.null(model.obj$best.iter)){
		n.trees<-model.obj$best.iter
	}else{
		n.trees<-model.obj$n.trees
	}
}


#############################################################################################
################################## Ask for rastLUT ##########################################
#############################################################################################

### If rastLUTfn is NULL, then the user selects file from pop-up browser.

if(is.null(rastLUTfn)){
	if( .Platform$OS.type=="windows"){
		LUT.available<-select.list(c("YES","NO"), title="rastLUT available?")
		if(LUT.available=="YES"){
			rastLUTfn<-choose.files(caption="Select raster look up table", filters = Filters["csv",], multi = FALSE)
		}else{
			build.rastLUT(	predList=predList,
						rastLUTfn=paste(MODELfn,"_rastLUT.csv",sep=""))	
		}
	}else{
		stop("you must provide a raster Look Up Table")
	}	
}

### Check if rastLUTfn file name is full path or basename

if(is.matrix(rastLUTfn)!=TRUE && is.data.frame(rastLUTfn)!=TRUE){
	if(identical(basename(rastLUTfn),rastLUTfn))
		{rastLUTfn<-file.path(folder,rastLUTfn)
	}
}

### if rastLUT is filename, read in lookup table

if(is.matrix(rastLUTfn)==TRUE || is.data.frame(rastLUTfn)==TRUE){
	rastLUT<-rastLUTfn
	rastLUT<-data.frame(rastLUT)
}else{
	rastLUT<-read.table(file=rastLUTfn,sep=",",header=FALSE,check.names=FALSE,stringsAsFactors=FALSE)
}


### Check that columns of rastLUT are correct format

if(is.factor(rastLUT[,1])){rastLUT[,1]<-as.character(rastLUT[,1])}
if(is.factor(rastLUT[,2])){rastLUT[,2]<-as.character(rastLUT[,2])}
if(is.factor(rastLUT[,3])){rastLUT[,3]<-as.numeric(as.character(rastLUT[,1]))}

if(is.list(rastLUT[,1])){stop("'rastLUT' is of incorrect format, the first collumn of your 'rastLUT' is a list")}
if(is.list(rastLUT[,2])){stop("'rastLUT' is of incorrect format, the second collumn of your 'rastLUT' is a list")}
if(is.list(rastLUT[,3])){stop("'rastLUT' is of incorrect format, the third collumn of your 'rastLUT' is a list")}


### Check that all predictors in predList are in rastLUT

pred.not.in.LUT<-!(predList%in%rastLUT[,2])
if(any(pred.not.in.LUT)){
	predNot<-paste(predList[pred.not.in.LUT]," ",sep="")
	stop("Predictors ",predNot,"from predList are not found in rastLUT")}


#############################################################################################
######################################### Make Map ##########################################
#############################################################################################

## Begin prediction

print("starting production prediction")


production.prediction(	model.obj=model.obj,
				model.type=model.type,
				rastLUT=rastLUT,
				#na.action=na.action,
				NA.ACTION=NA.ACTION,
				NAval=na.value,
				model.levels=model.levels,
				response.type=response.type,
				#response.name=response.name,
				keep.predictor.brick=keep.predictor.brick,
				map.sd=map.sd,
				OUTPUTfn=OUTPUTfn,			#path, name, extension
				OUTPUTfn.noext=OUTPUTfn.noext,	#path, name
				#OUTPUTpath=OUTPUTpath,			#path
				OUTPUTname=OUTPUTname,			#name
				OUTPUText=OUTPUText,			#extension
				quantiles=quantiles,
				n.trees=n.trees)

#print("finished production prediction")
#############################################################################################
############## If response/type=="categorical Write a key to levels #########################
#############################################################################################

if(response.type=="categorical"){
	#mapkey<-mapLevels(model.obj$y)

	#Ylev<-levels(model.obj$y)
	#Ylev<-model.obj$classes
	if(model.type=="RF"){ Ylev<-colnames(model.obj$votes)}
	#if(model.type=="QRF"){}#not needed as QRF always continuous
	if(model.type=="CF"){Ylev<-model.obj@responses@levels[[response.name]]}
	if(model.type=="SGB"){Ylev<-model.obj$classes}

	if(any(is.na(suppressWarnings(as.numeric(Ylev))))){
		mapkey<-data.frame(row=1:length(Ylev), category=Ylev,integercode=1:length(Ylev))
	}else{
		mapkey<-data.frame(row=1:length(Ylev), category=Ylev,integercode=as.numeric(Ylev))
	}

	MapKeyfn<-paste(OUTPUTfn.noext,"_key.csv",sep="")

	write.table(mapkey,file=MapKeyfn,sep=",",row.names=FALSE)
}
	
#############################################################################################
################################## Write a list of argumets #################################
#############################################################################################


A<-formals(model.mapmake)
A<-mget(names(A),ifnotfound="NULL",envir=as.environment(-1))

ARGfn<-paste(OUTPUTfn.noext,"_mapmake_arguments.txt",sep="")

if(is.matrix(rastLUTfn)==TRUE || is.data.frame(rastLUTfn)==TRUE){
	A$rastLUTfn<-"preloaded dataframe"
}

A$datestamp<-Sys.time()
A<-A[c(length(A),1:(length(A)-1))]

print(paste("ARGfn:",ARGfn))

capture.output(print(A),file=ARGfn)

#############################################################################################
####################################### Return nothing ######################################
#############################################################################################

}
