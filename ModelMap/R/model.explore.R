model.explore<-function(qdata.trainfn=NULL,
				folder=NULL,		
				predList=NULL,
				predFactor=FALSE,

				response.name=NULL,
				response.type=NULL,
				response.colors=NULL,	
				unique.rowname=NULL,

				OUTPUTfn=NULL,

			# Graphics Arguments
				device.type=NULL,	
				allow.default.graphics=FALSE,
				res=NULL,
				jpeg.res=72,
				MAXCELL=100000,
				device.width=NULL,
				device.height=NULL,
				units="in",
				pointsize=12,
				cex=1,
				#cex=par()$cex,

			# Raster arguments
				rastLUTfn=NULL,
				create.extrapolation.masks=FALSE,
				na.value=-9999,
				col.ramp=rainbow(101,start=0,end=.5),
				col.cat=palette()[-1]
){

NAval<-na.value
#MAXCELL <- 100000

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

#############################################################################################
################################# Select Response Type ######################################
#############################################################################################

## If model.obj null, ask for response.type

if(is.null(response.type)){
	response.type <- select.list(c("continuous","binary","categorical"), title="Select response type.")}
if(response.type=="" || is.null(response.type)){
	stop("'response.type' needed")}	

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
###################################### Load Data ############################################
#############################################################################################

## If training data is NULL, then the user selects file from pop-up browser.
if (is.null(qdata.trainfn)){
	if(.Platform$OS.type=="windows"){
		qdata.trainfn <- choose.files(caption="Select data file", filters = Filters["csv",], multi = FALSE)
		if(is.null(qdata.trainfn)){stop("")}
	}else{stop("to create a model or make validation predictions you must provide 'qdata.trainfn'")}
}

## Check if file name is full path or basename
qdata.basename<-NULL
if(is.matrix(qdata.trainfn)!=TRUE && is.data.frame(qdata.trainfn)!=TRUE){
	qdata.basename<-basename(qdata.trainfn)
	if(identical(basename(qdata.trainfn),qdata.trainfn)){
		qdata.trainfn<-file.path(folder,qdata.trainfn)}
}

## Read in training data
if(is.matrix(qdata.trainfn)==TRUE || is.data.frame(qdata.trainfn)==TRUE){
	qdata<-qdata.trainfn
}else{
	qdata<-read.table(file=qdata.trainfn,sep=",",header=TRUE,check.names=FALSE,as.is=TRUE)}


#############################################################################################
#################################### Pick Response ##########################################
#############################################################################################

## If the response variable is NULL, then the user selects variable from pop-up list.
if (is.null(response.name)){
	response.name <- select.list(names(qdata), title="Select response name.")
	if(response.name=="" || is.null(response.name)){
		stop("'response.name' is needed")}	
}

#print(paste("response.name =",response.name))

if(!response.name%in%names(qdata)){
	stop("'response.name' ",response.name,"must be a column name in 'qdata.trainfn'")}

if(response.type=="categorical"){
	qdata[,response.name]<-as.factor(qdata[,response.name])

	if(is.null(response.colors)){
		Ncat<-length(levels(qdata[,response.name]))
		if(Ncat==1){
			warning("Response contains only 1 category")}

		#if(Ncat==1){
		#	COLORS<-"grey50";warning("Response contains only 1 category")}
		#if(Ncat==2){
		#	COLORS<-c("gray20","gray80")}
		#if(Ncat>=3 && Ncat<=12){
		#	COLORS<-brewer.pal(Ncat,"Set3")}
		#if(Ncat>12){
		#	COLORS<-rep_len(brewer.pal(12,"Set3"),Ncat)}
		response.colors<-data.frame(	category=levels(qdata[,response.name]),
							colors=1:Ncat)
	}else{
		if(is.factor(response.colors$colors)){response.colors$colors<-as.character(response.colors$colors)}
	}
	
}

if(response.type=="binary"){
	qdata.y<-qdata[,response.name]
	if(!is.numeric(qdata.y)){
		stop("If 'response.type is 'Binary' then 'response.name' must be numeric")}
	qdata.y[qdata.y>0]<-1
	qdata.y[qdata.y<0]<-0
	qdata.y<-as.factor(qdata.y)
	qdata[,response.name]<-qdata.y
}

#############################################################################################
################################### Select Predictors #######################################
#############################################################################################

## This section gets the list of predictors from the user (either explicitly or through
##	pop-up window user selection. If the predictor list is NULL, allows user to  
##	select the predictors from collumn names of training data

#print("Select predictors")

if(is.null(predList)){

	## Presents list of possible predictors from raster lookup table to user for selection.
	predList = select.list(names(qdata), "Select predictors", multiple = TRUE)
	if(length(predList)==0 || is.null(predList)){
		stop("predictors must be selected")}

	
	## asks if any predictors are factors
	any.factors<-select.list(c("NO","YES"), title="any predictors catagorical?")
	if(any.factors=="YES"){
		predFactor<-select.list(predList, "Select catagorical predictors", multiple = TRUE)
		if(length(predFactor)==0 || is.null(predFactor)){
			predFactor<-FALSE
		}
	}	
}

predMissing<-predList[!predList%in%names(qdata)]
if(length(predMissing)>0){
	predMissing.paste<-paste(predMissing,collapse=" ")
	stop("Predictors: ", predMissing.paste," from 'predList' not found in 'qdata.trainingfn'")}

Npred<-length(predList)
#############################################################################################
############################## Select Factored Predictors ###################################
#############################################################################################

#print("Select factored predictors")
#print(paste("     predFactor =",predFactor)

factored.qdata<-sapply(qdata[,match(predList,names(qdata))],is.factor)
character.qdata<-sapply(qdata[,match(predList,names(qdata))],is.character)
factored.qdata<-factored.qdata|character.qdata

if(any(predFactor==FALSE)){
	if(any(factored.qdata)){
		fact.q<-paste(names(factored.qdata)[factored.qdata],collapse=" ")
		stop(	"predictors ",fact.q,
			" are catagorical predictors i.e. are non numeric such as factors or characters but are not included in 'predFactor' either add these predictors to 'predFactor' or correct the dataset"
			)
	}
}

if(!any(predFactor==FALSE)){

	factMissing<-predFactor[!predFactor%in%predList]
	if(length(factMissing)>0){
		factMissing.paste<-paste(factMissing,collapse=" ")
		stop("Factored Predictors: ", factMissing.paste," from 'predFactor' must be included in 'predList'")}

	if(any(!names(factored.qdata)[factored.qdata]%in%predFactor)){
		fact.q<-paste(names(factored.qdata)[factored.qdata][!names(factored.qdata)[factored.qdata]%in%predFactor],collapse=" ")
		stop(	"predictors ",fact.q,
			" are catagorical predictors i.e. are non-numeric such as factors or characters but are not included in 'predFactor' either add these predictors to 'predFactor' or correct the dataset"
			)
	}

	for(i in 1:length(predFactor)){
		qdata[,predFactor[i]]<-factor(qdata[,predFactor[i]])
	}
	
}

#############################################################################################
######################## Select unique row identifier #######################################
#############################################################################################

if (is.null(unique.rowname)){
	unique.rowname <- select.list(c(names(qdata),"row_index"), title="Select unique row identifier")	
	if(unique.rowname!="row_index" && unique.rowname!=""){
		if(anyDuplicated(qdata[,unique.rowname])){
			stop("'unique.rowname' ",unique.rowname," contains duplicated values")}
		rownames(qdata)<-qdata[,unique.rowname]
	}	
}else{
	if(!(unique.rowname%in%names(qdata))){
		warning("'unique.rowname' ",unique.rowname," not found in qdata, row index numbers will be used instead")
		unique.rowname<-FALSE
		rownames(qdata)<-1:nrow(qdata)
	}
	if(unique.rowname!=FALSE){
		if(anyDuplicated(qdata[,unique.rowname])){
			stop("'unique.rowname' ",unique.rowname," contains duplicated values")}
		rownames(qdata)<-qdata[,unique.rowname]
	}
}

########################################################################################
############################## check device type #######################################
########################################################################################

if(allow.default.graphics){
	warning(	"Use with caution.\n",
			"Do not move or close graphics device while 'model.explore()' is running.\n",
			"Closing graphic device before function finishes may crash R session.",
			immediate.=TRUE)
	device.type<-check.device.type(device.type)
}else{
	device.type<-check.device.type.nodefault(device.type)
}

if(is.null(res)){res<-jpeg.res}

#############################################################################################
######################### Drop unused columns of qdata ######################################
#############################################################################################

qdata<-qdata[,c(predList,response.name)]

#############################################################################################
########################### Omit rows with NA values ########################################
#############################################################################################

#qdata[,predList][qdata[,predList,drop=FALSE]==NAval]<-NA
#qdata[,response.name][qdata[,response.name]==NAval]<-NA

NA.pred<-apply(qdata[,predList,drop=FALSE],1,function(x){any(is.na(x))})
NA.resp<-is.na(qdata[,response.name])

if(any(NA.pred)){
	warning("Omiting ", sum(NA.pred), " datapoints with NA predictors")
	qdata <- qdata[!NA.pred,]
}

if(any(NA.resp)){
	warning("Omiting ", sum(NA.resp), " datapoints with NA response values")
	qdata <- qdata[!NA.resp,]
}

#############################################################################################
####################### Generate Raster Output File Names ###################################
#############################################################################################

## OUTPUTfn
if(is.null(OUTPUTfn)){
	OUTPUTfn<- paste("Exploratory_",response.type,"_",response.name,sep="")}

OUTPUTfn <- suppressWarnings(FNcheck(	OUTPUTfn=OUTPUTfn,
							folder=folder,
							ERROR.NAME="OUTPUTfn"))

### After FNcheck, output filename will have path, base and extension in every case

OUTPUTbase  <- basename(OUTPUTfn)				#name and extension, no path

OUTPUTsplit <- strsplit(OUTPUTbase,split="\\.")[[1]]
OUTPUTname <- OUTPUTsplit[1]  				#name, no extension or path
OUTPUText   <- paste(".",OUTPUTsplit[2],sep="")		#just extension

OUTPUTpath  <- dirname(OUTPUTfn)				#just path

OUTPUTfn.noext<-file.path(OUTPUTpath,OUTPUTname)	#path and name, no extension


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
						rastLUTfn=paste(OUTPUTfn.noext,"_rastLUT.csv",sep=""))	
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


########################################################################################
#################################### data range ########################################
########################################################################################

predContinuous<-predList[!predList%in%predFactor]

N.Continuous<-length(predContinuous)

if(any(predFactor==FALSE)){
	N.Factor<-0
}else{N.Factor<-length(predFactor)}

PRED.RANGE<-vector(mode = "list", length = length(predList))
names(PRED.RANGE)<-predList

if(N.Continuous>0){
	for(pred in predContinuous){
		PRED.RANGE[[pred]]<-range(qdata[,pred])}}
if(N.Factor>0){
	for(pred in predFactor){
		PRED.RANGE[[pred]]<-sort(as.numeric(as.character(levels(qdata[,pred]))))}}

########################################################################################
########################### Build raster stacks ########################################
########################################################################################

#require(raster)

RAST<-vector("list", 0)

for(p in predList){
	print(p)
	rastfn<-rastLUT[rastLUT[,2]==p,1]
	band<-  rastLUT[rastLUT[,2]==p,3]

	RAST[[p]]<-raster(rastfn,band=band)
}

RAST<-projfix(RAST,OUTPUTfn.noext=OUTPUTfn.noext)
RS<-stack(RAST)

#RB<-brick(RS,values=TRUE,filename=OUTPUTfn.brick,overwrite=TRUE)

########################################################################################
#################################### raster range ######################################
########################################################################################

valFun<-function(stratlayer,NAval){
	bs <- blockSize(stratlayer)
	vals <- {}
	for(i in 1:bs$n){
	  vals <- unique(c(vals, getValuesBlock(stratlayer, row=bs$row[i], nrows=bs$nrows[i])))
	}
	vals <- vals[!is.na(vals)]
	vals <- vals[vals!=NAval]
	return(sort(vals))
}

#rangeFun<-function(stratlayer,NAval){
#	bs <- blockSize(stratlayer)
#	ranges <- {}
#	for(i in 1:bs$n){
#        bval<-getValuesBlock(stratlayer, row=bs$row[i], nrows=bs$nrows[i])
#        bval<-bval[bval!=NAval]
#	  if(length(bval)>0){ranges <- range(c(ranges,bval))}
#	}
#	return(ranges)
#}


RAST.RANGE<-vector(mode = "list", length = length(predList))
names(RAST.RANGE)<-predList

#if(N.Continuous>0){
#	for(pred in predContinuous){
#		RAST.RANGE[[pred]]<-rangeFun(RAST[[pred]],NAval)
#		print(paste(pred," range: ",paste(RAST.RANGE[[pred]],collapse=" "),sep=""))}}
if(N.Factor>0){
	for(pred in predFactor){
		if(ncell(RAST[[pred]])>MAXCELL){
			RAST.RANGE[[pred]]<-valFun(RAST[[pred]],NAval)
		}else{
			RAST.RANGE[[pred]]<-unique(RAST[[pred]])
		}
		print(paste(pred," range: ",paste(RAST.RANGE[[pred]],collapse=" "),sep=""))
	}
}
print("raster ranges calculated")

########################################################################################
#################################### sample raster #####################################
########################################################################################

SAMP <- FALSE

if(ncell(RS)>MAXCELL){
	SAMP<-TRUE
	print("starting sampling")
	RB.samp<-sampleRegular(RS, size=MAXCELL, asRaster=TRUE)
	print("sampling done")
}else{
	RB.samp<-brick(RS)
}

########################################################################################
########################## NA values in sampled rasters ################################
########################################################################################

	#RB brick of all predictor layers, NA for all NAval in each layer
	#RL single layer, NA if any layer of RB is NA, else 1
	#
	#MB brick, mask for each layer, 1 if inside range, NA if RB is NA, or RB is outside range
	#ML single layer, TRUE if all layers of MB have data, FALSE if any layer is NA


	print("extrapolation mask from sampled rasters")
	NA.fun <- function(x) { x[x==NAval] <- NA; return(x) }
	RB.samp <- calc(RB.samp, NA.fun, overwrite=TRUE)
	if(Npred==1){RB.samp<-brick(RB.samp)}
	names(RB.samp)<-names(RS)

	RL.samp<-if(Npred>1){sum(RB.samp)}else{RB.samp}
	RL.samp[!is.na(RL.samp)]<-1

	MB.samp <- RB.samp
	if(N.Continuous>0){for(pred in predContinuous){
		#if(Npred>1){
			MB.samp[[pred]][MB.samp[[pred]]<PRED.RANGE[[pred]][1]]<-NA
			MB.samp[[pred]][MB.samp[[pred]]>PRED.RANGE[[pred]][2]]<-NA
		#}
		#if(Npred==1){
		#	MB.samp[MB.samp<PRED.RANGE[[pred]][1]]<-NA
		#	MB.samp[MB.samp>PRED.RANGE[[pred]][2]]<-NA
		#}			
	}}
	if(N.Factor>0){for(pred in predFactor){
		#if(Npred>1){
			MB.samp[[pred]][!MB.samp[[pred]]%in%PRED.RANGE[[pred]]]<-NA
		#}
		#if(Npred==1){
		#	MB.samp[!MB.samp%in%PRED.RANGE[[pred]]]<-NA
		#}
	}}
	MB.samp[!is.na(MB.samp)]<-1

	#ML.samp<-if(Npred>1){!any(is.na(MB.samp))}else{!is.na(MB.samp)} 
	ML.samp<-if(Npred>1){sum(MB.samp)}else{MB.samp}
	ML.samp[!is.na(ML.samp)]<-1

	print("finished extrapolation mask in sampled rasters")

########################################################################################
######################### check device height and width ################################
########################################################################################

if(is.null(device.height)){
	if(!is.null(device.width)){
		device.height<-device.width
	}else{
		device.height<-7
	}
}

########################################################################################
################################### CORR plot ##########################################
########################################################################################
if(N.Continuous>1){
	correlation.function(	qdata=qdata,
					predList=predList,
					predFactor=predFactor,
					MODELpredfn=OUTPUTfn.noext,

					device.type=device.type,
					res=res,
					device.width=device.height,
					device.height=device.height,
					units=units,
					pointsize=pointsize,
					cex=cex	
					)
}

########################################################################################
################################### pred plot ##########################################
########################################################################################

if(N.Continuous>0){
	#pred<-predContinuous[1]
	for(pred in predContinuous){
		explore.continuous(	qdata=qdata,
						response.name=response.name,
						response.type=response.type,
						response.colors=response.colors,
						pred=pred,
						pred.range=PRED.RANGE[[pred]],
						OUTfn=paste(OUTPUTfn.noext,"_",pred,sep=""),
						main=OUTPUTfn.noext,

						pred.rast=if(Npred==1){RB.samp[[1]]}else{RB.samp[[pred]]},
						pred.mask=if(Npred==1){MB.samp[[1]]}else{MB.samp[[pred]]},
						#MAXCELL=MAXCELL,
					
						col.ramp=col.ramp,
	
						device.type=device.type,
						res=res,
						device.width=if(is.null(device.width)){1.8*device.height}else{device.width},
						device.height=device.height,
						units=units,
						pointsize=pointsize,
						cex=cex
						)
	}
}

if(N.Factor>0){
	for(pred in predFactor){
		print("starting explore.categorical()")

		explore.categorical(	qdata=qdata,
						response.name=response.name,
						response.type=response.type,
						response.colors=response.colors,
						pred=pred,
						pred.range=PRED.RANGE[[pred]],
						rast.range=RAST.RANGE[[pred]],
						OUTfn=paste(OUTPUTfn.noext,"_",pred,sep=""),
						main=OUTPUTfn.noext,

						pred.rast=if(Npred==1){RB.samp[[1]]}else{RB.samp[[pred]]},
						pred.mask=if(Npred==1){MB.samp[[1]]}else{MB.samp[[pred]]},
						#MAXCELL=MAXCELL,
					
						col.cat=col.cat,

						device.type=device.type,
						res=res,
						device.width=if(is.null(device.width)){1.8*device.height}else{device.width},
						device.height=device.height,
						units=units,
						pointsize=pointsize,
						cex=cex
						)
	print("done explore.categorical()")
	}
}


########################################################################################
############################# overall mask graphic #####################################
########################################################################################

print("starting mask.graphic") 
mask.graphic(	RL=RL.samp,
			ML=ML.samp,
			OUTfn=paste(OUTPUTfn.noext,"_mask_all_predictors",sep=""),
			main="Extrapolation Mask - All Predictors",
			device.type=device.type,
			res=res,
			device.width=if(is.null(device.width)){0.9*device.height}else{0.5*device.width},
			device.height=device.height,
			units=units,
			pointsize=pointsize,
			cex=cex
)

print("done mask.graphic") 

########################################################################################
########################### create extrapolation masks #################################
########################################################################################

if(create.extrapolation.masks){
   if(SAMP){
   ###SAMP==TRUE###

		### Creat filename for full predictor raster bricks ###

		#OUTPUTfn.brick  <- rasterTmpFile(prefix=paste("raster_tmp_",OUTPUTname,"_brick_",sep=""))
		#OUTPUTfn.raster <- rasterTmpFile(prefix=paste("raster_tmp_",OUTPUTname,"_raster_",sep=""))
		MASKfn.brick    <- rasterTmpFile(prefix=paste("raster_tmp_",OUTPUTname,"_mask_",sep=""))
		MASKfn.raster   <- rasterTmpFile(prefix=paste("raster_tmp_",OUTPUTname,"_maskraster_",sep=""))

		#RB <- brick(RS[[1]],nl=nlayers(RS),filename=OUTPUTfn.brick,values=FALSE)
		#RL <- raster(RS[[1]])
		MB <- brick(RS[[1]],nl=nlayers(RS),values=FALSE, filename=MASKfn.brick, overwrite=TRUE)
		ML <- raster(RAST[[1]])

		#RB <- writeStart(RB, filename=OUTPUTfn.brick, overwrite=TRUE)
		#RL <- writeStart(RL, filename=OUTPUTfn.raster, overwrite=TRUE)
		MB <- writeStart(MB, filename=trim(MASKfn.brick),  overwrite=TRUE)
		ML <- writeStart(ML, filename=trim(MASKfn.raster), overwrite=TRUE)

		bs <- blockSize(MB)
		for (i in 1:bs$n) {
			v <- getValues(RS, row=bs$row[i], nrows=bs$nrows[i] )
			if(nlayers(RS)==1){colnames(v)<-predList} #maybe not needed?
	
			#v=RB#
			v[v==NAval] <- NA
	
			#v2=RL#	
			#v2<-if(Npred>1){apply(v,1,sum)}else{v}
			#v2[!is.na(v2)]<-1

			#v3=MB#
			v3<-v
			if(N.Continuous>0){for(pred in predContinuous){
				v3[,pred][v3[,pred]<PRED.RANGE[[pred]][1]]<-NA
				v3[,pred][v3[,pred]>PRED.RANGE[[pred]][2]]<-NA
			}}
			if(N.Factor>0){for(pred in predFactor){
				v3[,pred][!v3[,pred]%in%PRED.RANGE[[pred]]]<-NA
			}}
			v3[!is.na(v3)]<-1
	
			#v4=ML#
			#v4<-if(Npred>1){apply(v3,1,function(x){!any(is.na(x))})}else{!is.na(v3)}
			v4 <-if(Npred>1){apply(v3,1,sum)}else{v3}
			v4[!is.na(v4)]<-1


			#RB <- writeValues(RB, v,  bs$row[i])
			#RL <- writeValues(RL, v2, bs$row[i])
			MB <- writeValues(MB, v3, bs$row[i])	
			ML <- writeValues(ML, v4, bs$row[i])
		}
		#RB <- writeStop(RB)
		#RL <- writeStop(RL)
		MB <- writeStop(MB)
		ML <- writeStop(ML)
	
		#names(RB)<-names(RS)
		names(MB)<-names(RS)

		writeRaster(MB, filename=paste(OUTPUTfn.noext,"_mask.img",sep=""),overwrite=TRUE,datatype="INT1U",NAflag=0,suffix="names",bylayer=TRUE)
		writeRaster(ML, filename=paste(OUTPUTfn.noext,"_mask_all_predictors.img",sep=""),overwrite=TRUE,datatype="INT1U",NAflag=0)

		if(filename(MB)!=""){
			FILENAME<-filename(MB)
			rm(MB)
			file.remove(FILENAME)
			file.remove(grd2gri(FILENAME))
		}

		if(filename(ML)!=""){
			FILENAME<-filename(ML)
			rm(ML)
			file.remove(FILENAME)
			file.remove(grd2gri(FILENAME))
		}

	}else{ 
	###SAMP==FALSE###
		writeRaster(MB.samp, filename=paste(OUTPUTfn.noext,"_mask.img",sep=""),overwrite=TRUE,datatype="INT1U",NAflag=0,suffix="names",bylayer=TRUE)
		writeRaster(ML.samp, filename=paste(OUTPUTfn.noext,"_mask_all_predictors.img",sep=""),overwrite=TRUE,datatype="INT1U",NAflag=0)
	}
}


########################################################################################
}

