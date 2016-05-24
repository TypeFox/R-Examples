
#############################################################################################
#############################################################################################
##################################### build.rastLUT #########################################
#############################################################################################
#############################################################################################


build.rastLUT<-function(	imageList=NULL,
					predList=NULL,
					qdata.trainfn=NULL,
					rastLUTfn=NULL,
					folder=NULL
					){

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
	Filters<-rbind(Filters,csv=c("Comma-delimited files (*.csv)", "*.csv"))
}else{
	#stop("build.rastLUT() only works in windows")
}

#############################################################################################
################################ Select Output Folder #######################################
#############################################################################################

if(is.null(folder)){
	if(.Platform$OS.type=="windows"){
		folder<-choose.dir(default=getwd(), caption="Select output folder")
	}else{
		folder<-getwd()}
}

#############################################################################################
###################################### Load Data ############################################
#############################################################################################

if(is.null(predList)){


	## If training data is NULL, then the user selects file from pop-up browser.
	if (is.null(qdata.trainfn)){
		if(.Platform$OS.type=="windows"){
			qdata.trainfn <- choose.files(caption="Select data file", filters = Filters["csv",], multi = FALSE)
			if(is.null(qdata.trainfn)){stop("")}
		}else{
			stop("in nonwindows environment must provide 'predList' or 'qdata.trainfn'")
		}
	}

	## Check if file name is full path or basename
	if(is.matrix(qdata.trainfn)!=TRUE && is.data.frame(qdata.trainfn)!=TRUE){
		if(identical(basename(qdata.trainfn),qdata.trainfn)){qdata.trainfn<-paste(folder,"/",qdata.trainfn,sep="")}
	}

	## Read in training data
	if(is.matrix(qdata.trainfn)==TRUE || is.data.frame(qdata.trainfn)==TRUE){
		qdata<-qdata.trainfn
	}else{
		qdata<-read.table(file=qdata.trainfn,sep=",",header=TRUE,check.names=FALSE,as.is=TRUE)}
	## make predList
	predList<-names(qdata)

}



#############################################################################################
################################ Load Libraries #############################################
#############################################################################################

#r#equirergdal)
#r#equireraster)

#############################################################################################
######################################## Select images ######################################
#############################################################################################

if(is.null(imageList)){
	if(.Platform$OS.type=="windows"){
		imageList <- getRasts()
	}else{
		stop("in nonwindows environment must provide 'imageList'")
	}
}

## find number of bands in each stack ##

numbands <- rep(0,length(imageList))
for(i in 1:length(imageList)){
	rast <- imageList[i]
	if (!is.na(rast)){
		## Gets spatial information for raster
		sp.rast <- brick(rast)
		if(is.null(sp.rast)){
			stop("invalid raster")}
		## Gets number of bands for raster object
		numbands[i] <- nlayers(sp.rast)
	}
}

#############################################################################################
########################################## format LUT #######################################
#############################################################################################

rast<-rep(imageList,numbands)
#print(rast)

pred<-rep("not in data",length(rast))
#print(pred)

band<-NULL
for(i in 1:length(numbands)){band=c(band,seq(numbands[i]))}
#print(band)

rastLUT<-data.frame(rast=rast,pred=pred,band=band,stringsAsFactors=FALSE)

#print("whatever")
#print(rastLUT)
#print(rastLUT$rast[1])
#############################################################################################
################################### Select Predictors #######################################
#############################################################################################

for(i in 1:nrow(rastLUT)){

	#print(rastLUT$rast[i])
	print(paste("Select predictors for ",basename(rastLUT$rast[i]),": band ",rastLUT$band[i],sep=""))

	##filler to adjust title
	title<-paste(basename(rastLUT$rast[i]),": band ",rastLUT$band[i],sep="")
	nchar.title<-nchar(title)
	filler<-paste("not in training data",paste(rep(" ",2*nchar.title),collapse=""))

	## Presents list of possible predictors from raster lookup table to user for selection.
	pred <- select.list(	c(predList,filler),
					title=title, 
					multiple = FALSE)
	if(pred==""){
		stop("predictor must be selected")}
	if(pred==filler){
		pred<-"not in training data"}
	rastLUT$pred[i]<-pred

}

#############################################################################################
######################################## Return rastLUT #####################################
#############################################################################################

##create output filename

if(is.null(rastLUTfn)){
	if(is.matrix(qdata.trainfn)==TRUE || is.data.frame(qdata.trainfn)==TRUE || is.null(qdata.trainfn)){
		rastLUTfn<-paste(folder,"/rastLUT.csv",sep="")
	}else{
		rastLUTfn<-strsplit(basename(qdata.trainfn),".csv")[[1]]
		rastLUTfn<-paste(folder,"/",rastLUTfn,"_rastLUT.csv",sep="")
	}
}

##write table

write.table(	rastLUT,  
			file = rastLUTfn, 
			sep=",",
			append = FALSE, 
			row.names = FALSE,
			col.names = FALSE)

return(rastLUT)
}


