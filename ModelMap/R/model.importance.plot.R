
########## check model type #################


########################################################################
########################################################################
###########################  importance plot   #########################
########################################################################
########################################################################

model.importance.plot<-function(	model.obj.1=NULL, 
						model.obj.2=NULL, 
						model.name.1="Model 1", 
						model.name.2="Model 2",
						imp.type.1=NULL,
						imp.type.2=NULL,
						type.label=TRUE,
						class.1=NULL,
						class.2=NULL,
						quantile.1=NULL,
						quantile.2=NULL,
						col.1="grey",
						col.2="black", 
						scale.by="sum",
						sort.by="model.obj.1", 

						cf.mincriterion.1 = 0, 
						cf.conditional.1 = FALSE, 
       					cf.threshold.1 = 0.2, 
						cf.nperm.1 = 1,
						cf.mincriterion.2 = 0, 
						cf.conditional.2 = FALSE, 
       					cf.threshold.2 = 0.2, 
						cf.nperm.2 = 1,

						predList=NULL,
						folder=NULL,
						PLOTfn=NULL,
						
					# Graphics Arguments
						device.type=NULL,	
						res=NULL,
						jpeg.res=72,
						device.width=7,
						device.height=7,
						units="in",
						pointsize=12,
						cex=par()$cex,
						...){


#############################################################################################
################################### Add to Filters Table ####################################
#############################################################################################

## Adds to file filters to Cran R Filters table.
if(.Platform$OS.type=="windows"){
	Filters<-rbind(Filters,img=c("Imagine files (*.img)", "*.img"))
	Filters<-rbind(Filters,csv=c("Comma-delimited files (*.csv)", "*.csv"))}

###################################################################################
########################## Check Device Type ######################################
###################################################################################

device.type<-check.device.type(device.type)

if(is.null(res)){res<-jpeg.res}

###################################################################################
######################## Select Output Folder #####################################
###################################################################################

if(is.null(folder)){
	if(any(device.type%in%c("jpeg","pdf","postscript"))){
		if(.Platform$OS.type=="windows"){
			folder<-choose.dir(default=getwd(), caption="Select directory")
		}else{
			folder<-getwd()}
	}
}

###################################################################################
######################## check output filename ####################################
###################################################################################

if(is.null(PLOTfn)){PLOTfn<- paste(model.name.1,"_",model.name.2,sep="")}
if(identical(basename(PLOTfn),PLOTfn)){
	PLOTfn<-file.path(folder,PLOTfn)}

#####################################################################################
############################# check model.obj #######################################
#####################################################################################

if(is.null(model.obj.1)){
	if(is.null(MODELfn)){
		if(.Platform$OS.type=="windows"){
			MODELfn <- choose.files(caption="Select first model", filters = Filters["All",], multi = FALSE)
			if(is.null(MODELfn)){stop("must provide a model object")}
		}else{stop("must provide a model object")}
	}
	modelname<-load(MODELfn)
	if(length(modelname)!= 1){
		stop("file must contain single model object")}
	assign("model.obj.1",get(modelname))
}

if(is.null(model.obj.2)){
	if(is.null(MODELfn)){
		if(.Platform$OS.type=="windows"){
			MODELfn <- choose.files(caption="Select second model", filters = Filters["All",], multi = FALSE)
			if(is.null(MODELfn)){stop("must provide a model object")}
		}else{stop("must provide a model object")}
	}
	modelname<-load(MODELfn)
	if(length(modelname)!= 1){
		stop("file must contain single model object")}
	assign("model.obj.2",get(modelname))
}
#####################################################################################
####################### Check Model Type from model.obj #############################
#####################################################################################

model.type.1<-check.model.type(model.obj.1)
model.type.2<-check.model.type(model.obj.2)

if(model.type.1=="QRF"){if("QRF"%in%names(model.obj.1)){model.obj.1<-model.obj.1$QRF}}
if(model.type.2=="QRF"){if("QRF"%in%names(model.obj.2)){model.obj.2<-model.obj.2$QRF}}

#############################################################################################
########################### Load modelling packages #########################################
#############################################################################################

if(model.type.1=="CF"  || model.type.2=="CF"){REQUIRE.party()}
if(model.type.1=="QRF" || model.type.2=="QRF"){REQUIRE.quantregForest()}
if(model.type.1=="SGB" || model.type.1=="QSGB" || model.type.2=="SGB" || model.type.2=="QSGB"){REQUIRE.gbm()}

#####################################################################################
####################### should warnings be immeadiate ###############################
#####################################################################################
if(model.type.1=="CF" || model.type.2=="CF"){
	WARN.IM<-TRUE
}else{
	WARN.IM<-FALSE
}

#####################################################################################
######################### Check predList from model.obj #############################
#####################################################################################

predList.1<-check.predList(model.obj=model.obj.1,model.type=model.type.1)
predList.2<-check.predList(model.obj=model.obj.2,model.type=model.type.2)


if(!identical(sort(predList.1),sort(predList.2))){
	stop("'model.obj.1' contains different predictors than 'model.obj.2'")
}


if(!is.null(predList)){
	if(!identical(sort(predList),sort(predList.1))){
		stop("predList contains different predictors than models")
	}
}

#############################################################################################
######################## Extract  response.type from model.obj ##############################
#############################################################################################

#Note: response type is only used to print definition of importance used for the two models

response.type.1<-check.response.type(model.obj=model.obj.1,model.type=model.type.1,ONEorTWO="model.obj.1")
response.type.2<-check.response.type(model.obj=model.obj.2,model.type=model.type.2,ONEorTWO="model.obj.2")

#####################################################################################
###################### Extract Importace from model.obj #############################
#####################################################################################

###model 1###

if(model.type.1=="RF"){
	if(is.null(imp.type.1)){imp.type.1<-1}
	if(!is.null(class.1) && imp.type.1==2){
		warning("no class specific measure for 'imp.type.1=2' therefore importance type changed to 'imp.type.1=1'",immediate.=WARN.IM)
		imp.type.1<-1}
	IMP.1<-imp.extract.rf(model.obj.1,imp.type=imp.type.1,class=class.1)
}
if(model.type.1=="QRF"){
	if(is.null(imp.type.1)){imp.type.1<-1}
	if(!is.null(class.1)){
		warning("no class specific measure for QRF therefore 'class.1' ignored",immediate.=WARN.IM)
		class.1<-NULL}
	if(!"quantiles"%in%names(model.obj.1)){stop("QRF 'model.obj.1' was built with 'importance=FALSE' therefore no importances available")}
	quantile.names.1<-paste("quantile=",model.obj.1$quantiles)
	if(is.numeric(quantile.1)){quantile.1<-paste("quantile=",quantile.1)}
	if(!quantile.1%in%quantile.names.1){
		warning("quantile.1 (", quantile.1, ") not found in model.obj.1",immediate.=TRUE)
		quantile.1<-NULL}
	if(is.null(quantile.1)){quantile.1 <- select.list(quantile.names.1, title="Model 1", multiple = FALSE)}
	if(length(quantile.1)==0 || is.null(quantile.1) || !quantile.1%in%quantile.names.1){
		stop("must provide valid quantile for model.obj.1")}
	IMP.1<-imp.extract.qrf(model.obj.1,imp.type=imp.type.1,ONEorTWO="model.obj.1")
	IMP.1<-IMP.1[,c("pred",quantile.1)]
	names(IMP.1)<-c("pred","imp")
	IMP.1<-IMP.1[order(IMP.1$imp),]
}

if(model.type.1=="CF"){
	if(is.null(imp.type.1)){imp.type.1<-1}
	if(!is.null(class.1)){
		warning("no class specific measure for CF models therefore 'class.1' ignored",immediate.=WARN.IM)
		class.1<-NULL}
	if(imp.type.1==2 && response.type.1!="binary"){
		warning("AUC-based variables importances only available for binary response models therefor importance type changed to 'imp.type.1=1'",immediate.=WARN.IM)
		imp.type.1<-1}
	IMP.1<-imp.extract.cf(	model.obj.1,imp.type=imp.type.1,mincriterion=cf.mincriterion.1, 
					conditional=cf.conditional.1,threshold=cf.threshold.1,nperm=cf.nperm.1)
}
if(model.type.1=="SGB"){
	if(is.null(imp.type.1)){imp.type.1<-2}
	if(response.type.1=="categorical" && imp.type.1==1){
		imp.type.1<-2
		warning("'response.type' categorical is not yet supported for permutation importance, relative influence will be used instead",immediate.=WARN.IM)}
	if(!is.null(class.1)){
		warning("no class specific measure for SGB models therefore 'class.1' ignored",immediate.=WARN.IM)}
	print(paste("imp.type.1:",imp.type.1))
	IMP.1<-imp.extract.sgb(model.obj.1,imp.type=imp.type.1)
}

#print("IMP.1")
#print(IMP.1)

###model 2###

if(model.type.2=="RF"){
	if(is.null(imp.type.2)){imp.type.2<-1}
	if(!is.null(class.2) && imp.type.2==2){
		warning("no class specific measure for 'imp.type.2=2' therefore importance type changed to 'imp.type.2=1'",immediate.=WARN.IM)
		imp.type.2<-1}
	IMP.2<-imp.extract.rf(model.obj.2,imp.type=imp.type.2,class=class.2)
}
if(model.type.2=="QRF"){
	if(is.null(imp.type.2)){imp.type.2<-1}
	if(!is.null(class.2)){
		warning("no class specific measure for QRF therefore 'class.2' ignored",immediate.=WARN.IM)
		class.2<-NULL}
	if(!"quantiles"%in%names(model.obj.2)){stop("QRF 'model.obj.2' was built with 'importance=FALSE' therefore no importances available")}
	quantile.names.2<-paste("quantile=",model.obj.2$quantiles)
	if(is.numeric(quantile.2)){quantile.2<-paste("quantile=",quantile.2)}
	if(!quantile.2%in%quantile.names.2){
		warning("quantile.2 (", quantile.2, ") not found in model.obj.2",immediate.=TRUE)
		quantile.2<-NULL}
	if(is.null(quantile.2)){quantile.2 <- select.list(quantile.names.2, title="Model 2", multiple = FALSE)}
	if(length(quantile.2)==0 || is.null(quantile.2) || !quantile.2%in%quantile.names.2){
		stop("must provide valid quantile for model.obj.2")}
	IMP.2<-imp.extract.qrf(model.obj.2,imp.type=imp.type.2,ONEorTWO="model.obj.2")
	IMP.2<-IMP.2[,c("pred",quantile.2)]
	names(IMP.2)<-c("pred","imp")
	IMP.2<-IMP.2[order(IMP.2$imp),]
}
if(model.type.2=="CF"){
	if(is.null(imp.type.2)){imp.type.2<-1}
	if(!is.null(class.2)){
		warning("no class specific measure for CF models therefore 'class.2' ignored",immediate.=WARN.IM)
		class.2<-NULL}
	if(imp.type.2==2 && response.type.2!="binary"){
		warning("AUC-based variables importances only available for binary response models therefor importance type changed to 'imp.type.2=1'",immediate.=WARN.IM)
		imp.type.2<-1}
	IMP.2<-imp.extract.cf(	model.obj.2,imp.type=imp.type.2,mincriterion=cf.mincriterion.2, 
					conditional=cf.conditional.2,threshold=cf.threshold.2,nperm=cf.nperm.2)
}
if(model.type.2=="SGB"){
	if(is.null(imp.type.2)){imp.type.2<-2}
	if(response.type.2=="categorical" && imp.type.2==1){
		imp.type.2<-2
		warning("'response.type' categorical is not yet supported for permutation importance, relative influence will be used instead",immediate.=WARN.IM)}
	if(!is.null(class.2)){
		warning("no class specific measure for SGB models therefore 'class.2' ignored",immediate.=WARN.IM)}
	print(paste("imp.type.2:",imp.type.2))
	IMP.2<-imp.extract.sgb(model.obj.2,imp.type=imp.type.2)
}

#print("IMP.2")
#print(IMP.2)

IMP<-list(IMP1=IMP.1,IMP2=IMP.2)
############################################################################################
############################ Print Importance Types ########################################
############################################################################################

#print("print imp type")

if(is.null(class.1)){
	CLASS.1<-"Overall"
}else{
	CLASS.1<-paste("Class",class.1)
}

if(is.null(class.2)){
	CLASS.2<-"Overall"
}else{
	CLASS.2<-paste("Class",class.2)
}

if(model.type.1=="RF"){
	if(response.type.1%in%c("continuous")){
		if(imp.type.1==1){
			print(paste("model.obj.1 is a continuous RF model with", CLASS.1, "Importance measured by %IncMSE"))
			IMP.MEASURE.1<-"%IncMSE"}
		if(imp.type.1==2){
			print(paste("model.obj.1 is a continuous RF model with", CLASS.1, "Importance measured by IncNodePurity"))
			IMP.MEASURE.1<-"IncNodePurity"}
	}
	if(response.type.1%in%c("binary","categorical")){
		if(imp.type.1==1){
			print(paste("model.obj.1 is a", response.type.1, "RF model with", CLASS.1, "Importance measured by MeanDecreaseAccuracy"))
			IMP.MEASURE.1<-"MeanDecAccuracy"}
		if(imp.type.1==2){
			print(paste("model.obj.1 is a", response.type.1, "RF model with", CLASS.1, "Importance measured by MeanDecreaseGini"))
			IMP.MEASURE.1<-"MeanDecGini"}
	}
}

if(model.type.1=="QRF"){
	IMP.MEASURE.1<-paste(quantile.1)
	print(paste("model.obj.1 is continuous response QRF model with",IMP.MEASURE.1))}

if(model.type.1=="CF"){
	if(imp.type.1==1){
		if(cf.conditional.1 == FALSE){
			print(paste("model.obj.1 is a", response.type.1, "CF model with unconditional importance measured by MeanDecreaseAccuracy"))}
		if(cf.conditional.1 == TRUE){
			print(paste("model.obj.1 is a", response.type.1, "CF model with conditional importance measured by MeanDecreaseAccuracy"))}
		IMP.MEASURE.1<-"MeanDecAccuracy"	
	}
	if(imp.type.1==2){
		if(cf.conditional.1 == FALSE){
			print(paste("model.obj.1 is a", response.type.1, "CF model with unconditional importance measured by mean decrease in AUC"))}
		if(cf.conditional.1 == TRUE){
			print(paste("model.obj.1 is a", response.type.1, "CF model with conditional importance measured by mean decrease in AUC"))}
		IMP.MEASURE.1<-"MeanDecAUC"
	}
}

if(model.type.1=="SGB"){
	if(imp.type.1==1){
		if(response.type.1%in%c("continuous")){
			print("model.obj.1 is a continuous SGB model with importance measured by decrease in predictive performance with permutation")
			IMP.MEASURE.1<-"DecPredPerf"} 
		if(response.type.1%in%c("binary")){
			print("model.obj.1 is a binary SGB model with importance measured by decrease in predictive performance with permutation")
			IMP.MEASURE.1<-"DecPredPerf"}
	}
	if(imp.type.1==2){
		if(response.type.1%in%c("continuous")){
			print("model.obj.1 is a continuous SGB model with importance measured by decrease of squared error")
			IMP.MEASURE.1<-"DecSqError"} 
		if(response.type.1%in%c("binary")){
			print("model.obj.1 is a binary SGB model with importance measured by decrease in sum of squared error")
			IMP.MEASURE.1<-"DecSumSqError"}
		if(response.type.1%in%c("categorical")){
			print("model.obj.1 is a categorical SGB model with importance measured by decrease in sum of squared error")
			IMP.MEASURE.1<-"DecSumSqError"}
	}
}

if(model.type.2=="RF"){
	if(response.type.2%in%c("continuous")){
		if(imp.type.2==1){
			print(paste("model.obj.2 is a continuous RF model with", CLASS.2, "Importance measured by %IncMSE"))
			IMP.MEASURE.2<-"%IncMSE"}
		if(imp.type.2==2){
			print(paste("model.obj.2 is a continuous RF model with", CLASS.2, "Importance measured by IncNodePurity"))
			IMP.MEASURE.2<-"IncNodePurity"}
	}
	if(response.type.2%in%c("binary","categorical")){
		if(imp.type.2==1){
			print(paste("model.obj.2 is a", response.type.2, "RF model with", CLASS.2, "Importance measured by MeanDecreaseAccuracy"))
			IMP.MEASURE.2<-"MeanDecAccuracy"}
		if(imp.type.2==2){
			print(paste("model.obj.2 is a", response.type.2, "RF model with", CLASS.2, "Importance measured by MeanDecreaseGini"))
			IMP.MEASURE.2<-"MeanDecGini"}
	}
}

if(model.type.2=="QRF"){
	IMP.MEASURE.2<-paste(quantile.2)
	print(paste("model.obj.2 is continuous response QRF model with",IMP.MEASURE.2))}

if(model.type.2=="CF"){
	if(imp.type.2==1){
		if(cf.conditional.2 == FALSE){
			print(paste("model.obj.2 is a", response.type.2, "CF model with unconditional importance measured by MeanDecreaseAccuracy"))}
		if(cf.conditional.2 == TRUE){
			print(paste("model.obj.2 is a", response.type.2, "CF model with conditional importance measured by MeanDecreaseAccuracy"))}
		IMP.MEASURE.2<-"MeanDecAccuracy"	
	}
	if(imp.type.2==2){
		if(cf.conditional.2 == FALSE){
			print(paste("model.obj.1 is a", response.type.2, "CF model with unconditional importance measured by mean decrease in AUC"))}
		if(cf.conditional.2 == TRUE){
			print(paste("model.obj.1 is a", response.type.2, "CF model with conditional importance measured by mean decrease in AUC"))}
		IMP.MEASURE.2<-"MeanDecAUC"
	}
}

if(model.type.2=="SGB"){
	if(imp.type.2==1){
		if(response.type.2%in%c("continuous")){
			print("model.obj.2 is a continuous SGB model with importance measured by decrease in predictive performance with permutation")
			IMP.MEASURE.2<-"DecPredPerf"} 
		if(response.type.2%in%c("binary")){
			print("model.obj.2 is a binary SGB model with importance measured by decrease in predictive performance with permutation")
			IMP.MEASURE.2<-"DecPredPerf"}
	}
	if(imp.type.2==2){
		if(response.type.2%in%c("continuous")){
			print("model.obj.2 is a continuous SGB model with importance measured by decrease of squared error")
			IMP.MEASURE.2<-"DecSqError"} 
		if(response.type.2%in%c("binary")){
			print("model.obj.2 is a binary SGB model with importance measured by decrease in sum of squared error")
			IMP.MEASURE.2<-"DecSumSqError"}
		if(response.type.2%in%c("categorical")){
			print("model.obj.2 is a categorical SGB model with importance measured by decrease in sum of squared error")
			IMP.MEASURE.2<-"DecSumSqError"}
	}
}

#####################################################################################
############################ Scale Importace ########################################
#####################################################################################

#print("scale imp")

if(!scale.by %in% c("max","sum")){
	stop("scale.by must be either max or sum")}


IMP.1<-imp.scale(IMP.1,scale.by=scale.by)
IMP.2<-imp.scale(IMP.2,scale.by=scale.by)

if(scale.by=="sum"){
	MAX.IMP<-max(IMP.1$imp,IMP.2$imp)
	IMP.1$imp<-IMP.1$imp/MAX.IMP
	IMP.2$imp<-IMP.2$imp/MAX.IMP
}

if(!identical(sort(as.character(IMP.1$pred)),sort(as.character(IMP.2$pred)))){
	stop("models contain different predictors")
}

#if(!is.null(predList)){
#	if(!identical(sort(as.character(IMP.1$pred)),sort(predList))){
#		stop("predList contains different predictors than models")
#	}
#
#	if(!identical(sort(as.character(IMP.2$pred)),sort(predList))){
#		stop("predList contains different predictors than models")
#	}
#}

#####################################################################################
#################################### Sort Imp #######################################
#####################################################################################


if(!sort.by %in% c("predList","model.obj.1","model.obj.2")){
	stop("sort.by must be 'model.obj.1' or 'model.obj.2' or 'predList'")}



if(sort.by=="model.obj.1"){
	sort.by=IMP.1$pred
}else{
	if(sort.by=="model.obj.2"){
		sort.by=IMP.2$pred
	}else{
		if(sort.by=="predList"){
			sort.by=rev(predList)
		}
	}
}

IMP.1<-IMP.1[match(sort.by,IMP.1$pred),]
IMP.2<-IMP.2[match(sort.by,IMP.2$pred),]

#####################################################################################
################################### Make plot #######################################
#####################################################################################

	names.arg=IMP.1$pred
	NCHAR<-max(nchar(as.character(names.arg)))

#print("Starting graphics")
#if(!"none"%in%device.type){
for(i in 1:length(device.type)){

###################################################################
### Output filenames ###

initialize.device(	PLOTfn=PLOTfn,DEVICE.TYPE=device.type[i],
				res=res,device.width=device.width,device.height=device.height,
				units=units,pointsize=pointsize,cex=cex)

	op<-par(mar=par()$mar+c(0,(2*NCHAR/3)-2,0,0),cex=cex)

	barplot(-IMP.1$imp,horiz=TRUE,xlim=c(-1,1),las=1,col=col.1,axes=F,names.arg=names.arg,...)
	barplot(IMP.2$imp,horiz=TRUE,xlim=c(-1,1),las=1,col=col.2,add=TRUE,axes=F,...)

	IMP.TEXT<-paste(IMP.MEASURE.1,"  ",IMP.MEASURE.2)

	#mtext(IMP.TEXT,side=1, line=0, adj=.5, cex=cex)

	if(type.label){
		mtext(IMP.MEASURE.1,side=1, line=0, adj=0, cex=1*cex)
		mtext(IMP.MEASURE.2,side=1, line=0, adj=1, cex=1*cex)}

	mtext(model.name.1,side=1, line=1.5, adj=0, cex=1.3*cex)
	mtext(model.name.2,side=1, line=1.5, adj=1, cex=1.3*cex)

	par(op)
if(!device.type[i]%in%c("default","none")){dev.off()}
}
#}
#############################################################################

return(IMP)

#############################################################################
}


