#############################################
# arf3DS4 S4 MODEL FUNCTIONS  	     		#
# Wouter D. Weeda							#
# University of Amsterdam					#
#############################################

#[CONTAINS]
#newModel				[user]
#updateModelNames	
#saveModel				[user]
#saveModelBin			[user]
#saveModelBinSimple		[user]
#loadOptions			[user]
#saveOptions			[user]
#loadStart				[user]
#saveStart				[user]
#loadStart				[user]
#loadModel				[user]
#updateClass
#clearWarnings			[user]
#showModels				[user]


newModel <- 
function(modelname='defaultmodel',regions=1,subject='',condition='',type=c('gauss','simple'),options=new('options'),overwrite=T,experiment=NULL) 
#newModel makes a modeldirectory based on data and experiment information and a modelname
{
	#check experiment
	if(is.null(experiment)) {
		experiment <- try(get('.experiment',envir=.arfInternal),silent=T)
		if(attr(experiment,'class')=='try-error') stop('Experiment not loaded. Run loadExp first.')
	}
	
	#set separator
	sp <- .Platform$file.sep
	
	#load data from a condition
	if(subject=='') subject <- .experiment.subject.names(experiment)[1]
	if(condition=='') condition <- .experiment.condition.names(experiment)[1]
	arfdata <- loadRda(paste(.experiment.path(experiment),sp,.experiment.subjectDir(experiment),sp,subject,sp,.experiment.conditionDir(experiment),sp,condition,sp,.experiment.dataDir(experiment),sp,.experiment.dataRda(experiment),sep=''))
	
	#checks
	if(!class(arfdata)=='data') stop('Input must be of class \'data\'')
	if(modelname=='') stop('Modelname cannot be empty')
		
	#check if Averages exist (else create)
	if(!file.exists(.data.avgdatfile(arfdata)) | !file.exists(.data.avgWfile(arfdata)) | !file.exists(.data.avgtstatFile(arfdata))) {
		if(overwrite) {
			arfdata <- createAverages(arfdata,experiment)
		} else {
			warning('Averages do not exist, and overwrite is FALSE! Things might go wrong.')
		}
	} 
	
	#make new modelobject
	model <- new('model',arfdata)
		
	#set modelname
	.model.modelname(model) <- modelname 
			
	#path to modeldir
	path <- paste(.data.fullpath(arfdata),sp,.experiment.modelDir(experiment),sp,modelname,sep='')
	
	if(file.exists(path) & overwrite==F) stop('Model directory already exists')
	
	dir.create(path,showWarnings=F)
	.model.modelpath(model) <- path
	
	#path to modeldatadir
	dir.create(paste(path,sp,.experiment.modeldatDir(experiment),sep=''),showWarnings=F)
	.model.modeldatapath(model) <- paste(path,sp,.experiment.modeldatDir(experiment),sep='')
	
	#set filenames
	.model.residualFile(model) <- .experiment.residualFile(experiment)
	.model.derivativeFile(model) <- .experiment.derivativeFile(experiment)
	.model.weightFile(model) <- .experiment.weightFile(experiment)
	.model.modelFile(model) <- .experiment.modelRda(experiment)
	.model.optionsFile(model) <- .experiment.optionsRda(experiment)
	.model.modelDataFile(model) <- .experiment.modelDataFile(experiment)
	.model.startFile(model) <- .experiment.startRda(experiment)
	.model.logFile(model) <- .experiment.logFile(experiment)
	
	#set modeltype
	type = match.arg(type)
	.model.modeltype(model) <- type
	if(type=='simple') .model.params(model) <- 5
	if(type=='gauss') .model.params(model) <- 10
		
	#set number of regions
	.model.regions(model) <- regions
	
	#save model and options File and startvec
	startval <- .model.startval(model) <- rep(.options.start.vector(options),.model.regions(model))
	save(startval,file=paste(path,sp,.model.startFile(model),sep=''))
	save(options,file=paste(path,sp,.model.optionsFile(model),sep=''))
	save(model,file=paste(path,sp,.model.modelFile(model),sep=''))
	
	
	#updateModelNames in on dir up
	updateModelNames(dirname(.model.modelpath(model)))
		
	return(model)
}


updateModelNames <- 
function(path) 
#update ModelNames in a ModelNamesFile
{
		
	#list all dirs in path (minus the modelnames file)
	existingfiles <- list.files(path,full.names=F)
	filename <- list.files(path,'.Rda',full.names=T)
	existingfiles <- existingfiles[-grep('.Rda',existingfiles)]

	#save modelnames
	save(existingfiles,file=filename)	
	
}

saveModel <- 
function(arfmodel) save(arfmodel,file=paste(.model.modelpath(arfmodel),.Platform$file.sep,.model.modelFile(arfmodel),sep=''))
#save the model to the model.Rda

saveModelBin <- 
function(arfmodel,type=c('full','pos','neg','fpn','separate','sig')) 
#save the modelBinary
{
	
	
	#match type
	type <- match.arg(type)
	pos=neg=FALSE
	full=T
	if(type=='full') full=T
	if(type=='neg') neg=T
	if(type=='pos') pos=T
	if(type=='fpn') full=neg=pos=T
		
	sp = .Platform$file.sep
	
	#get Header info from avgdatfile
	headinf <- readHeader(getFileInfo(.model.avgdatfile(arfmodel)))
	
	#set fullpaths
	.nifti.header.fullpath(headinf) <- .model.modeldatapath(arfmodel)
	
	#write the Data to the modelNiftiFile
	if(full) {
		.nifti.header.filename(headinf) <- .model.modelDataFile(arfmodel)
		.model.fullmodelDataFile(arfmodel) <- headToName(headinf)
		writeData(headinf,.C('gauss',as.double(.model.estimates(arfmodel)),as.integer(.model.regions(arfmodel)*.model.params(arfmodel)),as.integer(.nifti.header.dims(headinf)[2]),as.integer(.nifti.header.dims(headinf)[3]),as.integer(.nifti.header.dims(headinf)[4]),as.double(numeric(.nifti.header.dims(headinf)[2]*.nifti.header.dims(headinf)[3]*.nifti.header.dims(headinf)[4])))[[6]])
	}
	
	if(pos)	{
		theta = matrix(.model.estimates(arfmodel),.model.params(arfmodel))
		which_neg = which(theta[10,]<0)
		
		if(length(which_neg)>0) theta = theta[,-which_neg]
		regs = ncol(theta)
		thetavec = as.vector(theta)
		
		.nifti.header.filename(headinf) <- paste(.model.modelDataFile(arfmodel),'_pos',sep='')
		writeData(headinf,.C('gauss',as.double(thetavec),as.integer(regs*.model.params(arfmodel)),as.integer(.nifti.header.dims(headinf)[2]),as.integer(.nifti.header.dims(headinf)[3]),as.integer(.nifti.header.dims(headinf)[4]),as.double(numeric(.nifti.header.dims(headinf)[2]*.nifti.header.dims(headinf)[3]*.nifti.header.dims(headinf)[4])))[[6]])
	}
	
	if(neg)	{
		theta = matrix(.model.estimates(arfmodel),.model.params(arfmodel))
		which_pos = which(theta[10,]>0)
		if(length(which_pos)>0) theta = theta[,-which_pos]
		regs = ncol(theta)
		thetavec = as.vector(theta)
		
		.nifti.header.filename(headinf) <- paste(.model.modelDataFile(arfmodel),'_neg',sep='')
		writeData(headinf,.C('gauss',as.double(thetavec),as.integer(regs*.model.params(arfmodel)),as.integer(.nifti.header.dims(headinf)[2]),as.integer(.nifti.header.dims(headinf)[3]),as.integer(.nifti.header.dims(headinf)[4]),as.double(numeric(.nifti.header.dims(headinf)[2]*.nifti.header.dims(headinf)[3]*.nifti.header.dims(headinf)[4])))[[6]])
	}
	
	if(type=='separate') {
		theta = matrix(.model.estimates(arfmodel),.model.params(arfmodel))
		regs = 1
		
		for(i in 1:ncol(theta)) {
			thetavec = as.vector(theta[,i])
			.nifti.header.filename(headinf) <- paste(.model.modelDataFile(arfmodel),'_region',i,sep='')
			writeData(headinf,.C('gauss',as.double(thetavec),as.integer(regs*.model.params(arfmodel)),as.integer(.nifti.header.dims(headinf)[2]),as.integer(.nifti.header.dims(headinf)[3]),as.integer(.nifti.header.dims(headinf)[4]),as.double(numeric(.nifti.header.dims(headinf)[2]*.nifti.header.dims(headinf)[3]*.nifti.header.dims(headinf)[4])))[[6]])
		}
		
	}
	
	if(type=='sig') {
		theta = matrix(.model.estimates(arfmodel),.model.params(arfmodel))
		wald = .model.wald(arfmodel)
		
		if(length(wald)!=0) {
			
			which_sig = which(.wald.pvalues(wald)[,4]<.05 & .wald.pvalues(wald)[,5]<.05)
		
			if(length(which_sig)>0) theta = theta[,which_sig]
			regs = ncol(theta)
			thetavec = as.vector(theta)
			
			.nifti.header.filename(headinf) <- paste(.model.modelDataFile(arfmodel),'_sig',sep='')
			writeData(headinf,.C('gauss',as.double(thetavec),as.integer(regs*.model.params(arfmodel)),as.integer(.nifti.header.dims(headinf)[2]),as.integer(.nifti.header.dims(headinf)[3]),as.integer(.nifti.header.dims(headinf)[4]),as.double(numeric(.nifti.header.dims(headinf)[2]*.nifti.header.dims(headinf)[3]*.nifti.header.dims(headinf)[4])))[[6]])	
			
			
		} else warning('Wald statistics not calculated! No models saved.')
		
		
	}
	
	
	return(arfmodel)
	
}

saveModelBinSimple <- 
function(arfmodel) 
#save the modelBinary
{
	
	#get Header info from avgdatfile
	headinf <- readHeader(getFileInfo(.model.avgdatfile(arfmodel)))

	#set fullpaths
	.nifti.header.fullpath(headinf) <- .model.modeldatapath(arfmodel)
	.nifti.header.filename(headinf) <- .model.modelDataFile(arfmodel)
	.model.fullmodelDataFile(arfmodel) <- headToName(headinf)
	
	#write the Data to the modelNiftiFile
	writeData(headinf,.C('simplegauss',as.double(.model.estimates(arfmodel)),as.integer(.model.regions(arfmodel)*.model.params(arfmodel)),as.integer(.nifti.header.dims(headinf)[2]),as.integer(.nifti.header.dims(headinf)[3]),as.integer(.nifti.header.dims(headinf)[4]),as.double(numeric(.nifti.header.dims(headinf)[2]*.nifti.header.dims(headinf)[3]*.nifti.header.dims(headinf)[4])))[[6]])
	
	return(arfmodel)
	
}


loadOptions <- 
function(arfmodel) return(loadRda(paste(.model.modelpath(arfmodel),.Platform$file.sep,.model.optionsFile(arfmodel),sep='')))
#loadOptions loads the options object

saveOptions <- 
function(options,arfmodel) save(options,file=paste(.model.modelpath(arfmodel),.Platform$file.sep,.model.optionsFile(arfmodel),sep=''))
#loadOptions loads the options object


loadStart <- 
function(arfmodel) return(loadRda(paste(.model.modelpath(arfmodel),.Platform$file.sep,.model.startFile(arfmodel),sep='')))
#loadStart loads the Start object

saveStart <- 
function(startval,arfmodel) save(startval,file=paste(.model.modelpath(arfmodel),.Platform$file.sep,.model.startFile(arfmodel),sep=''))  	
#save startingvalues

loadModel <- 
function(modelname,subject=NA,condition,experiment=NULL) 
#load a model based on subject and conditions or on a mnames object
{

	#check experiment
	if(is.null(experiment)) {
		experiment <- try(get('.experiment',envir=.arfInternal),silent=T)
		if(attr(experiment,'class')=='try-error') stop('Experiment not loaded. Run loadExp first.')
	}
	
	if(class(modelname)=='mnames') {
		if(is.na(subject)) num = 1 else num = subject
		experiment = .mnames.experiment(modelname)
		subject = .mnames.subject(modelname)
		condition = .mnames.condition(modelname)
		modelname = .mnames.mnames(modelname)[num]
	}
	
	sp <- .Platform$file.sep
	modname = paste(.experiment.path(experiment),sp,.experiment.subjectDir(experiment),sp,subject,sp,.experiment.conditionDir(experiment),sp,condition,sp,.experiment.modelDir(experiment),sp,modelname,sp,.experiment.modelRda(experiment),sep='')
	mod = loadRda(modname)
	
	return(mod)
}

updateClass <- 
function(object,...) 
#update a modelclass with elements
{

	new_object <- new(class(object),object,...)
	return(new_object)
	
}

showModels <-
function(subject,condition,experiment=NULL)
#show all models in for a subject and condition
{
	
	#check experiment
	if(is.null(experiment)) {
		experiment <- try(get('.experiment',envir=.arfInternal),silent=T)
		if(attr(experiment,'class')=='try-error') stop('Experiment not loaded. Run loadExp first.')
	}
	
	sp <- .Platform$file.sep
	modname = paste(.experiment.path(experiment),sp,.experiment.subjectDir(experiment),sp,subject,sp,.experiment.conditionDir(experiment),sp,condition,sp,.experiment.modelDir(experiment),sp,.experiment.modelnamesRda(experiment),sep='')
	mod = loadRda(modname)

	mnames = new('mnames')
	
	.mnames.experiment(mnames) <- experiment
	.mnames.subject(mnames) <- subject
	.mnames.condition(mnames) <- condition
	.mnames.mnames(mnames) <- mod
	
	return(mnames)
	
}

clearWarnings <- 
function(arfmodel,resetValid=T) 
#clears the warnings from a model and reset the valid object then save the object
{ 
	yn = readline('This will clear all warnings! Are you sure (y/n) ')
	
	if(substr(yn,1,1)=='y') {
		.model.warnings(arfmodel) = character(0)
		if(resetValid) .model.valid(arfmodel) = T
		saveModel(arfmodel)
	}
	return(arfmodel)
}

