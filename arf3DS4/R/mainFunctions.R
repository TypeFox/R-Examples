#############################################
# arf3DS4 MAIN FUNCTIONS                 	#
# Wouter D. Weeda							#
# University of Amsterdam					#
#############################################

#[CONTAINS]
#processModel		[user]
#processSeed		[user]
#searchRange		[user]
#minBIC				[user]

processModel <- 
function(arfmodel,options=loadOptions(arfmodel),dat=readData(.model.avgdatfile(arfmodel)),weights=readData(.model.avgWfile(arfmodel)),pr=T,printlevel=0,try.silen=T) 
#processModel does an entire arf-fitting pipeline 
{

	.model.warnings(arfmodel) <- character(0)
	.model.convergence(arfmodel) <- character(0)
	
	#fit model
	arfmodel <- fitModel(arfmodel,options=options,dat=dat,weights=weights,printlevel=printlevel,try.silen=try.silen)
	
	#check for bounded parameters
	arfmodel <- checkSolution(arfmodel,options,dat,thres=8)
	
	if(.model.valid(arfmodel)) {
		
		if(pr) cat('\n calculating variance/covariance matrix...')
		arfmodel = varcov(arfmodel)
		if(.model.valid(arfmodel)) {if(pr) {cat('ok\n')}} else {if(pr) {cat('fail\n')}} 
		
		if(.model.valid(arfmodel)) {
			if(pr) cat(' calculating wald statistics...')
			arfmodel = wald(arfmodel)
			if(.model.valid(arfmodel)) {if(pr) {cat('ok\n')}} else {if(pr) {cat('fail\n')}} 
		} else {
			if(pr) cat(' wald statistics not calculated.\n')
			
		}
	
	} 
	cat('\n')
	
	#save the model
	saveModel(arfmodel)
	
	if(pr) show(arfmodel)
	
	if(pr) cat(' \nmodelfile saved in \'',.model.modelpath(arfmodel),.Platform$file.sep,.model.modelFile(arfmodel),'\'\n',sep='')
	if(pr) cat(' arf process stopped at',as.character(Sys.time()),'\n\n')
	
	return(invisible(arfmodel))
}

processSeed <-
function(modelname='defaultmodel',seedreg,subject='',condition='',startmethod=c('default','simple'),grad=NULL,bound=NULL,pval=NULL,options=new('options'),pr=T,printlevel=0,try.silen=T,overwrite=T,experiment=NULL)
#process a sequence based on a seed number of regions
{
	
	if(is.null(experiment)) {
		experiment <- try(get('.experiment',envir=.arfInternal),silent=T)
		if(attr(experiment,'class')=='try-error') stop('Experiment not loaded. Run loadExp first.')
	}
	
	if(pr) cat('[',modelname,'] with',seedreg,'seed region(s) started at',as.character(Sys.time()),'\n\n')
	
	#make the model
	full_model = newModel(paste('full_',modelname,sep=''),regions=seedreg,subject=subject,condition=condition,type='gauss',options=options,overwrite=overwrite,experiment=experiment)
		
	#use simple starts or defaults
	startmethod = match.arg(startmethod[1],c('default','simple'))
	if(startmethod=='simple') {	
		#run a simple model
		.options.start.method(options) = 'rect'
		simple_model = newModel(paste('simple_',modelname,sep=''),regions=seedreg,subject=subject,condition=condition,type='simple',options=options,overwrite=overwrite,experiment=experiment)
		
		if(pr) cat(as.character(Sys.time()),'fitting simple model...')
		simple_model = fitModel(simple_model)
		if(pr) 	if(.model.valid(simple_model))	cat('ok\n') else cat('fail\n')
		
		if(.model.valid(simple_model)) {
			.model.startval(full_model) = .model.estimates(simple_model)
			.model.startval(full_model)[.model.startval(full_model)==0]=.01
			saveModel(full_model)
		} else if(pr) cat('simple model not valid, using normal starting values\n') 
	} 
	
	#fit the model
	full_model = fitModel(full_model,options=options)
	full_model <- checkSolution(full_model,options,thres=bound)
	
	#prune the model
	model = pruneModel(full_model,modelname,subject,condition,grad,bound,pval,options,overwrite,experiment)
		
	return(model)
	
}


fitRange <-
function(subject,condition,range,options=new('options'),modelprefix='searchmodel',modeltype=c('gauss','simple'),experiment=NULL)
#search a range of models for the one with the lowest BIC
{
	#experiment
	if(is.null(experiment)) {
		experiment <- try(get('.experiment',envir=.arfInternal),silent=T)
		if(attr(experiment,'class')=='try-error') stop('Experiment not loaded. Run loadExp first.')
	}
	
	#check for correct startmethod
	if(.options.start.method(options)=='load' | .options.start.method(options)=='use') stop('Search range only works with rect startingvalues')
	
	#matches
	modeltype = match.arg(modeltype[1],c('gauss','simple'))
		
		
	#set new sequence
	sequence = new('sequence')
	.sequence.regions(sequence) = range
	.sequence.fit(sequence) = .sequence.minimum(sequence) = as.numeric(rep(NA,length(.sequence.regions(sequence))))
	.sequence.valid(sequence) = .sequence.best(sequence) = rep(FALSE,length(.sequence.regions(sequence)))
	
	#perform fit loop
	p=1
	for(numreg in range) {
		model = newModel(paste(modelprefix,'_',numreg,'regs',sep=''),numreg,subject,condition,type=modeltype,options=options,experiment=experiment)
		model = fitModel(model,options=options)
		model <- checkSolution(model,options,thres=8)
		
		.sequence.fit(sequence)[p] = .model.fit(model)[1]
		.sequence.minimum(sequence)[p] = .model.minimum(model)
		
		if(.model.valid(model)) {
			.sequence.mnames(sequence)[p] = .model.modelname(model)
			.sequence.valid(sequence)[p] = .model.valid(model)
		}
		p=p+1
	}
				
	
	#which of the valid models has the minimum
	validfits = .sequence.fit(sequence)[which(.sequence.valid(sequence)==TRUE)]
	minimum = which.min(validfits)
	.sequence.best(sequence)[which(.sequence.valid(sequence)==TRUE)][minimum] = TRUE
	
	return(sequence)
		
}


minBIC <-
function(subject,condition) 
#show the minimal BIC returns a sequence object
{
	#load modelnames
	modobject <- showModels(subject,condition)
	modnames <- .mnames.mnames(modobject)
	
	#make new sequence object
	sequence = new('sequence')
	.sequence.mnames(sequence) <- modnames
	.sequence.fit(sequence) <-  .sequence.minimum(sequence) <- as.numeric(rep(NA,length(modnames)))
	.sequence.valid(sequence) <- .sequence.best(sequence) <- rep(FALSE,length(modnames))
	.sequence.current(sequence) <- as.numeric(NA)

	#fill in sequence object
	for(i in 1:length(modnames)) {
		mod = loadModel(modobject,i)
		.sequence.regions(sequence)[i] = .model.regions(mod)
		if(.model.valid(mod)) {
			.sequence.fit(sequence)[i] = .model.fit(mod)[1] 
			.sequence.minimum(sequence)[i] = .model.minimum(mod)
			.sequence.valid(sequence)[i]=TRUE
		} else {
			.sequence.valid(sequence)[i]=FALSE
		}
	}
	
	#which of the valid models has the minimum
	validfits = .sequence.fit(sequence)[which(.sequence.valid(sequence)==TRUE)]
	minimum = which.min(validfits)
	.sequence.best(sequence)[which(.sequence.valid(sequence)==TRUE)][minimum] = TRUE

	return(sequence)
}

