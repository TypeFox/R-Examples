#############################################
# arf3DS4 S4 FITMODEL FUNCTIONS				#
# Wouter D. Weeda							#
# University of Amsterdam					#
#############################################

#[CONTAINS]
#fitModel					[user]
#fitModelOptim
#fitSimpleModelOptim
#pruneModel					

fitModel <- 
function(arfmodel,options=loadOptions(arfmodel),dat=readData(.model.avgdatfile(arfmodel)),weights=readData(.model.avgWfile(arfmodel)),printlevel=0,try.silen=T) 
# fitModel is a wrapper for NLM and optim based on the options
{
	
	cat('[',.model.modelname(arfmodel),']\n')
	cat(' arf process for data',.model.name(arfmodel),'started',as.character(Sys.time()),'\n')
	cat(' fitting',.model.regions(arfmodel),'region(s)\n')
	
	if(.model.modeltype(arfmodel)=='simple') {
		type='simple'
		if(.model.params(arfmodel)!=5) stop('Modeltype - parameter mismatch!')
	} else {
		if(.model.modeltype(arfmodel)=='gauss') { 
			type='gauss'
			if(.model.params(arfmodel)!=10) stop('Modeltype - parameter mismatch!')
		} else stop('Type of model to fit cannot be found, check modeltype slot of arfmodel object.')
	}
		
	if(.options.start.method(options)=='rect') {
		if(type=='gauss') arfmodel <- determineStartRect(arfmodel)
		if(type=='simple') arfmodel <- determineStartRectSimple(arfmodel) 	
	}
	
	if(.options.start.method(options)=='simple') {
		startmodel <- arfmodel
		.model.modeltype(startmodel)='simple'
		.model.params(startmodel)=5
		startmodel <- determineStartRectSimple(arfmodel)
		startmodel = fitSimpleModelOptim(startmodel,options=options,dat=dat,weights=weights,printlevel=printlevel,try.silen=try.silen)
		.model.startval(arfmodel) <- .model.estimates(startmodel)
	}
	
	if(.options.start.method(options)=='load') .model.startval(arfmodel) <- loadStart(arfmodel)
		
	if(.options.min.routine(options)[1]=='optim') {
		if(type=='simple')	arfmodel = fitSimpleModelOptim(arfmodel,options=options,dat=dat,weights=weights,printlevel=printlevel,try.silen=try.silen) 
		if(type=='gauss')	arfmodel = fitModelOptim(arfmodel,options=options,dat=dat,weights=weights,printlevel=printlevel,try.silen=try.silen) 
	}
	
	#after fitting
	cat(' ',.model.convergence(arfmodel),'\n',sep='')
	cat(' <modelfit>\n')
	cat('  minimum:',round(.model.minimum(arfmodel)),'\n')
	cat('    BIC  :',round(.model.fit(arfmodel)[1]),'\n')
	cat('    RMSEA:',round(.model.fit(arfmodel)[2],1),'\n')
	return(arfmodel)	
}


fitModelOptim <- 
function(arfmodel,options=loadOptions(arfmodel),dat=readData(.model.avgdatfile(arfmodel)),weights=readData(.model.avgWfile(arfmodel)),printlevel=0,try.silen=T) 
# fitModelOptim calls the minimization routine (OPTIM)
{

	#set separator
	sp <- .Platform$file.sep
	
	#set modelobjects
	.options.min.routine(options)[1] <- 'optim'
	if(is.na(.options.min.routine(options)[2])) .options.min.routine(options)[2]='vpv'
	if(.model.modeltype(arfmodel)!='gauss') stop('Called fit to Gauss model for non-gauss model object.')
	if(.model.params(arfmodel)!=10) stop('Modeltype - parameter mismatch!')
	.model.valid(arfmodel) <- TRUE
	
	#start_time
	st_time <- Sys.time()
		
	#check if averages exist
	if(!file.exists(.model.avgdatfile(arfmodel))) stop('Averages do not exist, please run createAverages')
	if(!file.exists(.model.avgWfile(arfmodel))) stop('Averages do not exist, please run createAverages')
	
	#clear the warnings and deriv + residualfilres
	.model.warnings(arfmodel) <- character(0)
	if(file.exists(paste(.model.modeldatapath(arfmodel),sp,.model.residualFile(arfmodel),sep=''))) file.remove(paste(.model.modeldatapath(arfmodel),sp,.model.residualFile(arfmodel),sep=''))
	if(file.exists(paste(.model.modeldatapath(arfmodel),sp,.model.derivativeFile(arfmodel),sep=''))) file.remove(paste(.model.modeldatapath(arfmodel),sp,.model.derivativeFile(arfmodel),sep=''))
	
	#set analyticalgrad options
	if(.options.min.analyticalgrad(options)) {
		gradfunc=gradient.gauss		
		#if(.options.min.routine(options)[2]=='vpv') gradfunc=gradient.gauss else gradfunc=gradient.gauss.rpr
		angrad=FALSE
	} else {
		gradfunc=NULL
		angrad=FALSE
	}
	
	#set which method of model calculation to use 
	#if(.options.min.routine(options)[2]=='vpv') sumsofsquares = ssq.gauss else sumsofsquares = ssq.gauss.rpr
	sumsofsquares = ssq.gauss	

	#set boundaries in L-BFGS-B mode
	if(length(.options.opt.lower(options))==1 | length(.options.opt.upper(options))==1) {
		lowbound=-Inf
		upbound=Inf
	} else {
		#set location to maximal dim
		max_loc = c(.fmri.data.dims(dat)[2],.fmri.data.dims(dat)[3],.fmri.data.dims(dat)[4])
	
		#set width parameters to maxdim divided by the value given in the options
		max_width =  c(.fmri.data.dims(dat)[2],.fmri.data.dims(dat)[3],.fmri.data.dims(dat)[4]) / c(.options.opt.upper(options)[4],.options.opt.upper(options)[5],.options.opt.upper(options)[6])
		
		upbound = rep(c(max_loc,max_width,.options.opt.upper(options)[7:10]),.model.regions(arfmodel))
		lowbound = rep(.options.opt.lower(options),.model.regions(arfmodel))
	}

	#check startingvalues
	arfmodel <- validStart(arfmodel)
		
	#final check before fit
	if(!.model.valid(arfmodel)) {
		saveModel(arfmodel)	
		return(arfmodel)
	}
	
	#make progressElements
	progress = newProgressElement(arfmodel,options,lowbound,upbound)	

	#assign global variables
	assign('.gradient_latest',NA,envir=.arfInternal)
	assign('.theta_latest',NA,envir=.arfInternal)
	assign('.arf_error',numeric(0),envir=.arfInternal)
	
	#runoptim	
	optim.output <- try(suppressWarnings(optim(
						.model.startval(arfmodel),
						sumsofsquares,
						gradfunc,
						lower=lowbound,
						upper=upbound,
						datavec=.fmri.data.datavec(dat)[1:(.fmri.data.dims(dat)[2]*.fmri.data.dims(dat)[3]*.fmri.data.dims(dat)[4])],
						weightvec=.fmri.data.datavec(weights)[1:(.fmri.data.dims(weights)[2]*.fmri.data.dims(weights)[3]*.fmri.data.dims(dat)[4])],
						brain=.model.mask(arfmodel),
						np=.model.regions(arfmodel)*.model.params(arfmodel),
						dimx=.fmri.data.dims(dat)[2],
						dimy=.fmri.data.dims(dat)[3],
						dimz=.fmri.data.dims(dat)[4],
						ss_data=.model.ss(arfmodel),
						analyticalgrad=angrad,
						method=.options.opt.method(options),
						progress=progress,
						control=list(trace=printlevel,maxit=.options.min.iterlim(options)),
						hessian=T
					)),silent=try.silen)

	#end_time
	en_time <- Sys.time()

	# check for internal errors and set relevant arf model values
	if(is.null(attr(optim.output,'class'))) {
		if(optim.output$convergence==0) .model.convergence(arfmodel) <- paste('[optim] Optim converged in ',optim.output$counts[1],' iterations.',sep='')
		if(optim.output$convergence==1) {.model.convergence(arfmodel) <- '[optim] Iteration limit exceeded. No convergence.';.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),paste('[min] optim did not converge.',sep=''))}
		if(optim.output$convergence==10) .model.convergence(arfmodel) <- '[optim] Degeneracy of the Nelder-Mead Simplex'
		if(optim.output$convergence==51) .model.convergence(arfmodel) <- paste('[optim] BFGS raises warning:',gsub('\n','',optim.output$message),sep='')
		if(optim.output$convergence==52) .model.convergence(arfmodel) <-  paste('[optim] BFGS raises error:',gsub('\n','',optim.output$message),sep='')
		
		if(optim.output$convergence == 0) .model.valid(arfmodel) <- TRUE else .model.valid(arfmodel) <- FALSE
	
		#set model objects
		.model.minimum(arfmodel) <- optim.output$value
		.model.estimates(arfmodel) <- optim.output$par
		.model.hessian(arfmodel) <- optim.output$hessian
		.model.iterates(arfmodel) <- optim.output$counts[1]
		.model.sandwichmethod(arfmodel) <- .options.sw.type(options)
		.model.proctime(arfmodel)[1,1] <- as.numeric(difftime(en_time,st_time,units='sec'))

		if(.options.min.analyticalgrad(options)) .model.gradient(arfmodel) <- gradient.gauss(.model.estimates(arfmodel),.fmri.data.datavec(dat)[1:(.fmri.data.dims(dat)[2]*.fmri.data.dims(dat)[3]*.fmri.data.dims(dat)[4])],.fmri.data.datavec(weights)[1:(.fmri.data.dims(weights)[2]*.fmri.data.dims(weights)[3]*.fmri.data.dims(dat)[4])],.model.mask(arfmodel),.model.regions(arfmodel)*.model.params(arfmodel),.fmri.data.dims(dat)[2],.fmri.data.dims(dat)[3],.fmri.data.dims(dat)[4],.model.ss(arfmodel),analyticalgrad=T,progress=progress)
		
		if(.model.valid(arfmodel)) {
			#save the ModelBinary
			arfmodel <- saveModelBin(arfmodel)
			
			#save the weights in a binary file
			makeWeightsBin(arfmodel)	
			
			#make Derivatives 
			makeDerivs(arfmodel)
			
			#create residuals
			makeResiduals(arfmodel)
			
			#if(.options.min.analyticalgrad(options)) {
			#	df_fn <- paste(.model.modeldatapath(arfmodel),.Platform$file.sep,.model.derivativeFile(arfmodel),sep='')
			#	w_fn <- paste(.model.modeldatapath(arfmodel),.Platform$file.sep,.model.weightFile(arfmodel),sep='')
			#	n = .fmri.data.dims(weights)[2]*.fmri.data.dims(weights)[3]*.fmri.data.dims(weights)[4]
			#	p = .model.regions(arfmodel)*.model.params(arfmodel)
			#	hessian <- try(.C('approxHessian',as.integer(p),as.integer(.model.n(arfmodel)),as.character(df_fn),as.character(w_fn),as.double(numeric(p*p))),silent=try.silen)
			#	
			#	if(is.null(attr(hessian,'class'))) {
			#		hessian <- hessian[[5]]
			#		dim(hessian) = c(p,p)
			#		.model.hessian(arfmodel) <- hessian 	
			#	} else {
			#		.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),' [min] could not approximate Hessian with analytical gradient.')
			#		.model.valid(arfmodel) <- FALSE
			#		
			#	}
			#}
			
			#caluclate fits
			arfmodel = BIC(arfmodel,options=options)
			arfmodel = RMSEA(arfmodel,options=options)
			
		} else .model.warnings(arfmodel) <- c(.model.warnings(arfmodel),paste('[min] optim did not converge.',sep=''))
			
		
	} else {
		.model.estimates(arfmodel) <- get('.theta_latest',envir=.arfInternal)
		.model.gradient(arfmodel) <- get('.gradient_latest',envir=.arfInternal)
		pers <- get('.arf_error',envir=.arfInternal)
		
		if(length(pers)==0) {
			.model.convergence(arfmodel) <- '[optim] Internal error, no convergence.'
			.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),paste('[min] optim internal error: ',gsub('\n','',optim.output),sep=''))
			.model.proctime(arfmodel)[1,1] <- as.numeric(difftime(en_time,st_time,units='sec'))
			.model.valid(arfmodel) <- FALSE
		} else {
			if(pers$errtype=='persbound') {
				.model.convergence(arfmodel) <- '[arf] Persistent boundary error, no convergence.'
				.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),paste('[min] Persistent boundaries on regions: ',paste(pers$data,collapse=','),sep=''))
				.model.proctime(arfmodel)[1,1] <- as.numeric(difftime(en_time,st_time,units='sec'))
				.model.valid(arfmodel) <- FALSE
				arfmodel <- saveModelBin(arfmodel)
			}
			if(pers$errtype=='iterlim') {
				.model.convergence(arfmodel) <- '[arf] Iteration limit reached (ARF), no convergence.'
				.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),paste('[min] Iteration limit reached after ',paste(pers$data,collapse=','),' iterations.',sep=''))
				.model.proctime(arfmodel)[1,1] <- as.numeric(difftime(en_time,st_time,units='sec'))
				.model.valid(arfmodel) <- FALSE
				arfmodel <- saveModelBin(arfmodel)
			}
		}
	}
		
	if(!.model.valid(arfmodel)) .model.warnings(arfmodel) <- c(.model.warnings(arfmodel),.model.convergence(arfmodel)) 
	
	#clean up
	rm('.theta_latest','.gradient_latest','.gradit','.gradval','.objit','.bounded','.arf_error','.oldobj',envir=.arfInternal)
		
	#save the modelInfo
	saveModel(arfmodel)
	
	#return arf model object	
	return(invisible(arfmodel))
}

fitSimpleModelOptim <- 
function(arfmodel,options=loadOptions(arfmodel),dat=readData(.model.avgdatfile(arfmodel)),weights=readData(.model.avgWfile(arfmodel)),printlevel=0,try.silen=T) 
# fitModelOptim calls the minimization routine (OPTIM)
{
	
	#set separator
	sp <- .Platform$file.sep
	
	#set routine to optim
	.options.min.routine(options)[1] <- 'optim'
	if(is.na(.options.min.routine(options)[2])) .options.min.routine(options)[2]='vpv'
	if(.model.modeltype(arfmodel)!='simple') stop('Called fit to simple model for non-simple model object.')
	if(.model.params(arfmodel)!=5) stop('Modeltype - parameter mismatch!')
	.model.valid(arfmodel) <- TRUE

	
	#start_time
	st_time <- Sys.time()
	
	#check if averages exist
	if(!file.exists(.model.avgdatfile(arfmodel))) stop('Averages do not exist, please run createAverages')
	if(!file.exists(.model.avgWfile(arfmodel))) stop('Averages do not exist, please run createAverages')
	
	#clear the warnings and deriv + residualfilres
	.model.warnings(arfmodel) <- character(0)
	if(file.exists(paste(.model.modeldatapath(arfmodel),sp,.model.residualFile(arfmodel),sep=''))) file.remove(paste(.model.modeldatapath(arfmodel),sp,.model.residualFile(arfmodel),sep=''))
	if(file.exists(paste(.model.modeldatapath(arfmodel),sp,.model.derivativeFile(arfmodel),sep=''))) file.remove(paste(.model.modeldatapath(arfmodel),sp,.model.derivativeFile(arfmodel),sep=''))
	
	.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),'Simple Gaussmodel was fitted.')
	
	#set analyticalgrad options
	if(.options.min.analyticalgrad(options)) {
		#if(.options.min.routine(options)[2]=='vpv') gradfunc=gradient.simple else gradfunc=gradient.simple.rpr
		gradfunc=gradient.simple
		angrad=FALSE
	} else {
		gradfunc=NULL
		angrad=FALSE
	}
	
	
	#set which method of model calculation to use 
	#if(.options.min.routine(options)[2]=='vpv') sumsofsquares = ssq.simple else sumsofsquares = ssq.simple.rpr
	sumsofsquares = ssq.simple
	
	#set boundaries in L-BFGS-B mode
	if(length(.options.opt.lower(options))==1 | length(.options.opt.upper(options))==1) {
		lowbound=-Inf
		upbound=Inf
	} else {
		#set location to maximal dim
		max_loc = c(.fmri.data.dims(dat)[2],.fmri.data.dims(dat)[3],.fmri.data.dims(dat)[4])
		
		#set width parameters to maxdim divided by tphe value given in the options
		max_width =  c(.fmri.data.dims(dat)[2],.fmri.data.dims(dat)[3],.fmri.data.dims(dat)[4]) / c(.options.opt.upper(options)[4],.options.opt.upper(options)[5],.options.opt.upper(options)[6])
		
		upbound = rep(c(max_loc,max(max_width),.options.opt.upper(options)[10]),.model.regions(arfmodel))
		lowbound = rep(.options.opt.lower(options)[-(5:9)],.model.regions(arfmodel))
	}	
	
	#check startingvalues
	arfmodel <- validStart(arfmodel)
	
	#final check before fit
	if(!.model.valid(arfmodel)) {
		saveModel(arfmodel)	
		return(arfmodel)
	}
	
	#make progressWindow
	progress = NULL
	
	#runoptim	
	optim.output <- try(suppressWarnings(optim(
							.model.startval(arfmodel),
							sumsofsquares,
							gradfunc,
							lower=lowbound,
							upper=upbound,
							datavec=.fmri.data.datavec(dat)[1:(.fmri.data.dims(dat)[2]*.fmri.data.dims(dat)[3]*.fmri.data.dims(dat)[4])],
							weightvec=.fmri.data.datavec(weights)[1:(.fmri.data.dims(weights)[2]*.fmri.data.dims(weights)[3]*.fmri.data.dims(dat)[4])],
							brain=.model.mask(arfmodel),
							np=.model.regions(arfmodel)*.model.params(arfmodel),
							dimx=.fmri.data.dims(dat)[2],
							dimy=.fmri.data.dims(dat)[3],
							dimz=.fmri.data.dims(dat)[4],
							ss_data=.model.ss(arfmodel),
							analyticalgrad=angrad,
							method=.options.opt.method(options),
							progress=progress,
							control=list(trace=printlevel,maxit=.options.min.iterlim(options)),
							hessian=F
					)),silent=try.silen)
	
	#end_time
	en_time <- Sys.time()
	
	# check for internal errors and set relevant arf model values
	if(is.null(attr(optim.output,'class'))) {
		if(optim.output$convergence==0) .model.convergence(arfmodel) <- paste('Optim converged in ',optim.output$counts[1],' iterations.',sep='')
		if(optim.output$convergence==1) .model.convergence(arfmodel) <- 'Iteration limit exceeded. No convergence.';.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),paste('[min] optim did not converge.',sep=''))
		if(optim.output$convergence==10) .model.convergence(arfmodel) <- 'Degeneracy of the Nelder-Mead Simplex'
		if(optim.output$convergence==51) .model.convergence(arfmodel) <- paste('BFGS raises warning:',gsub('\n','',optim.output$message),sep='')
		if(optim.output$convergence==52) .model.convergence(arfmodel) <-  paste('BFGS raises error:',gsub('\n','',optim.output$message),sep='')
		
		if(optim.output$convergence == 0) .model.valid(arfmodel) <- TRUE else .model.valid(arfmodel) <- FALSE
		
		#set model essentials
		.model.estimates(arfmodel) <- optim.output$par
		.model.iterates(arfmodel) <- optim.output$counts[1]
		.model.proctime(arfmodel)[1,1] <- as.numeric(difftime(en_time,st_time,units='sec'))
			
		if(.model.valid(arfmodel)) {
			#save the ModelBinary
			arfmodel <- saveModelBinSimple(arfmodel)
			
			#set model objects
			.model.minimum(arfmodel) <- optim.output$value
			.model.estimates(arfmodel) <- rep(0,.model.regions(arfmodel)*10)
			
			for(i in 1:.model.regions(arfmodel)) {
				.model.estimates(arfmodel)[1+(10*(i-1))] <- optim.output$par[1+(5*(i-1))]
				.model.estimates(arfmodel)[2+(10*(i-1))] <- optim.output$par[2+(5*(i-1))]
				.model.estimates(arfmodel)[3+(10*(i-1))] <- optim.output$par[3+(5*(i-1))]
				.model.estimates(arfmodel)[4+(10*(i-1))] <- optim.output$par[4+(5*(i-1))]
				.model.estimates(arfmodel)[5+(10*(i-1))] <- optim.output$par[4+(5*(i-1))]
				.model.estimates(arfmodel)[6+(10*(i-1))] <- optim.output$par[4+(5*(i-1))]
				.model.estimates(arfmodel)[7+(10*(i-1))] <- 0
				.model.estimates(arfmodel)[8+(10*(i-1))] <- 0
				.model.estimates(arfmodel)[9+(10*(i-1))] <- 0
				.model.estimates(arfmodel)[10+(10*(i-1))] <- optim.output$par[5+(5*(i-1))]
			}
					
			if(.options.min.analyticalgrad(options)) .model.gradient(arfmodel) <- gradient.simple(.model.estimates(arfmodel),.fmri.data.datavec(dat)[1:(.fmri.data.dims(dat)[2]*.fmri.data.dims(dat)[3]*.fmri.data.dims(dat)[4])],.fmri.data.datavec(weights)[1:(.fmri.data.dims(weights)[2]*.fmri.data.dims(weights)[3]*.fmri.data.dims(dat)[4])],.model.mask(arfmodel),.model.regions(arfmodel)*.model.params(arfmodel),.fmri.data.dims(dat)[2],.fmri.data.dims(dat)[3],.fmri.data.dims(dat)[4],.model.ss(arfmodel),analyticalgrad=T,progress=progress)
		
			if(.model.valid(arfmodel)) {
				#caluclate fits
				arfmodel = BIC(arfmodel,options=options)
				arfmodel = RMSEA(arfmodel,options=options)
			}
		} else .model.warnings(arfmodel) <- c(.model.warnings(arfmodel),paste('[min] optim did not converge.',sep=''))
		
	} else {
		.model.convergence(arfmodel) <- 'Internal error, no convergence.'
		.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),paste('[min] optim internal error: ',gsub('\n','',optim.output),sep=''))
		.model.proctime(arfmodel)[1,1] <- as.numeric(difftime(en_time,st_time,units='sec'))
		.model.valid(arfmodel) <- FALSE
	}
	
	if(!.model.valid(arfmodel)) .model.warnings(arfmodel) <- c(.model.warnings(arfmodel),.model.convergence(arfmodel)) 

	#save the modelInfo
	saveModel(arfmodel)
	
	#return arf model object	
	return(invisible(arfmodel))
}

pruneModel <- 
#prune a given model until neat
function(arfmodel,modelname='defaultmodel',subject='',condition='',grad=NULL,bound=NULL,pval=NULL,options=new('options'),overwrite=T,experiment=NULL)
{
	if(is.null(experiment)) {
		experiment <- try(get('.experiment',envir=.arfInternal),silent=T)
		if(attr(experiment,'class')=='try-error') stop('Experiment not loaded. Run loadExp first.')
	}
	
	stop_prune = FALSE
	prune_num = 1
	
	pruned_model = arfmodel
	
	cat(as.character(Sys.time()),'pruning model with',.model.regions(arfmodel),'regions\n')
	
	while(!stop_prune) 
	{
		ests = matrix(.model.estimates(pruned_model),10)
	
		#check validness of model
		if(.model.valid(pruned_model)) stop_prune = TRUE
		
		if(!is.null(bound)) b_del = checkSolutionReturn(pruned_model,thres=bound) else b_del=numeric(0)
		if(!is.null(grad)) g_del = checkGradientReturn(pruned_model,absthres=grad) else g_del=numeric(0)
		if(!is.null(pval)) ns_del = checkNonSigReturn(pruned_model,alpha=pval) else ns_del=numeric(0)
		
		del = unique(c(b_del,g_del,ns_del))
		
		if(length(del)>0) {
			if(length(del)!=ncol(ests)) {
				ests = ests[,-del]
				ests = as.vector(ests)
				.options.start.method(options) = 'use'
				pruned_model = newModel(paste('pruned',prune_num,'_',modelname,sep=''),regions=length(ests)/10,subject=subject,condition=condition,type='gauss',options=options,overwrite=overwrite,experiment=experiment)
			
				.model.startval(pruned_model) = ests 
				saveModel(pruned_model)
				cat(as.character(Sys.time()),'*fitting pruned model with',.model.regions(pruned_model),'regions\n')
				pruned_model = fitModel(pruned_model,options=options)
				pruned_model <- checkSolution(pruned_model,options,thres=bound)
				
				
			} else {
				model = pruned_model
				cat(as.character(Sys.time()),'Pruned all regions (no regions left)\n')
				stop_prune = TRUE
			}
			
		} else {
			model = pruned_model 
			stop_prune = TRUE
		}
		
		prune_num = prune_num + 1
		
	} #end prune while
	
	return(pruned_model)
	
}

