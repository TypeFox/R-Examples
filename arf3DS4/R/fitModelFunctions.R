#############################################
# arf3DS4 S4 FITMODEL FUNCTIONS				#
# Wouter D. Weeda							#
# University of Amsterdam					#
#############################################

#[CONTAINS]
#ssq.gauss
#model.gauss
#gradient.gauss
#ssq.simple
#model.simple
#gradient.simple
#createAverages					[user]
#createAllAverages				[user]
#determineStartRect				[user]
#determineStartRectSimple		[user]
#regressAmplitudes
#isEdge
#fallOff
#fwhm.filter
#setMask
#validStart
#checkBound
#checkSolution
#checkSolutionReturn
#checkGradientReturn
#checkNonSigReturn
#setIsoContour					[user]
#persistentBound

ssq.gauss <- 
function(theta,datavec,weightvec,brain,np,dimx,dimy,dimz,ss_data,analyticalgrad,progress) 
##ssq.gauss returns the ssq of the full gauss model with an anlytical gradient attached
{
	
	if(length(theta[is.na(theta) | is.nan(theta) | theta==Inf | theta==-Inf])==0)  {
		ssqdat <- .C('ssqgauss',as.double(theta),as.double(datavec),as.double(weightvec),as.integer(brain),as.integer(np),as.integer(dimx),as.integer(dimy),as.integer(dimz),as.double(vector('numeric',1)))[[9]]
	} else ssqdat=ss_data
	
	if(is.nan(ssqdat) | ssqdat==Inf | is.na(ssqdat) | ssqdat==-Inf) ssqdat=ss_data
		
	assign('.theta_latest',theta,envir=.arfInternal)
	
	#get iteration data
	olgrad = get('.oldobj',envir=.arfInternal)
	gradobj = round(olgrad-ssqdat,0)
	objit = get('.objit',envir=.arfInternal)
	gradval = get('.gradval',envir=.arfInternal)
	bounded = get('.bounded',envir=.arfInternal)
	gradit = get('.gradit',envir=.arfInternal)
	
	#Progress Watcher
	if(!progress$disabled) writeProgress(ssqdat,theta,objit,gradobj,gradval,progress,bounded,gradit)
	
	#boundary checker
	boundvec = persistentBound(theta,progress$lower,progress$upper,progress$perslim)
	if(length(boundvec)>0) {
		ssqdat='killoptim'
		assign('.arf_error',list(errtype='persbound',data=boundvec),envir=.arfInternal)
	}
	
	#iterlimchecker
	if(objit>progress$iterlim) {
		ssqdat = 'killoptim'
		assign('.arf_error',list(errtype='iterlim',data=objit),envir=.arfInternal)
	}
	
	assign('.oldobj',ssqdat,envir=.arfInternal)
	assign('.objit',objit+1,envir=.arfInternal)

	return(invisible(ssqdat))	
	
}

model.gauss <-
function(theta,np,dimx,dimy,dimz)
#returns model estimate for full gaussmodel
{
	
	if(length(theta[is.na(theta) | is.nan(theta) | theta==Inf | theta==-Inf])==0)  {
		model <- .C('gauss',as.double(theta),as.integer(np),as.integer(dimx),as.integer(dimy),as.integer(dimz),as.double(rep(0,dimx*dimy*dimz)))[[6]]
	} else model=NA
	
	return(model)
	
}

gradient.gauss <- 
function(theta,datavec,weightvec,brain,np,dimx,dimy,dimz,ss_data,analyticalgrad,progress) 
#gradient returns the analytical gradient of the ssq to the thetaparameters
{

	if(length(theta[is.na(theta) | is.nan(theta) | theta==Inf | theta==-Inf])==0) {
		model <- .C('gauss',as.double(theta),as.integer(np),as.integer(dimx),as.integer(dimy),as.integer(dimz),as.double(rep(0,dimx*dimy*dimz)))[[6]]
		#model = truncDouble(model)
	} else model=NA
	
	if(length(model[is.na(model) | is.nan(model) | model==Inf | model==-Inf])==0) {
		grad <- .C('dfssq',as.integer(np),as.integer(brain),as.integer(dimx),as.integer(dimy),as.integer(dimz),as.double(theta),as.double(datavec),as.double(model),as.double(weightvec),as.double(vector('numeric',np)))[[10]]
		#grad = truncDouble(grad)
	} else grad=rep(1e+12,np) 
	
	#if((length(grad[is.na(grad) | is.nan(grad) | grad==Inf | grad==-Inf])!=0)) grad=rep(1e+12,np) 
	
	assign('.gradient_latest',grad,envir=.arfInternal)
	
	#progress Watcher (only assign gradients)
	assign('.gradval',grad,envir=.arfInternal)

	gradit = get('.gradit',envir=.arfInternal)
	assign('.gradit',gradit+1,envir=.arfInternal)
	
	return(grad)

}

truncDouble <- function(vec) 
{
	vec[vec > 0 & vec<.Machine$double.eps] = .Machine$double.eps
	vec[vec < 0 & abs(vec)<.Machine$double.neg.eps] = .Machine$double.eps*-1
	return(vec)
}


ssq.simple <- 
function(theta,datavec,weightvec,brain,np,dimx,dimy,dimz,ss_data,analyticalgrad,progress) 
##ssq.simple returns the ssq of the simple gauss model with an anlytical gradient attached
{
	
	if(length(theta[is.na(theta) | is.nan(theta) | theta==Inf | theta==-Inf])==0)  {
		ssqdat <- .C('simplessqgauss',as.double(theta),as.double(datavec),as.double(weightvec),as.integer(brain),as.integer(np),as.integer(dimx),as.integer(dimy),as.integer(dimz),as.double(vector('numeric',1)))[[9]]
	} else ssqdat=ss_data
		
	if(is.nan(ssqdat) | ssqdat==Inf | is.na(ssqdat) | ssqdat==-Inf) ssqdat=ss_data
	
	#progress Watcher
	
	
	return(invisible(ssqdat))	
	
}

model.simple <-
function(theta,np,dimx,dimy,dimz)
#returns model estimate for simple gaussmodel
{
	if(length(theta[is.na(theta) | is.nan(theta) | theta==Inf | theta==-Inf])==0)  {
		model <- .C('simplegauss',as.double(theta),as.integer(np),as.integer(dimx),as.integer(dimy),as.integer(dimz),as.double(rep(0,dimx*dimy*dimz)))[[6]]
	} else model=NA
	
	return(model)
	
}

gradient.simple <- 
function(theta,datavec,weightvec,brain,np,dimx,dimy,dimz,ss_data,analyticalgrad,progress) 
#gradient returns the analytical gradient of the ssq to the thetaparameters
{
	if(length(theta[is.na(theta) | is.nan(theta) | theta==Inf | theta==-Inf])==0) {
		model <- .C('simplegauss',as.double(theta),as.integer(np),as.integer(dimx),as.integer(dimy),as.integer(dimz),as.double(rep(0,dimx*dimy*dimz)))[[6]]
	} else model=NA
	
	if(length(model[is.na(model) | is.nan(model) | model==Inf | model==-Inf])==0) {
		grad <- .C('dfsimplessq',as.integer(np),as.integer(dimx),as.integer(dimy),as.integer(dimz),as.double(theta),as.double(datavec),as.double(model),as.double(weightvec),as.double(vector('numeric',np)))[[9]]
	} else grad=rep(1e+12,np) 

	#Progress Watcher
	
	
	return(grad)
	
}


createAverages <- 
function(arfdat,experiment=NULL) 
# createAverages averages of the data and weightfiles, set n mask and ss_data 
{
	
	#check experiment
	#if(is.null(experiment)) if(exists('.experiment')) experiment = .experiment else stop('Experiment not loaded. Run loadExp first.')
	if(is.null(experiment)) {
		experiment <- try(get('.experiment',envir=.arfInternal),silent=T)
		if(attr(experiment,'class')=='try-error') stop('Experiment not loaded. Run loadExp first.')
	}
	
	#set filesep 
	sp=.Platform$file.sep
	
	# add run data to avgdat and weightdat
	avgdat <- avgweight <- 0
	for(i in 1:.data.runs(arfdat)) {
		avgdat <- avgdat + .fmri.data.datavec(readData(.data.betafiles(arfdat)[i]))
		avgweight <- avgweight + .fmri.data.datavec(readData(.data.weightfiles(arfdat)[i]))
	}
	
	# divide avgdat by runnumber and avgweight by runnumber^2
	avgdat <- avgdat / .data.runs(arfdat)
	avgweight <- avgweight / .data.runs(arfdat)^2
	
	# get header info of first file (or reference file id supplied)
	headinf <- readHeader(getFileInfo(.data.betafiles(arfdat)[1])) 
	
	# write header and datainfo for datafiles
	filename <- paste(.data.fullpath(arfdat),sp,.experiment.dataDir(experiment),sp,.experiment.avgDir(experiment),sp,.experiment.avgdatFile(experiment),'.',.nifti.header.extension(headinf),sep='')
	
	if(.nifti.header.gzipped(headinf)==T) filename <- paste(filename,'.gz',sep='') 
	headinf <- newFile(filename,headinf)
	.nifti.header.descrip(headinf) <- 'Average data image (ARF)'
	.data.avgdatfile(arfdat) <- filename
	writeData(headinf,avgdat)
	
	# get header info of first file (or reference file id supplied)
	headinf <- readHeader(getFileInfo(.data.weightfiles(arfdat)[1])) 
	
	# write header and datainfo for weightfiles
	filename <- paste(.data.fullpath(arfdat),sp,.experiment.dataDir(experiment),sp,.experiment.avgDir(experiment),sp,.experiment.avgWFile(experiment),'.',.nifti.header.extension(headinf),sep='')
	
	if(.nifti.header.gzipped(headinf)==T) filename <- paste(filename,'.gz',sep='') 
	headinf <- newFile(filename,headinf)
	.nifti.header.descrip(headinf) <- 'Average weights image (ARF)'
	.data.avgWfile(arfdat) <- filename	
	writeData(headinf,avgweight)
	
	#save t-image
	filename <- paste(.data.fullpath(arfdat),sp,.experiment.dataDir(experiment),sp,.experiment.avgDir(experiment),sp,.experiment.avgtstatFile(experiment),'.',.nifti.header.extension(headinf),sep='')
	
	if(.nifti.header.gzipped(headinf)==T) filename <- paste(filename,'.gz',sep='') 
	headinf <- newFile(filename,headinf)
	.nifti.header.descrip(headinf) <- 'Average t-statistics image (ARF)'
	.data.avgtstatFile(arfdat) <- filename
	
	avgtstat <- avgdat/sqrt(avgweight)
	avgtstat[is.nan(avgtstat)]=0
	
	writeData(headinf,avgtstat)
	
	#define mask
	brain <- rep(1,length(avgtstat))
	nb <- which(avgtstat==0)
	if(length(nb)>0) brain[nb]=0
	.data.mask(arfdat) <- brain
	
	#define n
	.data.n(arfdat) <- sum(brain) 
	
	#get ssq_data
	.data.ss(arfdat) <- .C('ssqdata',as.double(avgdat),as.double(avgweight),as.integer(brain),as.integer(length(avgtstat)),as.double(numeric(1)))[[5]]
	
	#save arfdatafile with updated weights
	save(arfdat,file=paste(.data.fullpath(arfdat),sp,.experiment.dataDir(experiment),sp,.experiment.dataRda(experiment),sep=''))
	
	# return data class object	
	return(invisible(arfdat))
	
}

createAllAverages <- 
function(experiment=NULL) 
#create All averages is a wrapper to create all averages in an experiment
{
	
	#check experiment
	if(is.null(experiment)) {
		experiment <- try(get('.experiment',envir=.arfInternal),silent=T)
		if(attr(experiment,'class')=='try-error') stop('Experiment not loaded. Run loadExp first.')
	}
	
	#set filesep
	sp <- .Platform$file.sep
		
	#set subjects and conditions
	subs <- .experiment.subject.names(experiment)
	conds <- .experiment.condition.names(experiment)
	
	#run through all subs and conds
	for(subject in subs) {
		for(condition in conds) {
			
			#make filename
			filename <- paste(.experiment.path(experiment),sp,.experiment.subjectDir(experiment),sp,subject,sp,.experiment.conditionDir(experiment),sp,condition,sp,.experiment.dataDir(experiment),sp,.experiment.dataRda(experiment),sep='')
			
			#check and create, else warning
			if(file.exists(filename)) {
				createAverages(loadRda(filename))
			} else {
				warning(paste(filename,'does not exist. Avg not created'))
			}
		}
		
	}
	
	return(invisible(TRUE))
}



determineStartRect <- 
function(arfmodel,options=loadOptions(arfmodel)) 
# determineStartRect calculates starting values for regions (rectangular mode)
{
	
	.model.modeltype(arfmodel) <- 'gauss'
	.model.params(arfmodel) <- 10
	
	#load in fmriData
	fmridata <- readData(.model.avgtstatFile(arfmodel))
	
	#set theta to the default values (for all regions)
	theta <- rep(.options.start.vector(options),.model.regions(arfmodel))
	
	#set dimensions and read in data
	dimx <- .fmri.data.dims(fmridata)[2]
	dimy <- .fmri.data.dims(fmridata)[3]
	dimz <- .fmri.data.dims(fmridata)[4]
	data <- .fmri.data.datavec(fmridata)[1:(dimx*dimy*dimz)]
	
	mindim=c(1,1,1)
	maxdim=c(dimx,dimy,dimz)
	min_amp = .options.start.vector(options)[10]*-1
	max_amp = .options.start.vector(options)[10]
	
	#set dims of the data
	dim(data) <- c(dimx,dimy,dimz)
	
	#set location matrix
	location <- data.frame(x=rep(seq(1:dimx),times=dimz*dimy),y=rep(rep(seq(1:dimy),each=dimx),times=dimz),z=rep(seq(1:dimz),each=dimx*dimy))
	
	#set mask
	mask <- .model.mask(arfmodel)
	dim(mask) <-  c(dimx,dimy,dimz)
	
	for(reg in 1:.model.regions(arfmodel)) {
		
		#retrieve location in x,y,z
		m <- location[which.max(abs(data)),]
		if(data[m$x,m$y,m$z]<0) theta[10+(10*(reg-1))]=min_amp else theta[10+(10*(reg-1))]=max_amp
		
		#set maximum locations
		theta[1+(10*(reg-1))] <- m$x
		theta[2+(10*(reg-1))] <- m$y
		theta[3+(10*(reg-1))] <- m$z
				
		#caluclatefalloff
		xf <- fallOff(data[,m$y,m$z],.options.start.maxfac(options))
		yf <- fallOff(data[m$x,,m$z],.options.start.maxfac(options))
		zf <- fallOff(data[m$x,m$y,],.options.start.maxfac(options))
		
		#set width in x,y and z dirs
		theta[4+(10*(reg-1))] <- round(mean(xf))
		theta[5+(10*(reg-1))] <- round(mean(yf))
		theta[6+(10*(reg-1))] <- round(mean(zf))
		
		#check widths for zeroes
		theta[4+(10*(reg-1))][theta[4+(10*(reg-1))]<=0]=1 
		theta[5+(10*(reg-1))][theta[5+(10*(reg-1))]<=0]=1 
		theta[6+(10*(reg-1))][theta[6+(10*(reg-1))]<=0]=1 
		
		#zero tha data
		xvec <- (m$x-xf[1]):(m$x+xf[2])
		yvec <- (m$y-yf[1]):(m$y+yf[2])
		zvec <- (m$z-zf[1]):(m$z+zf[2])
		
		rmx <- which(xvec<mindim[1] | xvec>maxdim[1])
		rmy <- which(yvec<mindim[2] | yvec>maxdim[2])
		rmz <- which(zvec<mindim[3] | zvec>maxdim[3])
		
		if(length(rmx)>0 & length(rmx)<length(xvec)) xvec <- xvec[-rmx]
		if(length(rmy)>0 & length(rmy)<length(yvec)) yvec <- yvec[-rmy]
		if(length(rmz)>0 & length(rmz)<length(zvec)) zvec <- zvec[-rmz]	
		
		data[xvec,yvec,zvec]=0

	}
	
	#save startvalues to arfmodel
	.model.startval(arfmodel) <- theta
	
	#regress amplitudes
	theta[((1:.model.regions(arfmodel))*10)]=regressAmplitudes(arfmodel,'start')
	
	#save startingvalues
	.model.startval(arfmodel) <- theta	
	saveStart(.model.startval(arfmodel),arfmodel)
	
	#save startmap
	.fmri.data.fullpath(fmridata) <- .model.modeldatapath(arfmodel)
	.fmri.data.filename(fmridata) <- 'startmap'
	.fmri.data.intent_name(fmridata) <- 'start_search_map'
	writeData(fmridata,as.vector(data))
	
	#save model
	saveModel(arfmodel)
	
	return(invisible(arfmodel))
	
}

determineStartRectSimple <- 
function(arfmodel,options=loadOptions(arfmodel)) 
# determineStartRect calculates starting values for regions (rectangular mode)
{
	
	.model.modeltype(arfmodel) <- 'simple'
	.model.params(arfmodel) <- 5
	
	#load in fmriData
	fmridata <- readData(.model.avgtstatFile(arfmodel))
	
	#set theta to the default values (for all regions)
	theta <- rep(.options.start.vector(options),.model.regions(arfmodel))
	newstart <- .model.regions(arfmodel)*5
	
	#set dimensions and read in data
	dimx <- .fmri.data.dims(fmridata)[2]
	dimy <- .fmri.data.dims(fmridata)[3]
	dimz <- .fmri.data.dims(fmridata)[4]
	data <- .fmri.data.datavec(fmridata)[1:(dimx*dimy*dimz)]
	
	mindim=c(1,1,1)
	maxdim=c(dimx,dimy,dimz)
	min_amp = .options.start.vector(options)[10]*-1
	max_amp = .options.start.vector(options)[10]
	
	#set dims of the data
	dim(data) <- c(dimx,dimy,dimz)
	
	#set location matrix
	location <- data.frame(x=rep(seq(1:dimx),times=dimz*dimy),y=rep(rep(seq(1:dimy),each=dimx),times=dimz),z=rep(seq(1:dimz),each=dimx*dimy))
	
	mask = .model.mask(arfmodel)
	dim(mask)=c(dimx,dimy,dimz)
	
	for(reg in 1:.model.regions(arfmodel)) {
		
		#retrieve location in x,y,z
		m <- location[which.max(abs(data)),]
		if(data[m$x,m$y,m$z]<0) theta[10+(10*(reg-1))]=min_amp else theta[10+(10*(reg-1))]=max_amp
				
		#set maximum locations
		theta[1+(10*(reg-1))] <- m$x
		theta[2+(10*(reg-1))] <- m$y
		theta[3+(10*(reg-1))] <- m$z
		
		#caluclatefalloff
		xf <- fallOff(data[,m$y,m$z],.options.start.maxfac(options))
		yf <- fallOff(data[m$x,,m$z],.options.start.maxfac(options))
		zf <- fallOff(data[m$x,m$y,],.options.start.maxfac(options))
		
		#set width in x,y and z dirs
		theta[4+(10*(reg-1))] <- round(mean(xf))
		theta[5+(10*(reg-1))] <- round(mean(yf))
		theta[6+(10*(reg-1))] <- round(mean(zf))
		
		#check widths for zeroes
		theta[4+(10*(reg-1))][theta[4+(10*(reg-1))]<=0]=1 
		theta[5+(10*(reg-1))][theta[5+(10*(reg-1))]<=0]=1 
		theta[6+(10*(reg-1))][theta[6+(10*(reg-1))]<=0]=1 
				
		#zero tha data
		xvec <- (m$x-xf[1]):(m$x+xf[2])
		yvec <- (m$y-yf[1]):(m$y+yf[2])
		zvec <- (m$z-zf[1]):(m$z+zf[2])
		
		rmx <- which(xvec<mindim[1] | xvec>maxdim[1])
		rmy <- which(yvec<mindim[2] | yvec>maxdim[2])
		rmz <- which(zvec<mindim[3] | zvec>maxdim[3])
		
		if(length(rmx)>0 & length(rmx)<length(xvec)) xvec <- xvec[-rmx]
		if(length(rmy)>0 & length(rmy)<length(yvec)) yvec <- yvec[-rmy]
		if(length(rmz)>0 & length(rmz)<length(zvec)) zvec <- zvec[-rmz]	
				
		data[xvec,yvec,zvec]=0
				
		newstart[1+(5*(reg-1))]=theta[1+(10*(reg-1))]
		newstart[2+(5*(reg-1))]=theta[2+(10*(reg-1))]
		newstart[3+(5*(reg-1))]=theta[3+(10*(reg-1))]
		newstart[4+(5*(reg-1))]=mean(theta[4+(10*(reg-1))],theta[5+(10*(reg-1))],theta[6+(10*(reg-1))])
		newstart[5+(5*(reg-1))]=theta[10+(10*(reg-1))]
			
	}
	
	#regressamps
	.model.startval(arfmodel) <- theta
	theta[((1:.model.regions(arfmodel))*10)]=regressAmplitudes(arfmodel,'start')
	newstart[5+(5*((1:.model.regions(arfmodel))-1))]=theta[((1:.model.regions(arfmodel))*10)]
		
	#save startingvalues
	.model.startval(arfmodel) <- newstart
	saveStart(.model.startval(arfmodel),arfmodel)
	
	#save startmap
	.fmri.data.fullpath(fmridata) <- .model.modeldatapath(arfmodel)
	.fmri.data.filename(fmridata) <- 'startmap'
	.fmri.data.intent_name(fmridata) <- 'start_search_map'
	writeData(fmridata,as.vector(data))
	
	#save model
	saveModel(arfmodel)
	
	return(invisible(arfmodel))
	
}

regressAmplitudes <-
function(arfmodel,which='start')
{
	
	which = match.arg(which,c('start','estimates'))
	
	#get Header info from avgdatfile
	headinf <- readHeader(getFileInfo(.model.avgdatfile(arfmodel)))
	n = .nifti.header.dims(headinf)[2]*.nifti.header.dims(headinf)[3]*.nifti.header.dims(headinf)[4]
	regs = 1:.model.regions(arfmodel)
	
	#make model design matrix
	X = matrix(NA,n,length(regs))
	
	if(.model.modeltype(arfmodel)=='simple') .model.params(arfmodel)=10
	if(which=='estimates') theta = matrix(.model.estimates(arfmodel),.model.params(arfmodel))
	if(which=='start') theta = matrix(.model.startval(arfmodel),.model.params(arfmodel))
		
	#set theta amplitudes to one (all)
	theta[10,]=rep(1,.model.regions(arfmodel))
	
	p=1
	for(i in regs) {
		thetavec = as.vector(theta[,i])
		X[,p] = .C('gauss',as.double(thetavec),as.integer(.model.params(arfmodel)),as.integer(.nifti.header.dims(headinf)[2]),as.integer(.nifti.header.dims(headinf)[3]),as.integer(.nifti.header.dims(headinf)[4]),as.double(numeric(.nifti.header.dims(headinf)[2]*.nifti.header.dims(headinf)[3]*.nifti.header.dims(headinf)[4])))[[6]]
		p=p+1
	}
	
	
	funcdata = readData(.model.avgdatfile(arfmodel))	
	funcvolume = .fmri.data.datavec(funcdata)
	dim(funcvolume) = c(.fmri.data.dims(funcdata)[2],.fmri.data.dims(funcdata)[3],.fmri.data.dims(funcdata)[4])
	
	
	Xp = solve(t(X)%*%X)%*%t(X)
	y = as.vector(funcvolume)
	b = Xp%*%y
	
	return(b)
	
}


isEdge <-  
function(mask,xvec,yvec,zvec)
#check if a cube-region touches non-brain voxels
{
	if(sum(as.vector(mask[xvec,yvec,zvec]))==length(as.vector(mask[xvec,yvec,zvec]))) return(FALSE) else return(TRUE)
		
}

fallOff <- 
function(vec,fwhm=2)
#calculates mean falloff of a vector
{
	#smooth vec with fwhm filter
	vec <- fwhm.filter(vec,fwhm)
	
	#determine max of vector and falloff amount
	m <- which.max(vec)
	maxval <- max(vec)
	falloffval <- maxval/2
	
	#set min and max-dim elements
	maxdim <- length(vec)
	mindim <- 1
	
	#check falloff to the right
	i=0
	while(vec[m+i]>falloffval) {
		i=i+1;		
		if((m+i)>=maxdim) break()
	}
	if(i>1) i=i-1
	
	#check falloff to the left
	j=0
	while(vec[m-j]>falloffval) {
		j=j+1;		
		if((m-j)<=mindim) break()
	}
	if(j>1) j=j-1
	
	#return left and right values
	return(c(j,i))
}

fwhm.filter <- 
function(vec,fwhm) 
#smooth a datavector (1D) using FWHM filter (used in startval detect)
{
	
	
	fl=50
	filt=dnorm(1:100,mean=50,sd=fwhm)
	
	fdat <- convolve(vec,filt,type='open')
	vec <- fdat[-c(1:(fl-1),(length(fdat)-(fl-1)):length(fdat))]
	
	return(vec)
	
}

setMask <- 
function(arfmodel) 
#set mask attributes if necessary
{
	
	avgtstat <- .fmri.data.datavec(readData(.model.avgtstatFile(arfmodel)))
	
	#define mask
	brain <- rep(1,length(avgtstat))
	nb <- which(avgtstat==0)
	if(length(nb)>0) brain[nb]=0
	.model.mask(arfmodel) <- brain
	
	#define n
	.model.n(arfmodel) <- sum(brain)
	
	#get ssq_data
	.model.ss(arfmodel) <- .C('ssqdata',as.double(.fmri.data.datavec(readData(.model.avgdatfile(arfmodel)))),as.double(.fmri.data.datavec(readData(.model.avgWfile(arfmodel)))),as.integer(brain),as.integer(length(avgtstat)),as.double(numeric(1)))[[5]]
	
	return(invisible(arfmodel))
}

validStart <- 
function(arfmodel) 
#check if startvalues are valid
{
	
	theta <- .model.startval(arfmodel)
	mess = character(0)
	onlywarn=TRUE
		
	dat=readData(.model.avgtstatFile(arfmodel))
	mask = .model.mask(arfmodel)
	dim(mask)=c(.fmri.data.dims(dat)[2],.fmri.data.dims(dat)[3],.fmri.data.dims(dat)[4])
	
	if(length(.model.startval(arfmodel))!=(.model.regions(arfmodel)*.model.params(arfmodel))) {
		mess=c(mess,'[startval] Vector of startingvalues are not a multiple of regions*parameters')
		onlywarn=FALSE
	} else {
		for(reg in 1:.model.regions(arfmodel)) {
			st_loc = theta[(1:3)+(.model.params(arfmodel)*(reg-1))]
			if(.model.modeltype(arfmodel)=='gauss') st_sd = theta[(4:6)+(.model.params(arfmodel)*(reg-1))]  else st_sd=rep(theta[(4)+(.model.params(arfmodel)*(reg-1))],3)
			if(.model.modeltype(arfmodel)=='gauss') st_rho = theta[(7:9)+(.model.params(arfmodel)*(reg-1))] else st_rho=rep(0,3)
			if(.model.modeltype(arfmodel)=='gauss') st_amp = theta[(10)+(.model.params(arfmodel)*(reg-1))] else st_amp=theta[(5)+(.model.params(arfmodel)*(reg-1))]
			
			sigma = matrix(c(st_sd[1]^2,st_sd[1]*st_sd[2]*st_rho[1],st_sd[1]*st_sd[3]*st_rho[2],st_sd[1]*st_sd[2]*st_rho[1],st_sd[2]^2,st_sd[2]*st_sd[3]*st_rho[3],st_sd[1]*st_sd[3]*st_rho[2],st_sd[2]*st_sd[3]*st_rho[3],st_sd[3]^2),3,3)
			det_sig=det(sigma)
			
			if(!is.na(det_sig) & !is.nan(det_sig) & det_sig!=Inf & det_sig!=-Inf) {
				#if(det_sig<1 & det_sig>0) mess = c(mess,paste('[startval] Region ',reg,' has a small volume (',round(det_sig,2),').',sep=''))
				if(det_sig<=0) mess = c(mess,paste('[startval] Region ',reg,' has a zero volume (',round(det_sig,2),').',sep=''))
				
			} else {
				mess = c(mess,paste('[startval] Region ',reg,' returns NA/NaN/Inf/-Inf as determinant.',sep=''))
				onlywarn=FALSE
			}
			
			if(st_loc[1]<1 | st_loc[2]<1 | st_loc[3]<1) {
				mess=c(mess,paste('[startval] Region ',reg,' has a location smaller than 1',sep=''))
			} 
			
			if(st_amp==0) {
				mess=c(mess,paste('[startval] Region ',reg,' has an amplitude of zero.'))
				onlywarn=FALSE
			}
		}
		
	}
	
	if(length(mess)>0) {
		.model.warnings(arfmodel) <- c(.model.warnings(arfmodel),mess)
		if(!onlywarn) .model.valid(arfmodel) <- FALSE else .model.valid(arfmodel) <- TRUE
	} else {
		.model.valid(arfmodel) <- TRUE
	}
	
	if(!.model.valid(arfmodel)) .model.convergence(arfmodel) <- 'Starting values are not valid. Minimization routine not started.'
	
	return(invisible(arfmodel))
	
}


checkBound <- 
function(arfmodel,lowbound,upbound,thres=6) 
#check if parameters are on the bound
{
	
	estimates = matrix(.model.estimates(arfmodel),.model.params(arfmodel))

	for(i in 1:.model.params(arfmodel)) {
		regs = which(round(estimates[i,],thres)<=round(lowbound[i],thres) | round(estimates[i,],thres)>=round(upbound[i],thres))
		
		if(length(regs)>0) {
			mess = paste('[optim] Parameter',i,'is at boundary for region(s)',paste(regs,collapse=","))
			.model.warnings(arfmodel) = c(.model.warnings(arfmodel),mess)
		}
		
	}
	
	return(arfmodel)
	
}


checkSolution <- 
function(arfmodel,options=loadOptions(arfmodel),dat=readData(.model.avgdatfile(arfmodel)),thres=6) 
#check the solution for boundaries
{
	
	#set boundaries in L-BFGS-B mode
	if(length(.options.opt.lower(options))==1 | length(.options.opt.upper(options))==1) {
		lowbound=-Inf
		upbound=Inf
	} else {
		#set location to maximal dim
		max_loc = c(.fmri.data.dims(dat)[2],.fmri.data.dims(dat)[3],.fmri.data.dims(dat)[4])
		
		#set width parameters to maxdim divided by tphe value given in the options
		max_width =  c(.fmri.data.dims(dat)[2],.fmri.data.dims(dat)[3],.fmri.data.dims(dat)[4]) / c(.options.opt.upper(options)[4],.options.opt.upper(options)[5],.options.opt.upper(options)[6])
		
		upbound = rep(c(max_loc,max_width,.options.opt.upper(options)[7:10]),.model.regions(arfmodel))
		lowbound = rep(.options.opt.lower(options),.model.regions(arfmodel))
	}
	
	arfmodel = checkBound(arfmodel,lowbound,upbound,thres=thres)
	regstodel = checkSolutionReturn(arfmodel,options,dat,thres)
	
	if(length(regstodel)>0) .model.valid(arfmodel) = FALSE
	
	return(arfmodel)

}


checkSolutionReturn <- 
function(arfmodel,options=loadOptions(arfmodel),dat=readData(.model.avgdatfile(arfmodel)),thres=8) 
#check the solution for boundaries
{
	
	#set boundaries in L-BFGS-B mode
	if(length(.options.opt.lower(options))==1 | length(.options.opt.upper(options))==1) {
		lowbound=-Inf
		upbound=Inf
	} else {
		#set location to maximal dim
		max_loc = c(.fmri.data.dims(dat)[2],.fmri.data.dims(dat)[3],.fmri.data.dims(dat)[4])
		
		#set width parameters to maxdim divided by tphe value given in the options
		max_width =  c(.fmri.data.dims(dat)[2],.fmri.data.dims(dat)[3],.fmri.data.dims(dat)[4]) / c(.options.opt.upper(options)[4],.options.opt.upper(options)[5],.options.opt.upper(options)[6])
		
		upbound = rep(c(max_loc,max_width,.options.opt.upper(options)[7:10]),.model.regions(arfmodel))
		lowbound = rep(.options.opt.lower(options),.model.regions(arfmodel))
	}
	
	#fill matrix (params, by regs with bounded params)	
	estimates = matrix(.model.estimates(arfmodel),.model.params(arfmodel))
	whichwhat = matrix(0,.model.params(arfmodel),.model.regions(arfmodel))
	
	for(i in 1:.model.params(arfmodel)) {
		regs = which(round(estimates[i,],thres)<=round(lowbound[i],thres) | round(estimates[i,],thres)>=round(upbound[i],thres))
		
		if(length(regs)>0) {
			whichwhat[i,regs] = 1 
		}
		
	}
	
	#check whichwhat
	regstodel = which(apply(whichwhat,2,sum)>0)	

	return(regstodel)	
	
}

checkGradientReturn <- 
function(arfmodel,absthres=1000) 
#check the solution for boundaries
{

	#fill matrix (params, by regs with bounded params)	
	gradient = matrix(.model.gradient(arfmodel),.model.params(arfmodel))
	whichwhat = matrix(0,.model.params(arfmodel),.model.regions(arfmodel))
	
	for(i in 1:.model.params(arfmodel)) {
		regs = which(abs(gradient[i,])>=absthres)
		if(length(regs)>0) {
			whichwhat[i,regs] = 1 
		}
		
	}
	
	#check whichwhat
	regstodel = which(apply(whichwhat,2,sum)>0)	
	
	return(regstodel)	
	
}


checkNonSigReturn <-
function(arfmodel,alpha=.05) 
#check for non_sigs
{
	
	if(length(.model.varcov(arfmodel))==0) {
		arfmodel = varcov(arfmodel)
		arfmodel = wald(arfmodel)
	}
	
	if(.model.valid(arfmodel)) {
		ns_del = which(.wald.pvalues(.model.wald(arfmodel))[,4]>=alpha)
	} else 	ns_del = numeric(0)
		
	return(ns_del)
	
	
}


setIsoContour <- 
function(arfmodel,conf.int=95)
{
	#set confidence interval and Chi-square threshold
	CI = (conf.int/100)
	Chi = qchisq(CI,3)
	
	#read dat, set dims, and create new 0-filled object
	arfdat = readData(.model.avgdatfile(arfmodel))
	dimx = .fmri.data.dims(arfdat)[2]
	dimy = .fmri.data.dims(arfdat)[3]
	dimz = .fmri.data.dims(arfdat)[4]
	
	
	totvec = numeric(0)
	
	#read in estimates in a matrix
	ests = matrix(.model.estimates(arfmodel),.model.params(arfmodel))
	
	#for each region set the confidence interval based on X'C-1X <= Chi3
	for(reg in 1:.model.regions(arfmodel)) {
	
		newdat = array(0,c(dimx,dimy,dimz))
		
		theta = ests[,reg]
		Sigma = matrix(c(theta[4]^2,theta[4]*theta[5]*theta[7],theta[4]*theta[6]*theta[8],theta[4]*theta[5]*theta[7],theta[5]^2,theta[5]*theta[6]*theta[9],theta[4]*theta[6]*theta[8],theta[5]*theta[6]*theta[9],theta[6]^2),3)
		SI = solve(Sigma)
		
		#cat(det(Sigma),'\n')
		
		for(z in 1:dimz) {
			for(y in 1:dimy) {
				for(x in 1:dimx) {
					
					X = c(x-theta[1],y-theta[2],z-theta[3]) 
					C = t(X)%*%SI%*%X
					if(C<=Chi) newdat[x,y,z] = 1
					
				}
			}
		}
		
		totvec = c(totvec,as.vector(newdat))
	}

	#set new data object attributes (path and name)
	.fmri.data.dims(arfdat)[1] = 4
	.fmri.data.dims(arfdat)[5] = .model.regions(arfmodel)
	.fmri.data.fullpath(arfdat) = .model.modeldatapath(arfmodel)
	.fmri.data.filename(arfdat) = paste('iscontours_CI',as.character(conf.int),sep='')
	.fmri.data.descrip(arfdat) = 'ARF regions Isocontours'
	.fmri.data.datavec(arfdat) = as.vector(totvec)
	
	writeData(arfdat)
	
	return(arfdat)
	
}

persistentBound <-
function(theta,lower,upper,itlim)
#check for persistent boundary violation of theta 7,8,9
{
	oldbounds = get('.bounded',envir=.arfInternal)

	estvec = matrix(theta,10)
	bounded = rep(NA,length(oldbounds))
	
	bmat = matrix(NA,10,length(theta)/10)
	for(param in 1:dim(estvec)[1]) {
		bmat[param,] = as.numeric(estvec[param,]<=lower[param] | estvec[param,]>=upper[param])
	}
	
	newbounds = as.numeric(apply(bmat,2,any))
	
	for(i in 1:length(oldbounds)) if(newbounds[i]>0) bounded[i]=oldbounds[i]+newbounds[i] else bounded[i]=0

	assign('.bounded',bounded,envir=.arfInternal)
	
	return(which(bounded>itlim))
}
