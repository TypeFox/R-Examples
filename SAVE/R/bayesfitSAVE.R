##########################################################################
## Bayesfit Method
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 3, April 2013.
##
## Copyright (C) 2013-present Jesus Palomo, Gonzalo Garcia-Donato, 
##							  and Rui Paulo
##    
##########################################################################
bayesfit.SAVE <- function (object, prior,mcmcMultmle=1,
		     prob.prop=0.5, method=2, n.iter, nMH=20, n.burnin=0,
         	 n.thin=1,verbose=FALSE){
    	 	
	# To deprecate the unused parameters and to include in the 
	# call() all the default parameters not used in the call
	# to the function
	dprct <- .deprecate.parameters(call=sys.call(sys.parent(1)))
	object@bayesfitcall <- as.call(dprct)
	
	#####
	#Write the file for the prior
	if (!is.null(object@calibrationnames) && length(object@calibrationnames)!=0){
		if ((length(prior)/6)<length(object@calibrationnames))
          {stop("Not enough prior distributions\n")}
		if ((length(prior)/6)>length(object@calibrationnames))
          {stop("Too many prior distributions\n")}
		prior.matrix<- matrix(0, ncol=length(object@calibrationnames), nrow=5)
		colnames(prior.matrix)<- object@calibrationnames
		for (i.name in object@calibrationnames){
			thisplace<- which(prior==i.name)
			prior.matrix[,i.name]<- as.numeric(prior[thisplace+1:5])
		}
		object@prior <- t(prior.matrix)
		write.table(object@prior, file=paste(object@wd,"/bounds.dat",sep=""),
				col.names=F, row.names=F)
	} else {
		if (length(prior)!=0) print ("The specified priors for the calibration parameters will be deprecated")
		aux <- which(names(object@bayesfitcall)=="prior")
		object@bayesfitcall <- object@bayesfitcall[-aux]
		prior.matrix <- matrix(0,0,0)
		object@prior <- prior.matrix
		#write.table(object@prior, file=paste(object@wd,"/bounds.dat",sep=""),
		#		col.names=F, row.names=F)
	}

    ## Beginning of bayesfit
    
    ####################
	#load some auxiliary functions
	duplicates <- function(dat)
	{
		s <- do.call("order", as.data.frame(dat))
		if(dim(as.matrix(dat))[2]==1){
			non.dup <- !duplicated(as.matrix(dat[s]))
		}
		else{
			non.dup <- !duplicated(as.matrix(dat[s,]))
		}
		orig.ind <- s[non.dup]
		first.occ <- orig.ind[cumsum(non.dup)]
		#first.occ[non.dup] <- NA
		first.occ[order(s)]
	}
	###############end of loading auxiliary functions.

    #####
    #Field data:
    #Select those columns corresponding to the controllable inputs:
    # x<- field.data[,object@controllablenames]
    #Extract the unique part of the field design matrix as required by bayesfit:
	if(!object@constant.controllables){
    	x.unique <- as.data.frame(unique(object@df))
    	names(x.unique) <- object@controllablenames
    	my.original <- duplicates(object@df)
    	yf.ordered <- numeric(0)
    	nreps <- NULL
    	for (i in unique(my.original)){
        	yf.ordered <- c(yf.ordered,object@yf[my.original==i])
        	nreps <- c(nreps, sum(my.original==i))
    	}
	}
	else{
		x.unique<- as.data.frame(1)
		nreps<- length(object@yf)
		yf.ordered<- object@yf
	}
		
    
	#Write to the files
	#file field_data.dat:
	write.table(yf.ordered, file=paste(object@wd,"/field_data.dat", sep=""),
                col.names=F, row.names=F)
    
	#file field_nreps.dat:
	write.table(nreps, file=paste(object@wd,"/field_nreps.dat", sep=""),
                col.names=F, row.names=F)
    
	#file field_unique.dat:
	write.table(x.unique, file=paste(object@wd,"/field_unique.dat",sep=""),
                col.names=F, row.names=F)
	#####
	#####
	#Model data:
	write.table(as.data.frame(object@dm), file=paste(object@wd,"/model_inputs.dat",sep=""),
                col.names=F, row.names=F)
    
	write.table(as.vector(object@ym),
                file=paste(object@wd,"/model_data.dat",sep=""),
                col.names=F, row.names=F)

    #####
    #####
    #Mean response:
    write.table(object@xm,file=paste(object@wd,"/mcmc.field.design.M.matrix.dat",sep=""),
                col.names=F, row.names=F)
    write.table(object@xf,file=paste(object@wd,"/mcmc.field.design.F.matrix.dat",sep=""),
                col.names=F, row.names=F)
    #####
    if (!object@constant.controllables){dthetaM <- length(object@calibrationnames) + length(object@controllablenames)}
    else{dthetaM <- length(object@calibrationnames)}

    dthetaM <- dthetaM*2+1
    dthetaL <- ncol(object@xm)

    if (!object@constant.controllables){dthetaF <- 2*length(object@controllablenames)+2}
    else{dthetaF <- 2}

    if(length(object@mle$thetaM)!=dthetaM)
        {#object@mle<- NULL;
            stop("Dimension of MLE for the covariance of the computer model is incorrect\n")}
    if(length(object@mle$thetaL)!=dthetaL)
        {#object@mle<- NULL;
            stop("Dimension of MLE for the mean of the computer model is incorrect\n")}
    if(length(object@mle$thetaF)!=dthetaF)
        {#object@mle<- NULL;
            stop("Dimension of MLE for the second stage parameters is incorrect\n")}

    #####
    #Write the MLE to the appropriate files
    write.table(object@mle$thetaM, file=paste(object@wd,"/thetaM_mle.dat", sep=""),
                col.names=F, row.names=F)
    write.table(object@mle$thetaL, file=paste(object@wd,"/thetaL_mle.dat", sep=""),
                col.names=F, row.names=F)
    if (!object@constant.controllables){
        write.table(object@mle$thetaF, file=paste(object@wd,"/thetaF_mle.dat", sep=""),
                    col.names=F, row.names=F)}
    else{
        thetaFF<- c(object@mle$thetaF[1],1.0,2.0,object@mle$thetaF[2])
        write.table(thetaFF, file=paste(object@wd,"/thetaF_mle.dat", sep=""),
                    col.names=F, row.names=F)
    }


	#####
	#Values of the parameters:
	#total number of inputs:
    if (!object@constant.controllables){
        numInputs<- length(c(object@controllablenames,object@calibrationnames))
		#dimension of field design unique inputs (NF in Rui''s notation)
		sizeField<- dim(object@xf)[1]
    }
    else{
        numInputs<- length(object@calibrationnames)
		#dimension of field design unique inputs (NF in Rui''s notation)
		sizeField<- 1
    }
    #number of calibration inputs:
	numCalibration<- length(object@calibrationnames)
	#dimension of the linear model for the mean of the GP prior
	#        (q in C notation)
	numPModel<- dim(object@xm)[2]
	#dimension of code design (NM in C notation)

	sizeData<- length(object@ym)
#Probability of sampling from the prior
	object@method <- method
	object@mcmcMultmle <- mcmcMultmle

	# load the C code
	if(!is.loaded("bayesfit")) {
		lib.file <- file.path(paste("bayesfit",
                                            .Platform$dynlib.ext, sep=""))
		dyn.load(lib.file)
		cat(" -Loaded ", lib.file, "\n")
	}

    #Call to the fucnction:
	if ((object@method !=1) && (object@method != 2)){
		stop ("Wrong type of method introduced as parameter")
	}else
	output <- .C('bayesfit',as.integer(verbose),
                     as.integer(numInputs),as.integer(numCalibration),
                     as.integer(numPModel),as.integer(sizeData),
                     as.integer(sizeField),as.double(prob.prop),
                     as.double(object@mcmcMultmle),as.integer(object@method),
                     as.integer(n.iter),as.integer(nMH),
                     as.character(object@wd))
#cat('The results can be found on ',object@wd,'\n')
#	system(paste('ls ',object@wd,'*.out',sep=''))

	if (!is.null(object@calibrationnames) && length(object@calibrationnames)!=0){
		#info: Only in this case there are M-H steps
		cat("Acceptance rate:",scan(file=paste(object@wd,"rate.out",sep="/"), 
			quiet=T)[1],"\n")
	}
     	 	
    #//////
	post.thetaF<- read.table(file=paste(object@wd,"thetaF.out",sep="/"), header=F, 
                col.names=c("lambdaB","lambdaF"))
	
	if (!is.null(object@calibrationnames) && length(object@calibrationnames)!=0){
		post.calparams<- read.table(file=paste(object@wd,"ustar.out",sep="/"), 
			 header=F, col.names=object@calibrationnames)
		auxparams<- cbind(post.calparams,post.thetaF)
	}else auxparams<- post.thetaF
    #burnin and thining:
    auxparams<- auxparams[-(1:n.burnin),]
	auxparams<- auxparams[seq(from=1,to=dim(auxparams)[1], by=n.thin),]
	row.names(auxparams)<- 1:(dim(auxparams)[1])
	object@mcmcsample <- as.matrix(auxparams)
	
	unlink(paste0(object@wd,'/*'))

	return(object)
}

if(!isGeneric("bayesfit")) {
  setGeneric(name = "bayesfit",
             def = function(object, prior,mcmcMultmle=1,
		     prob.prop=0.5, method=2, n.iter, nMH=20, n.burnin=0,
         	 n.thin=1,verbose=FALSE,...) standardGeneric("bayesfit")
             )
}

setMethod("bayesfit", "SAVE", 
          function(object,prior, n.iter,...) {
			 #print(match.call(expand.dots=T))
			 #deprctcall <- .deprecate.parameters()
			 #print(as.call(deprctcall))
			 #object@bayesfitcall <- as.call(deprctcall)
			 bayesfit.SAVE(object=object, prior=prior,
			 n.iter = n.iter, mcmcMultmle=mcmcMultmle,
		     prob.prop=prob.prop, method=method, nMH=nMH, n.burnin=n.burnin,
         	 n.thin=n.thin,verbose=verbose)
          }
          )
          

