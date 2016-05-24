##########################################################################
## Predictreality Method
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 3, April 2013.
##
## Copyright (C) 2013-present Jesus Palomo, Gonzalo Garcia-Donato, 
##							  and Rui Paulo
##    
##########################################################################

predictreality.SAVE <- function(object, newdesign=NULL, n.burnin=0, n.thin=1, tol=1E-10,verbose=F){

	if ((dim(object@mcmcsample)[1]-n.burnin)<=2*n.thin){stop("Burnin and/or thinning too large for the size of your mcmc sample\n")}
	postsamples<- object@mcmcsample
	
	if (object@constant.controllables && (!is.null(newdesign))){
		stop("In the situation with constant controllable inputs, newdesign should be set to NULL.") 
	}

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
###############end of loading auxiliary funcs.	
	
#####
#Field data:
#Select those columns corresponding to the controllable inputs:
#  x<- object@df
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
#write to the files:
write.table(object@dm, file=paste(object@wd,"/model_inputs.dat",sep=""),
                    col.names=F, row.names=F)
write.table(object@ym, file=paste(object@wd,"/model_data.dat",sep=""),
              col.names=F, row.names=F)
#####

#####
#Mean responses:
  write.table(object@xm,
              file=paste(object@wd,"/predictionsII.design.M.old.matrix.dat",sep=""),
              col.names=F, row.names=F)
  write.table(object@xf,
              file=paste(object@wd,"/predictionsII.design.F.old.matrix.dat",sep=""), 
              col.names=F, row.names=F)
#####

# Design at which to predict
if(!object@constant.controllables){
  x.new <- as.data.frame(newdesign[,object@controllablenames])
  names(x.new)<- object@controllablenames
  write.table(x.new, file=paste(object@wd,"/inputs_real.dat",sep=""),
              col.names=F, row.names=F)
  #Mean responses:
  meannew <- model.matrix(object@meanformula, x.new)
  write.table(meannew,
              file=paste(object@wd,"/predictionsII.design.M.new.matrix.dat",sep=""),
              col.names=F, row.names=F)}
else{
  write.table(1, file=paste(object@wd,"/inputs_real.dat",sep=""),
              col.names=F, row.names=F)  
  write.table(1,file=paste(object@wd,"/predictionsII.design.M.new.matrix.dat",sep=""),
			              col.names=F, row.names=F)
	}  
#####

#####
#Write the MLE to the appropriate files
  write.table(object@mle$thetaM, file=paste(object@wd,"/thetaM_mle.dat", sep=""),
              col.names=F, row.names=F)
  write.table(object@mle$thetaL, file=paste(object@wd,"/thetaL_mle.dat", sep=""),
              col.names=F, row.names=F)	
  write.table(object@mle$thetaF, file=paste(object@wd,"/thetaF_mle.dat", sep=""),
              col.names=F, row.names=F)

#####
	  
#####
#Write the files related with the priors and calibration parameters
  n.iter <- dim(postsamples)[1]
  howmanycal <- dim(postsamples)[2]-2
  if (length(object@calibrationnames) != 0) {
	#####
	#Write the file for the prior
	write.table(object@prior, file=paste(object@wd,"/bounds.dat",sep=""),
			  col.names=F, row.names=F)
	  
	# calibration parameters which we are supposed to sample from posterior
	learn.names<- object@calibrationnames
  	learn <- rep(0,length(object@calibrationnames))
  	j <- 1
  	for(i in learn.names){
    	learn <- learn + (learn.names[j]== object@calibrationnames)
    	j <- j+1
  	}

  	write.table(learn, file=paste(object@wd,"/learn.dat",sep=""),
              col.names=F, row.names=F)
	write.table(postsamples[,1:howmanycal], file=paste(object@wd,"/ustar.dat",sep=""),
			  col.names=F, row.names=F)
  }
  
  write.table(postsamples[,-(1:howmanycal)],
              file=paste(object@wd,"/thetaF.dat",sep=""),
              col.names=F, row.names=F)

#####
#Values of the parameters:
	if (!object@constant.controllables){
    	numInputs<- length(c(object@controllablenames,object@calibrationnames))
		#dimension of field design unique inputs (NF in Rui''s notation)
		sizeField<- dim(object@xf)[1]
		sizeNewData <- dim(x.new)[1]
		sizeNewField<- dim(x.new)[1]		
	}
	else{
    	numInputs<- length(object@calibrationnames)
		#dimension of field design unique inputs (NF in Rui''s notation)
		sizeField<- 1
		sizeNewData <- 1
		sizeNewField<- 1				
	}

	#total number of inputs:
	#number of calibration inputs:
  	numCalibration<- length(object@calibrationnames)
	#dimension of the linear model for the mean of the GP prior
	#        (q in C notation)
	numPModel<- dim(object@xm)[2]
	#dimension of code design (Nold in C's notation)
	sizeData<- length(object@ym)
	#dimension of new code design (Nnew in C's notation)
	#tolerance for the pivoting algorithm
	tolerance <- tol
	workingPath<- object@wd

# load the C code
  if(!is.loaded("predict_reality")) {
    lib.file <- file.path(paste("predict_reality",
                                .Platform$dynlib.ext, sep=""))
    dyn.load(lib.file)
    cat(" -Loaded ", lib.file, "\n")
  }
  
#Call to the function:
  output <- .C('predict_reality',
			 as.integer(verbose),as.integer(numInputs),as.integer(numCalibration),
                 as.integer(numPModel),as.integer(sizeData),
                 as.integer(sizeField),as.integer(sizeNewData),as.integer(sizeNewField), as.double(tolerance),
                 as.integer(n.iter),
                 as.integer(n.burnin),as.integer(n.thin),
                 as.character(workingPath))
  #cat('The results can be found on ',workingPath,'\n')
  #system(paste('ls ',workingPath,'*.out',sep=''))
  
	#Results are being stored in an object called results
  results<- new("predictreality.SAVE")
  samples <- read.table(file=paste(object@wd,"real.dat",sep="/"), header=F)
  #first half of samples is prediction and second half is bias

  # To deprecate the unused parameters and to include in the 
  # call() all the default parameters not used in the call
  # to the function
  dprct <- .deprecate.parameters(call=sys.call(sys.parent(1)))
  #print(dprct)
  results@predictrealitycall<- as.call(dprct)
  results@modelpred<- as.data.frame(samples[,1:(dim(samples)[2]/2)])
  results@biaspred<- as.data.frame(samples[,((dim(samples)[2]/2)+1):dim(samples)[2]])

if (!object@constant.controllables){
  colnames(results@modelpred)<- rownames(newdesign)
  colnames(results@biaspred)<- rownames(newdesign)
  results@newdesign<- newdesign
}
else{
  colnames(results@modelpred)<- 1
  colnames(results@biaspred)<- 1
}

	unlink(paste0(object@wd,'/*'))
	return(results)
}

if(!isGeneric("predictreality")) {
  setGeneric(name = "predictreality",
             def = function(object, newdesign=NULL,
                           n.burnin=0, n.thin=1, tol=1E-10, verbose=FALSE, ...) 
			 standardGeneric("predictreality")
             )
}

setMethod("predictreality", "SAVE", 
          function(object, newdesign,...) {
            predictreality.SAVE(object=object, newdesign=newdesign ,n.burnin=n.burnin, n.thin=n.thin, 
								tol=tol,verbose=verbose)
          }
          )

setMethod("show","predictreality.SAVE",
        function(object){
		smx <- summary.predictreality.SAVE(object)
		show.summary.predictreality.SAVE (smx) }
)


summary.predictreality.SAVE<- function(object){
	result<- new("summary.predictreality.SAVE")
	result@call<- object@predictrealitycall
	mysummary<- function(x){
		round(c(mean(x),sd(x),median(x),quantile(x,probs=c(.025,.975))),3)
	}

	real<- object@modelpred+object@biaspred

	result@biascorr<- matrix(0, nrow=dim(real)[2], ncol=5)
	result@biascorr<- t(apply(real, MARGIN=2, FUN=mysummary))
	dims <- dim(real)

	#Symmetric bounds are used:
	tmp <- real - matrix(result@biascorr[,1],ncol=dims[2],nrow=dims[1],byrow=T)
	tmp <- apply(tmp,2,abs)
	tau.real <- apply(tmp,2,quantile,0.95)
	result@biascorr[,4]<- round(result@biascorr[,1]-tau.real,3)
	result@biascorr[,5]<- round(result@biascorr[,1]+tau.real,3)
		
	rownames(result@biascorr)<- colnames(object@biaspred)
	colnames(result@biascorr)<- c("Mean", "SD", "Median", "2.5%", "97.5%")

	result@biaspred<- matrix(0, nrow=dim(real)[2], ncol=5)
	result@biaspred<- t(apply(object@biaspred, MARGIN=2, FUN=mysummary))
	rownames(result@biaspred)<- colnames(object@biaspred)
	colnames(result@biaspred)<- c("Mean", "SD", "Median", "2.5%", "97.5%")
	return(result)
}

setMethod("summary","predictreality.SAVE",
	function(object){ summary.predictreality.SAVE (object) }
)


show.summary.predictreality.SAVE<- function(object){
	cat("\n")
	cat("---------------\n")
	cat("Call to predictreality:\n")
	print(object@call)
	cat("---------------\n")
	cat("Summary of results:\n")
	cat("---Bias corrected prediction-----|---------Bias function-----------\n")
	print(cbind(object@biascorr, object@biaspred))
}

setMethod("show","summary.predictreality.SAVE",
        function(object){show.summary.predictreality.SAVE (object) }
)

plot.predictreality.SAVE<- function(x, option="biascorr", ...){
    par(mfrow=c(1,1))
	if (option=="biascorr"){
		summaries<- summary.predictreality.SAVE(x)@biascorr
		plot(x=1:dim(summaries)[1], y=as.vector(summaries[,"Mean"]), xlab="Input points", 
				ylim=c(min(summaries[,"2.5%"])*.95,max(summaries[,"97.5%"])*1.05),ylab="Bias corrected prediction", type="n")
		points(x=1:dim(summaries)[1], y=as.vector(summaries[,"Mean"]), cex=0.8)
		for (i in 1:(dim(summaries)[1])){
			lines(x=c(i,i),y=c(summaries[i,"2.5%"],summaries[i,"97.5%"]))
		}
	}else
	 if (option=="biasfun"){
		summaries<- summary.predictreality.SAVE(x)@biaspred
		plot(x=1:dim(summaries)[1], y=as.vector(summaries[,"Mean"]), xlab="Input points", 
		ylim=c(min(summaries[,"2.5%"])*.95,max(summaries[,"97.5%"])*1.05),ylab="Bias function", type="n")
		points(x=1:dim(summaries)[1], y=as.vector(summaries[,"Mean"]), cex=0.8)
		for (i in 1:(dim(summaries)[1])){
			lines(x=c(i,i),y=c(summaries[i,"2.5%"],summaries[i,"97.5%"]))
		}
	 }
	 else{stop("Invalid option\n")}
}

setMethod("plot",signature(x="predictreality.SAVE", y="missing"), 
	function(x, option="biascorr", ...) {
	plot.predictreality.SAVE(x = x, option = option, ...)
}
)