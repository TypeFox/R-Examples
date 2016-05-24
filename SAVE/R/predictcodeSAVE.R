##########################################################################
## Predictcode Method
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 3, April 2013.
##
## Copyright (C) 2013-present Jesus Palomo, Gonzalo Garcia-Donato, 
##							  and Rui Paulo
##    
##########################################################################

predictcode.SAVE <- function(object, newdesign, n.iter, sampledraws, tol,verbose=FALSE){
#####
#Model data:
#write to the files:
  write.table(object@ym,
              file=paste(object@wd,"/model_data.dat",sep=""),
              col.names=F, row.names=F)
  write.table(object@dm, file=paste(object@wd,"/model_inputs.dat",sep=""),
              col.names=F, row.names=F)

#####
#New design
  if (!object@constant.controllables){ #Case 1
	x.new <- as.data.frame(newdesign[,c(object@controllablenames,object@calibrationnames)])
  	names(x.new) <- c(object@controllablenames,object@calibrationnames)
  } else{ #Case 2: Constant controllable inputs
	x.new <- as.data.frame(newdesign[,object@calibrationnames])
    names(x.new) <- object@calibrationnames
  }
#write to the files:
  write.table(x.new, file=paste(object@wd,"/inputs_pure.dat",sep=""),
              col.names=F, row.names=F)
#####

#####
#Mean response for both new and old designs:
  write.table(model.matrix(object@meanformula, as.data.frame(object@dm)),
              file=paste(object@wd,"/predictionsI.design.old.matrix.dat",sep=""),
              col.names=F, row.names=F)
	  
  write.table(model.matrix(object@meanformula,x.new),
              file=paste(object@wd,"/predictionsI.design.new.matrix.dat",sep=""),
              col.names=F, row.names=F)
#####

#####
#Write the MLE to the appropriate files
  write.table(object@mle$thetaM, file=paste(object@wd,"/thetaM_mle.dat", sep=""),
              col.names=F, row.names=F)
  write.table(object@mle$thetaL, file=paste(object@wd,"/thetaL_mle.dat", sep=""),
              col.names=F, row.names=F)	
	  
#####
#Values of the parameters:
	#total number of inputs:
    if (!object@constant.controllables){
        numInputs<- length(c(object@controllablenames,object@calibrationnames))
    }
    else{
        numInputs<- length(object@calibrationnames)
    }
	#dimension of the linear model for the mean of the GP prior
	#        (q in C notation)
	numPModel<- dim(object@xm)[2]
	#dimension of code design (Nold in C's notation)
	sizeData<- length(object@ym)
	#dimension of new code design (Nnew in C's notation)
	sizeNewData <- dim(x.new)[1]
	#length of simulation
	if (!sampledraws){n.iter<- 1} #never used, only mean is returned
	sizeSimulation <- n.iter
	#tolerance for the pivoting algorithm
	tolerance <- tol
	# Parameter to indicate whether stage I parameters were estimated with  mcmc or 	with mle  
	ismleormcmc <- 1
	workingPath<- object@wd

	# load the C code
	if(!is.loaded('predict_code')) {
		lib.file <- file.path(paste("predict_code",
                                            .Platform$dynlib.ext, sep=""))
		dyn.load(lib.file)
		cat(" -Loaded ", lib.file, "\n")
	}
	#Call to the function:
	output <- .C('predict_code',as.integer(verbose),
			as.integer(numInputs),
                   as.integer(numPModel),as.integer(sizeData),
                     as.integer(sizeNewData),
			as.double(tolerance),as.integer(sizeSimulation),
			as.integer(ismleormcmc),as.integer(sampledraws),
                     as.character(workingPath))
	#cat('The results can be found on ',workingPath,'\n')
	#system(paste('ls ',workingPath,'*.dat',sep=''))
	#Results are being stored in an object called results
	#cat("Creating the object results\n")
    results<- new("predictcode.SAVE")
	#cat("Created the object results\n")

	# To deprecate the unused parameters and to include in the 
	# call() all the default parameters not used in the call
	# to the function
	dprct <- .deprecate.parameters(call=sys.call(sys.parent(2)))
	results@predictcodecall<- as.call(dprct)
    results@mle<- object@mle
	#print(sampledraws)
	if (sampledraws) {
		#cat('About to read path_pure.dat\n')
		tryCatch({
			sampl<-read.table(file=paste(object@wd,"path_pure.dat",sep=""),header=F)
					 },
			error= function(ex) {
			stop ("An error occurred when sampled in predictcode. Nothing to be read. Pure samples")
			}
		)
	}else { # If not samples are requested we set the sample mean
		#cat("Setting the mean as samples\n")
		tryCatch({
			sampl <- t(scan(file=paste(object@wd,"mean_vector.dat",sep=""),quiet=T))					 },
			error= function(ex) {
			stop ("An error occurred when sampled in predictcode. Nothing to be read. Mean vector")
			}
		)
	}
	
	names(sampl) <- rownames(newdesign) #gon:20-2-13
	#cat("is samples a data.frame ") 
	#print(is.data.frame(sampl))
	results@samples <- as.data.frame(sampl)
    results@newdesign <- newdesign
	results@modelmean<- as.vector(scan(file=paste(object@wd,"mean_vector.dat",sep=""),quiet=T))
	results@covmat<- as.matrix(read.table(file=paste(object@wd,"cov_mat.dat",sep="")),header=F)

	unlink(paste0(object@wd,'/*'))

	return(results)
}

if(!isGeneric("predictcode")) {
  setGeneric(name = "predictcode",
             def = function(object,newdesign,n.iter=1000,sampledraws=T, tol=1E-10,...) standardGeneric("predictcode")
             )
}

setMethod("predictcode", "SAVE",
	#signature(object="SAVE",newdesign="data.frame",n.iter="missing",sampledraws="missing"),
          definition=function(object, newdesign,n.iter=1000,sampledraws=T,tol=1E-10,verbose=FALSE) {
            predictcode.SAVE(object = object, newdesign=newdesign, n.iter=n.iter,
				sampledraws=sampledraws,tol=tol,verbose=verbose)
          }
          )

setMethod("show","predictcode.SAVE",
        function(object){
		smx <- summary.predictcode.SAVE(object)
		show.summary.predictcode.SAVE (smx) }
)


summary.predictcode.SAVE<- function(object){
	result<- new ("summary.predictcode.SAVE")
	result@call<- object@predictcodecall
	# NULL assignment is not permuted ---- error ----------
	#result@modelpred<- NULL
	#we work on the theoretical results
		result@summariesmodelpred<- matrix(0, nrow=length(object@modelmean), ncol=4)
		theomean<- object@modelmean
		theovar<- diag(object@covmat)
		
		result@summariesmodelpred[,1]<- theomean
		result@summariesmodelpred[,2]<- sqrt(theovar)
		result@summariesmodelpred[,3]<- theomean-1.96*sqrt(theovar)
		result@summariesmodelpred[,4]<- theomean+1.96*sqrt(theovar)
		result@summariesmodelpred<- round(result@summariesmodelpred,3)
		rownames(result@summariesmodelpred)<- colnames(object@samples)
		colnames(result@summariesmodelpred)<- c("Mean", "SD", "2.5%", "97.5%")
	
	return(result)
}

setMethod("summary","predictcode.SAVE",
	function(object){ summary.predictcode.SAVE (object) }
)


show.summary.predictcode.SAVE<- function(object){
	cat("\n")
	cat("---------------\n")
	cat("Call to predictcode:\n")
	print(object@call)
	cat("---------------\n")
	cat("Summary of results:\n")
	if (!is.null(object@summariesmodelpred) || length(object@summariesmodelpred!=0)){
		cat("---Model prediction-----\n")
		print(object@summariesmodelpred)
	}
}

setMethod("show","summary.predictcode.SAVE",
        function(object){show.summary.predictcode.SAVE (object) }
)


plot.predictcode.SAVE<- function(x, ...){
	summaries<- summary.predictcode.SAVE(x)@summariesmodelpred
	plot(x=1:dim(summaries)[1], y=as.vector(summaries[,"Mean"]), xlab="Input points", 
		ylim=c(min(summaries[,"2.5%"])*.95,max(summaries[,"97.5%"])*1.05),ylab="Emulation", type="n")
	points(x=1:dim(summaries)[1], y=as.vector(summaries[,"Mean"]), cex=0.8)
	for (i in 1:(dim(summaries)[1])){
		lines(x=c(i,i),y=c(summaries[i,"2.5%"],summaries[i,"97.5%"]))
		}
}

setMethod("plot",signature(x="predictcode.SAVE",y="missing"), 
	function(x, ...) {
	plot.predictcode.SAVE(x = x, ...)
}
)

#show.predictcode.SAVE<- function(x, ...){
#	show.summary.predictcode.SAVE(x)
#}
