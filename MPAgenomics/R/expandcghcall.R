#
# This is the ExpandCGHcall function from CGHcall package (GPL License).
# 
# Expands result from \link{CGHcall} function to CGHcall object.
#
# @title Expands result fron CGHcall to CGHcall object.
# 
# @param listcall  List object; output of function \link{CGHcall}.
# @param inputSegmented a list.See Details section in \link{CGHcall} function.
# @param digits Number of decimal digits to be saved in the resulting call object. Allows for saving storage space.
# @param divide Number of batches to divide the work load in. Larger values saves memory, but requires more computing time.
# @param memeff  When set to TRUE, memory efficient mode is used: results are written in batches to multiple external files. If FALSE, one output object is provided.
# @param fileoutpre	Only relevant when memeff=TRUE. Define prefix for output file names.
# @param CellularityCorrectSeg	If TRUE, corrects segmented and normalized values for cellularity as well.
# @param verbose If TRUE, print some details.
#
# @return A list containing the same elemnts as inputSegmented and
# \describe{
#   \item{calls}{A matrix, of the same size as inputSegmented$copynumber matrix, containing the label of each point.
#   -2=double loss, -1=loss, 0=normal, 1=gain, 2=amplification.}
#   \item{probdloss}{(if you have ran \link{CGHcall} with nclass=5) A matrix of the same size as inputSegmented$copynumber matrix. It contains the probability for each segmented copynumber to be a double loss.}
#   \item{probloss}{A matrix of the same size as inputSegmented$copynumber matrix. It contains the probability for each segmented copynumber to be a loss.}
#   \item{probdnorm}{A matrix of the same size as inputSegmented$copynumber matrix. It contains the probability for each segmented copynumber to be normal.}
#   \item{probdgain}{A matrix of the same size as inputSegmented$copynumber matrix. It contains the probability for each segmented copynumber to be a gain.}
#   \item{probdamp}{(if you have ran \link{CGHcall} with nclass=4 or 5) A matrix of the same size as inputSegmented$copynumber matrix. It contains the probability for each segmented copynumber to be an amplification.}
# }
# 
#  
# @details  It allows more memory efficient handling of large data objects. 
# If R crashes because of memory problem, we advise to set memeff = TRUE and increase the value of divide. 
# When multiple files are output (in case of memeff=TRUE) the function combine may be used to combine CGHcall objects.
#
#
# @author Sjoerd Vosse & Mark van de Wiel
# @references Mark A. van de Wiel, Kyung In Kim, Sjoerd J. Vosse, Wessel N. van Wieringen, Saskia M. Wilting and Bauke Ylstra. CGHcall: calling aberrations for array CGH tumor profiles. Bioinformatics, 23, 892-894.
# 
# @export
ExpandCGHcall <- function(listcall,inputSegmented, digits=3,divide=4, memeff = FALSE, fileoutpre="Callobj_",CellularityCorrectSeg=TRUE,verbose=TRUE)
{

  timeStarted <- proc.time()
  posteriorfin2 <- listcall[[1]];nclone<-listcall[[2]];nctot <- listcall[[3]];nclass <- listcall[[4]];regionsprof<-listcall[[5]]
  cellularity <- listcall[[7]]
  
  adjustForCellularity <- function(matrix, cellularity,pmode) {
    if(verbose)
    {
      if(pmode=="seg") cat("Adjusting segmented data for cellularity ... \n") else cat("Adjusting normalized data for cellularity ... \n")
    }
        result  <- c();
        adjustCellularity <- function(value, cellularity) {
            corrected   <- (2^value / cellularity - (1 - cellularity) / cellularity)
            if (corrected < 2^(-5)) {
                corrected <- 2^value;
            }
            new.value   <- log2(corrected)
            return(new.value)
        }
        for (i in 1:ncol(matrix)) {
            if(verbose)
              cat("Cellularity sample", i, ": ", cellularity[i], "\n");
            if (cellularity[i] < 1) {
                new.column  <- sapply(matrix[,i], adjustCellularity, cellularity[i]);
                result      <- cbind(result, new.column);
            } else {
                result      <- cbind(result, matrix[,i]);
            }
        }
        colnames(result) <- colnames(matrix)
        return(result);
    }  
  
  #digits=3;divide=4; memeff = FALSE; fileoutpre="Callobj_"
  if(CellularityCorrectSeg) {
    inputSegmented$copynumber <- adjustForCellularity(inputSegmented$copynumber,cellularity,pmode="seg")
    inputSegmented$segmented <- adjustForCellularity(inputSegmented$segmented,cellularity,pmode="norm")
    }
  inputSegmented$copynumber<-round(inputSegmented$copynumber,digits)
  inputSegmented$segmented<-round(inputSegmented$segmented,digits)
  if (divide > nctot) divide <- nctot
  nperturn <- floor(nctot/divide)
  
  for(part in 1:divide)
  {
    #print(part)
    if (part < divide) whprof <- ((part-1)*nperturn+1):(part*nperturn) else whprof <- ((part-1)*nperturn+1):nctot
    #IS <- inputSegmented[,whprof]
    IS=list(copynumber=inputSegmented$copynumber[,whprof],segmented=inputSegmented$segmented[,whprof],
            startPos=inputSegmented$startPos,chromosome=inputSegmented$chromosome,featureNames=inputSegmented$featureNames,
            sampleNames=inputSegmented$sampleNames[whprof])
    nc <- length(whprof)

    dataprob<-array(0,c(nclone,nclass*nc))
    gc()#print(gc())
    for (k in 1:nc) 
    {      
        post        <- (posteriorfin2[posteriorfin2[,1]==whprof[k],,drop=FALSE])[,-1,drop=FALSE] #add ,drop=FALSE
        regionsk    <- (regionsprof[regionsprof[,1]==whprof[k],,drop=FALSE])[,-1,drop=FALSE] #add ,drop=FALSE
        nregk       <- nrow(post)
        probs       <- c()
        for (i in (1:nregk)) 
        {
            regl    <- regionsk[i,2]-regionsk[i,1]+1
            togeth  <- post[i,(1:nclass)]
            probs   <- c(probs,rep(togeth,regl))
        }
        probs <- round(probs,digits)
        allprobs    <- matrix(probs, ncol=nclass, byrow=TRUE)
      #  datk        <- datareg[[1]][,k]
        dataprob[,(nclass*(k-1)+1):(nclass*k)] <- allprobs
        rm(probs,post,regionsk,togeth);
        gc()#print(gc())
    }
    gc()#print(gc())
   
    #ncolscl     <- ncol(normalizedData) #25/11/09 moved upward for efficiency reasons
    classify.res        <- array(NA,c(nclone,nc))
    for (i in 1:nc) 
    {
        inc <- 0
        if (nclass==5) 
        {
            prob.dl.ind    <- dataprob[,(nclass*(i-1)+1)]
            ticksdl        <- which(prob.dl.ind >= 0.5)
            ltdl              <- length(ticksdl)
            inc <- 1
        }    
   #     genomdat        <- dataprob[,(ncpp*(i-1)+1)] #removed 16/7/10
        prob.loss.ind   <- dataprob[,(nclass*(i-1)+1+inc)]
        prob.none.ind   <- dataprob[,(nclass*(i-1)+2+inc)]
        prob.gain.ind   <- dataprob[,(nclass*(i-1)+3+inc)]
        if (nclass>=4) 
        {
            prob.amp.ind    <- dataprob[,(nclass*(i-1)+4+inc)]
            ticksamp        <- which(prob.amp.ind > 0.5)
            lt              <- length(ticksamp)
        }        
        
        if (nclass==3) {
            classify.res[(1:nclone),i] <- 
            (2*(as.numeric(prob.gain.ind>prob.loss.ind))-1)*(as.numeric(apply(cbind(prob.loss.ind,prob.gain.ind),1,max)>prob.none.ind))
        }
        
        if (nclass==4) {
            prob.gain.amp <- prob.gain.ind + prob.amp.ind
            classify.res[(1:nclone),i] <- (as.numeric(prob.amp.ind>0.5)+1)*(2*(as.numeric(prob.gain.amp>prob.loss.ind))-1)*(as.numeric(apply(cbind(prob.loss.ind,prob.gain.amp),1,max)>prob.none.ind))
        }    
        if (nclass==5) {
            prob.gain.amp <- prob.gain.ind + prob.amp.ind
            prob.l1.l2 <- prob.dl.ind + prob.loss.ind
            classify.res[(1:nclone),i] <- 
   (as.numeric(prob.dl.ind>0.5)+1)*(as.numeric(prob.amp.ind>0.5)+1)*(2*(as.numeric(prob.gain.amp>prob.l1.l2))-1)*(as.numeric(apply(cbind(prob.l1.l2,prob.gain.amp),1,max)>prob.none.ind))
        }            
    }
    gc()#print(gc())
    ncolprob <- ncol(dataprob)
    neworder<-as.vector(sapply(1:nclass,function(x)seq(x,ncolprob,by=nclass)))
    dataprob   <- dataprob[,neworder]
    calls       <- .assignNames(classify.res, IS)
    rm(classify.res);
    gc()#print(gc());
    if(nclass==5){
    probdloss <- .assignNames(dataprob[,1:nc,drop=FALSE], IS)
    dataprob <- dataprob[,-(1:nc),drop=FALSE]
    gc()#print(gc())
    }
    
    probloss <- .assignNames(dataprob[,1:nc,drop=FALSE], IS)
    dataprob <- dataprob[,-(1:nc),drop=FALSE]
    gc()#print(gc())
    probnorm <- .assignNames(dataprob[,1:nc,drop=FALSE], IS)
    dataprob <- dataprob[,-(1:nc),drop=FALSE]
    gc()#print(gc())
    probgain <- .assignNames(dataprob[,1:nc,drop=FALSE], IS)
    dataprob <- dataprob[,-(1:nc),drop=FALSE]
    gc()#print(gc())
    if(nclass>=4) probamp <- .assignNames(dataprob[,1:nc,drop=FALSE], IS)
    rm(dataprob);
    gc()#print(gc())
    if (nclass == 5) {assayData <-list(calls=calls,probdloss=probdloss, probloss=probloss,probnorm=probnorm,probgain=probgain,probamp=probamp)} 
    if (nclass == 4) {assayData <-list(calls=calls,probloss=probloss,probnorm=probnorm,probgain=probgain,probamp=probamp)} 
    if (nclass == 3) {assayData <-list(calls=calls,probloss=probloss,probnorm=probnorm,probgain=probgain) }
    rm(probloss,probnorm,probgain); if(nclass>=4) rm(probamp);if(nclass==5) rm(probdloss)
    gc()#print(gc())

    #result0  <- CGHcall:::.callFromSeg(IS, assayData)  
    result0=c(IS,assayData)
    data        <- data.frame(Cellularity=cellularity[whprof])
    dimLabels   <- c("sampleNames", "sampleInfo")
    metadata    <- data.frame(labelDescription=c("Proportion of tumor cells"), row.names=c("Cellularity"))
    #pd <- new("AnnotatedDataFrame", data=data, dimLabels=dimLabels, varMetadata=metadata)   
    #sampleNames(pd) <- colnames(copynumber(IS))
    #phenoData(result0) <- pd 
   
    
    if(!memeff)
    {  
        if (part==1) 
          {result <- result0} 
        else 
        {
          #result <- combine(result,result0)  #merge 2 objet
          for(i in 1:length(result))
          {
            result[[i]]=cbind(result[[i]],result0[[i]])
          }            	
        }
    } 
    else 
    {
        fileout <- paste(fileoutpre,"profiles",whprof[1],"to",whprof[length(whprof)],".Rdata",sep="")  
        save(result0,file=fileout)
    }
    rm(assayData,result0);
    gc()#print(gc())
  }
  
  if(verbose)
    cat("FINISHED!\n")
  timeFinished <- round((proc.time() - timeStarted)[1] / 60)
  if(verbose)
    cat("Total time:", timeFinished, "minutes\n")
  if (memeff) print("Results printed to separate files") else return(result)
}
