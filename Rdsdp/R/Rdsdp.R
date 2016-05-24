# Rdsdp.R
# Interface to the DSDP semidefinite programming library by Steve Benson, Yinyu Ye
# 
#
# Created: 08 December 2014
# Author: Zhisu Zhu


dsdp <- function(A, b, C, K, OPTIONS=NULL){
    data_filename = tempfile(pattern = "inputdata", tmpdir = tempdir(), fileext = ".dat-s")
    sedumi2sdpa(data_filename, A, b, C, K)

    options_filename=""
    if(!is.null(OPTIONS)){
      options_filename = tempfile(pattern = "inputoptions", tmpdir = tempdir(), fileext = ".options")
      write.options.file(options_filename, OPTIONS) 
    }

    # release memory
    remove(A,b,C)

    # solve the problem by dsdp5
   result = dsdp.readsdpa(data_filename,options_filename)

   if(file.exists(data_filename)) file.remove(data_filename)
   if(file.exists(options_filename)) file.remove(options_filename)

 
   STATS=as.list(result[[3]])
   # Mapping integer soltype back to DSDPSolutionType
   
   solutiontype=NULL
   if (STATS[[1]]==1) solutiontype="PDFeasible"
   else if (STATS[[1]]==3) solutiontype="Unbounded"
   else if (STATS[[1]]==4) solutiontype="Infeasible"
   STATS[[1]]=solutiontype
   
   # Assign the name to STATS output, should be consistent with Rreadsdpa.c 
   if(length(STATS)>3) names(STATS)=c("stype", "dobj","pobj","r","mu","pstep","dstep","pnorm")
   else names(STATS)=c("stype", "dobj","pobj")

   list(X=result[[1]],y=result[[2]],STATS=STATS)
   
}

dsdp.readsdpa <- function(sdpa_filename, options_filename="") {
  if(file.exists(sdpa_filename)) result <- .Call("dsdp", sdpa_filename, options_filename, PACKAGE="Rdsdp")
  else stop('File does not exist!')
  result
}



validate.data <- function(C,A,b,K) {

}


write.options.file <- function(options_filename, OPTIONS)
{
  write(paste(paste("-",names(OPTIONS),sep=""),OPTIONS), file=options_filename)
}

sedumi2sdpa <- function(filename,A,b,C,K)
{
  # Input A,b,C,K in sedumi format

  if(file.exists(filename)) file.remove(filename)
  
  # write mDim
  mDim=dim(A)[1]
  cat(dim(A)[1], file=filename, append=TRUE)
  cat("\n", file=filename, append=TRUE)
  
  # write nBlock
  nBlock = NULL
  blkSizes = NULL
  if(!is.null(K$l) && K$l != 0) 
  {
    nBlock = 1+length(K$s)
    blkSizes = c(-K$l,K$s)
  }else
  {
    nBlock = length(K$s)
    blkSizes = c(K$s)
  }
  cat( nBlock, file=filename, append=TRUE)
  cat("\n", file=filename, append=TRUE)
  
  # write blockStruct
  cat( blkSizes, file=filename, append=TRUE)
  cat("\n", file=filename, append=TRUE)
  
  # write b
  cat(as(b, "matrix"), file = filename, append=TRUE)
  cat("\n", file=filename, append=TRUE)
  
  # write C
  # block one for linear
  outputMat=NULL
  matnum=0;
  blocknum=0;
  curind = 1;
  currow=-C;
  if(!is.null(K$l) && K$l != 0){
    blocknum = blocknum+1;
    nzind = which(currow[c(1:K$l)]!=0, arr.ind=TRUE)
    if(length(nzind)>0){
      tmpoutput = cbind(matnum, blocknum, row=nzind, col=nzind, x=currow[nzind])
      outputMat=rbind(outputMat,tmpoutput)      
    }
    curind=curind+K$l;
  }
  for(conesize in K$s){
    blocknum = blocknum+1;
    blockmat = matrix(currow[curind:(curind+conesize*conesize-1)],conesize,conesize,byrow=TRUE)
    nzind = which(blockmat!=0, arr.ind=TRUE)
    if(dim(nzind)[1]>0){
      nzind = nzind[nzind[,"row"]<=nzind[,"col"],]
      if(is.null(dim(nzind)) ) nzind = t(as.matrix(nzind,byrow=FALSE))
      tmpoutput = cbind(matnum,blocknum,nzind,x=blockmat[nzind])
      outputMat=rbind(outputMat,tmpoutput)
    }
    # write(t(as.matrix(tmpoutput)), file = filename, append=TRUE)
    curind=curind+conesize*conesize;
  }
  
  # write A
  for (matnum in 1:mDim){
    blocknum=0;
    curind = 1;
    currow = A[matnum,]
    if(!is.null(K$l) && K$l != 0){
      blocknum = blocknum+1;
      nzind = which(currow[c(1:K$l)]!=0,arr.ind=TRUE)
      if(length(nzind)>0){
        tmpoutput = cbind(matnum, blocknum, row=nzind, col=nzind, x=currow[nzind])
        outputMat=rbind(outputMat,tmpoutput)
      }
      curind=curind+K$l;
    }
    for(conesize in K$s){
      blocknum = blocknum+1;
      blockmat = matrix(currow[curind:(curind+conesize*conesize-1)],conesize,conesize,byrow=TRUE)
      nzind = which(blockmat!=0, arr.ind=TRUE)
      if(dim(nzind)[1]>0){
        nzind = nzind[nzind[,"row"]<=nzind[,"col"],]
        if(is.null(dim(nzind)) ) nzind = t(as.matrix(nzind,byrow=FALSE))
        tmpoutput = cbind(matnum,blocknum,nzind,x=blockmat[nzind])
        outputMat=rbind(outputMat,tmpoutput)
        #write(t(as.matrix(tmpoutput)), file = filename, append=TRUE)       
      }
      curind=curind+conesize*conesize;
    }
  }
  # outputMat[,5]=-outputMat[,5]
  write.table(as.matrix(outputMat), file = filename, append=TRUE,row.names=FALSE,col.names=FALSE)
}
