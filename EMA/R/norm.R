#####
##
## Normalisation functions
##
#####

expFilter <- function(data, threshold=3.5, p=0.01, graph=TRUE){
  ##vecMax <- apply(data, 1, max)
 
  if (graph){
    par(font.lab=2) 
    hist(data, breaks=200, main="Data distribution", xlab="Expression Level")
    abline(v=threshold, col=2)
    mtext(threshold, at=threshold, side=1, las=1, col=2, cex=0.6)

  }
  vecL <- apply(data,1, function(x){length(which(x>threshold))})
  cat("Keep probes with at least",ceiling(ncol(data)*p),"sample(s) with an expression level higher than",threshold,"\n")
  data <- data[which(vecL>=ceiling(ncol(data)*p)),]
  ##data <- data[which(vecMax>threshold),]
}

normAffy <- function(filenames, celfile.path, method=c("GCRMA","RMA","MAS5"), cdfname = NULL, rmaffx=TRUE, fast=TRUE){
    ## check if arguments are ok
    if( missing(filenames) & missing(celfile.path) )
        stop("** FAILURE : 'filenames' or 'celfile.path' arguments are empty")
    if(class(celfile.path)!="character")
        stop("** FAILURE : 'celfile.path' argument is not a character class")
   
    celfile.path <- file.path(celfile.path)
    if (missing(filenames)){
        filenames <- list.celfiles(path=celfile.path, full.names=FALSE)
    }
    else{
        filenames<-as.vector(filenames)
    }
    print (filenames)
    method <- match.arg(method)
    
    if(method == "RMA"){
        print("--- RMA normalization ---")
        normData <- justRMA(filenames=filenames, celfile.path=celfile.path, cdfname=cdfname)
        normData <- exprs(normData)
    }
    if (method == "MAS5"){
        rawData <- ReadAffy(filenames=filenames, celfile.path=celfile.path, cdfname=cdfname)
        print("--- MAS5 normalization ---")
        normData <- mas5(rawData)
        normData <- log2(exprs(normData))
    }
    if (method == "GCRMA"){
        print("--- justGCRMA normalization ---")
        normData <- justGCRMA(filenames=filenames, celfile.path=celfile.path, type="affinities", cdfname=cdfname, fast=fast)
        normData <- exprs(normData)
    }
    
    if(rmaffx){
        normData=normData[-grep("AFFX", rownames(normData)),]
    }
    
    return(normData)
}

