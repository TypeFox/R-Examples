options("warn"=-1)

msg <- file("message.log", open="wt")
sink(msg,type="message",append=TRUE)
#sink(msg,type="output")

library(Biobase)
library(methods)

cmd_args <- commandArgs(trailingOnly=TRUE)

makeEset <-function(cmd_args){
  set <- unlist(strsplit(cmd_args[1],"\\."))[2]

  index <- grep("ID_REF",cmd_args)
  probes <- scan(cmd_args[index],sep="\t",what="complex",skip=3)
  #print(cmd_args)
  cmd_args <- cmd_args[-index]
	
  len <- length(cmd_args)
  ncols <- as.numeric(cmd_args[len-1])
  #print(ncols)
  flag <- as.numeric(cmd_args[len])
	
  nrows <- unlist(strsplit(cmd_args[len-2],"-"))[1]
  colnames <- scan(cmd_args[1],what="complex",nlines=1,skip=1,sep="\t")

  eset <- new("ExpressionSet")
  for(i in 1:(length(cmd_args)-3)){
		
    name <- unlist(strsplit(cmd_args[i],"\\."))[1]
    data <- scan(cmd_args[i],skip=3,sep="\t",what="numeric")
    data <- matrix(data,byrow=TRUE,ncol=ncols)
	
    ## if the values are characters such as in ABS_CALL all values will be NA
    data1 <- apply(data,2,as.numeric)
    r <- dim(data)[1]
    c <- dim(data)[2]
	
    ## if all values are not NA then the data will be numeric else it will be character data
    if(length(which(is.na(data1)))!=r*c){
      data <- data1
    }
	
    if(nrow(data)==nrows)
      {
        colnames(data) <- colnames	
        rownames(data) <- probes	
        
        if(name=="VALUE"){
          assayDataElement(eset,"exprs") <- data
        }else{
          assayDataElement(eset,name) <- data
        }
      }
  }
  
  cat(paste("Creating feature data",date()))
  if(!is.na(match("gplnew.annotation",list.files()))){
    gpl <- scan("gplnew.annotation",what="complex",sep="\t",quote="") ### Annotation provided by submitter
    ncolumns <- scan("gplnew.annotation",what="complex",sep="\t",quote="",nlines=1)
    gpl <- matrix(gpl,byrow=TRUE,ncol=length(ncolumns))
  }else{
    gpl <- scan("gplnew",what="complex",sep="\t",quote="") ### Annotation provided by NCBI GPLxxx.annot
    gpl <- matrix(gpl,byrow=TRUE,ncol=10)
  }
  colnames(gpl) <- gpl[1,]
  gpl <- gpl[-1,]
	
  probes <- data.frame(probes)
  colnames(probes) <- "ID"
  rownames(probes) <- probes[,1]
  i <- match(probes[,1],gpl[,1])
  gplannotation <- cbind(probes,gpl[i,-1])
  ## some GPLs have duplicate column names
  ## see http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL8227 for an example
  colnames(gplannotation) <- make.unique(colnames(gplannotation))
  fData(eset) <- gplannotation
	
  if(flag){ 
    ## if number of sets is more than 1
    compendiumESet <- eset
    for(j in (flag-1):1){
      load(paste("esetset",j,".RData",sep=""))
      compendiumESet <- c(eset,compendiumESet)
    }
    save(compendiumESet,file="eset.RData")
  }else{
    compendiumESet=eset
    save(eset,file=paste("eset",set,".RData",sep=""))
    save(compendiumESet,file="eset.RData")
  }	
}

possibleError <- tryCatch(
  makeEset(cmd_args),
  error = function(e) e
)
print(possibleError)
