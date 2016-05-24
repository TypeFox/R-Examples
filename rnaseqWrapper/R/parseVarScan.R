## Parse VarScan ouputs

parseVarScan <- function(
  file, ## a tab separated filename with the varscan output or a data.frame/matrix holding the information
  sampleNames = NULL, ## Character vector with the sample names
  ignoreIndels = TRUE ## Logical, should indels be omitted from the final output
  ){
  
  ## Read in varscan:
  if(class(file) == "matrix" | class(file) == "data.frame"){
    fullVarScan <- file
  } else if (class(file) == "character"){
    if( length(file) ==1){
      fullVarScan <-  read.table(file,header=TRUE,stringsAsFactors=FALSE,sep="\t")      
    } else{
      stop("length(file) > 1")
    }
  } else{
    stop("file must be a data.frame, matrix, or filename to read in")
  }

  
  ## Remove indels
  if(ignoreIndels){
    matchesIndel <- grep("[+-]",fullVarScan$Var)
    if(length(matchesIndel > 0)){
      limVarScan <- fullVarScan[-matchesIndel,]
    } else {
      limVarScan <- fullVarScan 
    }
  } else{
    limVarScan <- fullVarScan
  }
 
  
  
  ###########################################
  ###### Break up the pool information ######
  ###########################################
  allInfo <- as.list(data.frame(do.call(rbind,strsplit(as.character(limVarScan$Cons.Cov.Reads1.Reads2.Freq.P.value),":")),stringsAsFactors=FALSE))
  names(allInfo) <- unlist(strsplit("Cons.Cov.Reads1.Reads2.Freq.Pvalue",".",fixed=TRUE))
  allInfo$Freq <- gsub("%","",allInfo$Freq,fixed=TRUE)
  
  suppressWarnings(tempNumerics <- lapply(allInfo,as.numeric) )
  propNa <- lapply(tempNumerics,function(x) {sum(is.na(x))/length(x)})
  
  for(k in 1:length(allInfo)){
    if(propNa[k] < 0.9){
      allInfo[[k]] <- tempNumerics[[k]]
    }
  }
  
  combinedAllInfo <- as.data.frame(allInfo)
  
  
  ###########################################
  ##### Break up the sample information #####
  ###########################################
  sampInfo <- as.list(data.frame(do.call(rbind,strsplit(as.character(limVarScan$Cons.Cov.Reads1.Reads2.Freq.P.value.1)," ")),stringsAsFactors=FALSE))
  
  if(is.null(sampleNames)){
    sampleNames <- c(letters,LETTERS,paste(rep(letters,each=26),LETTERS,sep=""))[1:length(sampInfo)]
  }
  names(sampInfo) <- sampleNames
  
  
  brokenSampInfo <- lapply(
    sampInfo,function(y) {
      broken <- as.list(data.frame(do.call(rbind,strsplit(y,":")),stringsAsFactors=FALSE))
      names(broken) <- unlist(strsplit("Cons.Cov.Reads1.Reads2.Freq.Pvalue",".",fixed=TRUE))
      broken$Freq <- gsub("%","",broken$Freq,fixed=TRUE)
      
      suppressWarnings(tempNumerics <- lapply(broken,as.numeric))
      propNa <- lapply(tempNumerics,function(x) {sum(is.na(x))/length(x)})
      
      for(k in 1:length(broken)){
        if(propNa[k] < 0.9){
          broken[[k]] <- tempNumerics[[k]]
        }
      }
      
      #     hist(broken$Freq)
      out <- as.data.frame(broken)
      return(out)
      
    })
  
  combinedSampInfo <- do.call(cbind,brokenSampInfo)
  #
  
  
  ###########################################
  ###### Combine the broken up methods ######
  ###########################################
  
  parsedVarScan <- cbind(limVarScan[,-c(5,6,11)],combinedAllInfo,combinedSampInfo)
  
  return(parsedVarScan)
  
}