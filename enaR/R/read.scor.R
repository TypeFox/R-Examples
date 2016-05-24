#' read.scor --- SCOR formatted file into R
#' in multiple formats
#' INPUT = file path
#' OUTPUT = network model in chosen format
#' S. Borrett and M. Lau | July 2011
#' ------------------------------------

read.scor <- function(file,from.file=TRUE,warn=FALSE){
  if (from.file){text <- readLines(file,warn=warn)}else{text <- file} # read in file
                                        #Partition the meta-data
  meta <- text[1]
                                        #Retrieve the number of node (n)
  n <- as.numeric(sub(' ','',substr(text[2],1,3)))
                                        #Determine the number of living nodes and create the living vector
  n.live <- as.numeric(sub(' ','',substr(text[2],4,6)))
  living <- c(rep(TRUE,n.live),rep(FALSE,(n-n.live)))
                                        #Retrieve vertex names
  vertex.names <- str_trim(text[3:(2+n)])
                                        # find negative ones (delimiters)
  br <- grep(pattern="( -1)", x=text)
  if (length(br) != 5){warning('Possible error in SCOR formatting')} # check expected number

  ## STORAGE
  B <- text[(3+n):(2+2*n)] # first cut at getting biomass data

  vertex.no <- as.numeric(sapply(B,function(x) (substr(x,1,3))))
  bm <- as.numeric(sapply(B,function(x) scifix(substr(x,5,nchar(x)))))
  storage <- data.frame("vertex"=vertex.no,"value"=bm) # final storage data

  ## INPUT

  if((br[2]-br[1])>1){
    inpt <- text[(br[1]+1):(br[2]-1)]
                                        #Condense into a function
    vertex.no <- as.numeric(sapply(inpt,function(x) substr(x,1,3)))
    z <- as.numeric(sapply(inpt,function(x) scifix(substr(x,5,nchar(x)))))
    inputs <- data.frame("vertex"=vertex.no,"value"=z) # final storage data

  } else {
    inputs <- NA
  }

  ## EXPORT

  if((br[3]-br[2])>1){
    export <- text[(br[2]+1):(br[3]-1)]
    vertex.no <- as.numeric(sapply(export,function(x) substr(x,1,3)))
    expt <- as.numeric(sapply(export,function(x) scifix(substr(x,5,nchar(x)))))
    exports <- data.frame("vertex"=vertex.no,"value"=expt) # final storage data

  } else {
    exports <- NA
  }


  ## RESPIRATION
  if((br[4]-br[3])>1){
    resp <- text[(br[3]+1):(br[4]-1)]
    vertex.no <- as.numeric(sapply(resp,function(x) substr(x,1,3)))
    resp <- as.numeric(sapply(resp,function(x) scifix(substr(x,5,nchar(x)))))
    respiration <- data.frame("vertex"=vertex.no,"value"=resp) # final storage data

  }else{
    respiration <- NA
  }


  ## FLOWS
                                        # assume there must be internal flows
  flows <- text[(br[4]+1):(br[5]-1)]
  strt <- as.numeric(sapply(flows,function(x) substr(x,1,3)))
  stp <- as.numeric(sapply(flows,function(x) substr(x,4,6)))
  value <- as.numeric(sapply(flows,function(x) scifix(substr(x,7,nchar(x)))))
  flows <- data.frame("tail"=strt,"head"=stp,"value"=value)
                                        #network data output type.
                                        # Note that the other output types depend on this sub-routine.

                                        #convert the flows to a matrix
  flow.mat <- array(0,dim=c(n,n))
  rownames(flow.mat) <- colnames(flow.mat) <- vertex.names
  for (i in seq(along=flows$tail)){
    flow.mat[flows$head[i],flows$tail[i]] <- flows$value[flows$tail==flows$tail[i]&flows$head==flows$head[i]]
  }
                                        #transpose flow matrix
  flow.mat <- t(flow.mat)
                                        #Vectorize the inputs
  input <- numeric(n)
  input[inputs$vertex] <- inputs$value
                                        #vectorize respiration and exports
  res <- numeric(n)
  exp <- numeric(n)

  if (any(is.na(exports)) == FALSE&any(is.na(respiration)) == FALSE){
    res[respiration$vertex] <- respiration$value
    exp[exports$vertex] <- exports$value
  }else if (any(is.na(exports))&any(is.na(respiration)) == FALSE){
    res[respiration$vertex] <- respiration$value
    exp <- rep(0,n)
  }else if (any(is.na(exports)) == FALSE&any(is.na(respiration))){
    exp[exports$vertex] <- exports$value
    res <- rep(0,n)
  }
  stor <- numeric(n)
  stor[storage$vertex] <- storage$value
                                        #produce an output vector given the values of export and respiration
  if (any(is.na(res)) == FALSE & any(is.na(exp))){
    output <- res
  }else if (any(is.na(res)) & any(is.na(exp)) == FALSE){
    output <- exp
  }else{output <- exp + res} #Outputs = respiration + exports for nea data type
                                        #introduce variables into the network format
  x <- pack(flow=flow.mat,input=input,export=exp,respiration=res,output=output,storage=stor,living=living)

  return(x)
}
