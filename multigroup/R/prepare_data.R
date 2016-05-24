#' @title Preparing data
#' 
#' @description 
#' To preparing data for doing the analysis
#' 
#' @param The input parameters of the main function
#' @return Prepared data. Variables are centered within each group and block. Other preprocessing are voluntary
#' @export
#' @keywords internal
prepare_data= function(Data, Group, nBlock,   Block.name=NULL,
                       ScaleGroup=FALSE, ScaleDataA=FALSE, ScaleDataB=FALSE, norm=FALSE){
  
  rownames(Data) = Group                 #----rownames of data=groups
  M = length(levels(Group))              #----number of groups: M
  n = as.vector(table(Group))            #----number of individuals in each group
  N = sum(n)                             #----number of individuals
  K = length(nBlock)                     #----number of blocks
  if(is.null(Block.name)) {Block.name=paste("Block", 1:K, sep="")}
  P = nBlock
  PP = sum(P)
  
  #----------- split Data into K blocks
  var.block = 0
  K.Data = list()
  for (k in 1:K) {
    K.Data[[k]] = Data[, (var.block + 1):(var.block + nBlock[k])]
    var.block = var.block + nBlock[k]   # variables in blocks
    
    if(norm == TRUE){
      K.Data[[k]]= scale(K.Data[[k]], center=TRUE, scale=FALSE)
      K.Data[[k]] = normM(K.Data[[k]])
    }
    
    if(ScaleDataB == TRUE){
      K.Data[[k]] = scale(K.Data[[k]], center=FALSE, scale=TRUE)
    }
    
  } 
  
  
  #----------- split each block to M groups of individuals
  for (k in 1:K) {
    K.Data[[k]] = transfer_Group(K.Data[[k]], Group)
    for (m in 1:M){
      K.Data[[k]][[m]] = scale(K.Data[[k]][[m]], center=TRUE, scale=ScaleGroup)
    }
  }
  
  #----------------------- concatenated dataset by row as groups for each block
  concat.Data = list()
  for(k in 1:K){ 
    concat.Data [[k]] = K.Data[[k]][[1]]
    if(M>1){
      for(m in 2:M){
        concat.Data[[k]] = rbind(concat.Data[[k]], K.Data[[k]][[m]])
      }
    }
    
    
    
    if(ScaleDataA== TRUE){
      concat.Data[[k]] = scale(concat.Data[[k]], center=FALSE,scale=TRUE)
      K.Data[[k]]      = transfer_Group( concat.Data[[k]], Group)
    }
    
    
  }
  
  
  concat.block.Data = concat.Data[[1]]  #concat.block.Data=[X1|...|XK]
  if(K>1){
    for(k in 2:K){
      concat.block.Data = cbind(concat.block.Data, concat.Data[[k]])
    }
  }
  
  res=list()
  res$concat.block.Data=concat.block.Data
  res$K.Data= K.Data
  res$concat.Data = concat.Data
  return(res)
}
