#Function: ghap.pca
#License: GPLv3 or later
#Modification date: 2 Feb 2016
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute principal components from a kinship matrix

ghap.pca<-function(haplo, K, npc=2){
  
  #Check if haplo is a GHap.haplo object
  if (class(haplo) != "GHap.haplo") {
    stop("Argument haplo must be a GHap.haplo object.")
  }
  
  #Check if kinship matrix is symmetrical
  if(identical(colnames(K),rownames(K)) == FALSE){
    stop("Names in rows and columns must be identical.")
  }
  
  #Check if names in the kinship matrix match with the GHap.haplo object
  if (length(which(colnames(K) %in% haplo$id)) != ncol(K)) {
    stop("All ids in the kinship matrix must be present in the GHap.haplo object.")
  }else{
    ids <- rep(NA, times = ncol(K))
    for (i in 1:length(ids)) {
      ids[i] <- which(haplo$id == colnames(K)[i])
    }
    pop <- haplo$pop[ids]
  }
  
  #Singular value decomposition
  svdK <- svd(K)
  
  #Eigenvalues
  eigenval <- svdK$d[1:npc]
  propvar <- eigenval/sum(eigenval)
  
  #Output
  results <- NULL
  results$eigenvec <- data.frame(pop,rownames(K),svdK$u[,1:npc],stringsAsFactors = FALSE)
  colnames(results$eigenvec) <- c("POP","ID",paste("PC",1:npc,sep=""))
  results$eigenval <- eigenval[1:npc]
  results$propvar <- propvar[1:npc]
  
  #Return output
  return(results)
  
}