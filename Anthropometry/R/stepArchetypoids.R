stepArchetypoids <- function(numArchoid,nearest="cand_ns",data,ArchObj){
  
  N = dim(data)[1]
  
  ai <- archetypes::bestModel(ArchObj[[numArchoid]])
 
  if(is.null(archetypes::parameters(ai))){
    stop("No archetypes computed")  
  }else{
    ras <- rbind(archetypes::parameters(ai),data)
    dras <- dist(ras, method = "euclidean", diag = FALSE, upper = TRUE, p = 2)
    mdras <- as.matrix(dras)
    diag(mdras) = 1e+11
  }
  
  if(nearest == "cand_ns"){
    ini_arch <- sapply(seq(length=numArchoid),nearestToArchetypes,numArchoid,mdras) 
    
    if( all(ini_arch > numArchoid) == FALSE){
      k=1
      neig <- knn(data, archetypes::parameters(ai), 1:N, k=k)
      indices1 <- attr(neig, "nn.index")
      ini_arch <- indices1[,k]
      
      while(any(duplicated(ini_arch))){
        k=k+1  
        neig <- knn(data, archetypes::parameters(ai), 1:N, k=k)
        indicesk <- attr(neig, "nn.index")
        
        dupl <- anyDuplicated(indices1[,1])
        ini_arch <- c(indices1[-dupl,1],indicesk[dupl,k])
      }
    }
    
  }else if(nearest == "cand_alpha"){
    ini_arch <- apply(coef(ai, "alphas"), 2, which.max) 
  }else if(nearest == "cand_beta"){
    ini_arch <- c()
    for (j in 1:numArchoid){
      ini_arch[j] <- which.max(ai$betas[j,])
    }
  }else{
    stop("The nearest vector must be cand_ns, cand_alpha or cand_beta")
  }
  
  res <- archetypoids(numArchoid,data,huge=200,step=TRUE,init=ini_arch)
  cat("Done!") 
  return(list(cases = res[[1]], rss = res[[2]], archet_ini = ini_arch, alphas = res[[4]]))
}
