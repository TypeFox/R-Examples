SortvarClust <- function(data,
                         nbCluster,
                         lambda,
                         rho,
                         nbCores)
{
  # check data parameter
  if(missing(data)){
    stop("data is missing !")
  } 
  if(is.matrix(data) == FALSE && is.data.frame(data) == FALSE) 
    stop(paste(sQuote("data"), "must be a matrix"))
  
  
  # check nbCluster parameter
  if(missing(nbCluster)){
    stop("nbCluster is missing!")
  }
  if(sum(!is.wholenumber(nbCluster))){
    stop("nbCluster must contain only integer!")
  }
  if(sum(nbCluster < 1)){ 
    stop(paste(sQuote("nbCluster"), "must be an integer greater than 0!"))
  }
  
  
  # check lambda parameter
  if(missing(lambda)){
    stop("lambda is missing!")
  }
  if(is.vector(lambda) == FALSE || length(lambda) <= 1){ 
    stop(paste(sQuote("lambda"), "must be a vector with length >= 2"))
  }
  if (sum(lambda<=0)){
    stop("lambda must be greater than 0!")
  }
  
  
  # check rho parameter
  if(missing(rho)){
    stop("rho is missing!")
  }
  if(is.vector(rho) == FALSE){ 
    stop(paste(sQuote("rho"), "must be a vector"))
  }
  if(sum(rho<=0)){
    stop("rho must be greater than 0!")
  }
 
  # check nbCores 
  nb.cpus <- detectCores(all.tests = FALSE, logical = FALSE)
  if(missing(nbCores))
  {
    if(nb.cpus > 1)
      nbCores <- 2
    if(nb.cpus == 1)
      nbCores <- 1
  }
  
  
  data <- as.matrix(scale(data, TRUE, TRUE))
  n <- as.integer(nrow(data))
  p <- as.integer(ncol(data))
  K <- as.integer(nbCluster)
  
  
  VarRole <- array(NA,dim=c((length(lambda)*length(rho)), p, length(nbCluster))) 
  VarRole <- ClusteringEMGlasso(data,nbCluster,lambda,rho, nbCores)
  ## Calcul de la matrice O de taille length(nCluster) * p
  Matrix0 <- matrix(0, nrow=length(nbCluster), ncol=p)
  for (k in 1:length(nbCluster))
    Matrix0[k,]<- colSums(VarRole[,,k])    
  
  
  ## Ordre des variables
  OrderVariable <- matrix(NA, nrow=length(nbCluster),ncol=p)
  for (k in 1:length(nbCluster))
    OrderVariable[k,] <- sort.int(Matrix0[k,],decreasing=TRUE,index.return=TRUE)$ix
  
  return(OrderVariable)    
}
