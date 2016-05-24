SortvarLearn <- function(data,
                         knownlabels,
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
  
  
  # check whether the knownlabels is missing
  if ( missing(knownlabels)){
    stop("labels are missing!")
  }
  
  
  if(missing(knownlabels) || length(knownlabels)==0){
    warning("knownlabels are missing, intrumental initialization without output effect")
    knownlabels <- rep(1, nrow(data)) 
    
  }
  
  
  if(min(knownlabels) <= 0 || length(knownlabels) != nrow(data)){
    stop("Each observation in knownLabels must have a valid cluster affectation !")
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
  nbCluster <- as.integer(max(knownlabels))
  
  VarRole <- matrix(NA,(length(lambda)*length(rho)), p) 
  VarRole <- DiscriminantAnalysisGlasso(data, nbCluster, lambda, rho, knownlabels = knownlabels, nbCores)
  var.role.sum <- colSums(VarRole) 
  OrderVariable <- sort.int(var.role.sum,decreasing=TRUE,index.return=TRUE)$ix
  
  return(OrderVariable)    
}
