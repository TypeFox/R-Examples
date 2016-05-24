#Function: ghap.simpheno
#License: GPLv3 or later
#Modification date: 2 Feb 2016
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Simulate phenotypes

ghap.simpheno<-function(
  haplo,
  K,
  h2,
  g2,
  major=NULL,
  seed=NULL
){
  
  #Check if haplo is a GHap.haplo object
  if(class(haplo) != "GHap.haplo"){
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
    ids <- rep(NA, times = nrow(K))
    for (i in 1:length(ids)) {
      ids[i] <- which(haplo$id == rownames(K)[i])
    }
  }

  
  # Simulate major haplotypes
  if(is.null(major) == FALSE){
    if(g2 < 0 | g2 > 1){
      stop("Argument g2 must be between zero and one.")
    }
    X <- as.matrix(haplo$genotypes[major,ids])
    if(length(major) > 1){
      X <- t(X)
      X <- scale(X, center = TRUE, scale = FALSE)
    }else{
      X <- scale(X, center = TRUE, scale = FALSE)
    }
    if(is.null(seed) == FALSE){
      set.seed(seed)
    }
    b <- rnorm(length(major),mean=0,sd=1)
    Xb <-  X%*%b
    varg <- var(Xb)*(1-g2)/g2
  }else{
    b <- NA
    Xb <- 0
    varg <- 1
  }
  
  # Simulate uncorrelated random effects
  if(is.null(seed) == FALSE){
    set.seed(seed)
  }
  g <- rnorm(haplo$nsamples.in, mean=0, sd=sqrt(varg))
  
  # Cholesky decomposition of relationship matrix
  svdK <- svd(K)
  L <- svdK$u %*% diag(sqrt(svdK$d))
  
  # Make random effects correlated by K
  # Now g = polygenic effect (sum of genome-wide haplotype effects) - major effect
  g <- as.vector(L%*%g)
  
  # Simulate breeding value
  u <-  Xb + g
  varu <- var(u)
  
  # Simulate residuals
  vare <- varu*((1-h2)/h2)
  if(is.null(seed) == FALSE){
    set.seed(seed)
  }
  e <- rnorm(haplo$nsamples.in, 0, sqrt(vare))
  
  # Simulate phenotypes
  y <- u + e
  
  #Return output
  sim <- NULL
  sim$h2 <- h2
  sim$g2 <- g2
  sim$major <- major
  sim$major.effect <- b
  sim$u <- as.vector(u)
  names(sim$u) <- colnames(K)
  sim$varu <- as.vector(varu)
  sim$vare <- as.vector(vare)
  sim$data <- data.frame(y,colnames(K))
  colnames(sim$data) <- c("phenotype","individual")
  return(sim)
  
}