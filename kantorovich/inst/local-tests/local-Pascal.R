context("Pascal - source")

Pascal_Mn <- function(n){
  M <- matrix(0, nrow=n+1, ncol=n+2)
  for(i in 1:(n+1)){
    M[i,][c(i,i+1)] <- 1
  }
  return(M)
}

centralKernels <- function(Mn.fun, N){
  L <- Kernels <- vector("list", N)
  # initialisation
  k <- 0
  M <- Mn.fun(k)
  m <- nrow(M); n <- ncol(M)
  if(m != 1) stop("M0 must have only one row")
  dims0 <-  as.vector(as.bigz(M))
  Kernels[[k+1]] <- matrix(as.character(dims0), dimnames=list(1:n, 1:m))
  for(k in 1:N){
    M <- Mn.fun(k)
    m <- nrow(M); n <- ncol(M)
    S <- apply(M, 2, function(x) which(x!=0))
    dims <- as.vector(dims0%*%M)
    P <- lapply(1:n, function(i){
      as.character(dims0[S[[i]]]*M[S[[i]],i]/dims[i])
    })
    Kernels[[k+1]] <- matrix("0", nrow=n, ncol=m, dimnames=list(1:n,1:m))
    for(i in 1:n){
      Kernels[[k+1]][i,][S[[i]]] <- P[[i]]
    }
    dims0 <- dims
  }
  return(Kernels)
}

N <- 3
ckernels <- centralKernels(Pascal_Mn, N)

RHO <- lapply(ckernels, function(kernel) matrix("", nrow=nrow(kernel), ncol=nrow(kernel)))
RHO[[1]] <- (diag(2) + 1) %% 2
for(k in 1:N){
  diag(RHO[[k+1]]) <- "0"
  K <- nrow(RHO[[k+1]])
  kernel <- ckernels[[k+1]]
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      RHO[[k+1]][i,j] <- RHO[[k+1]][j,i] <-
        as.character(kantorovich(as.bigq(kernel[i,]), as.bigq(kernel[j,]), dist = RHO[[k]]))
    }
  }
}

