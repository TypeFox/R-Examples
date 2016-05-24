rotateKurtosis <- function(evec, KT, i=1, j=1, k=1, l=1) {
  
  if ((i < 1) | (i > 3)) stop("rotateKurtosis: index i out of range")
  if ((j < 1) | (j > 3)) stop("rotateKurtosis: index j out of range")
  if ((k < 1) | (k > 3)) stop("rotateKurtosis: index k out of range")  
  if ((l < 1) | (l > 3)) stop("rotateKurtosis: index l out of range")

  if (length(dim(evec)) != 3) stop("rotateKurtosis: dimension of direction array is not 3")
  if (dim(evec)[1] != 3) stop("rotateKurtosis: length of direction vector is not 3")
  if (dim(evec)[2] != 3) stop("rotateKurtosis: number of direction vectors is not 3")

  if (length(dim(KT)) != 2) stop("rotateKurtosis: dimension of kurtosis array is not 2")
  if (dim(KT)[1] != 15) stop("rotateKurtosis: kurtosis tensor does not have 15 elements")
  if (dim(KT)[2] != dim(evec)[3]) stop("rotateKurtosis: number of direction vectors does not match number of kurtosis tensors")
  
  nvox <- dim(KT)[2]
  
  ## we create the full symmetric kurtosis tensor for each voxel ...
  W <- defineKurtosisTensor(KT)
  ## ..., i.e., dim(W) == c( 3, 3, 3, 3, nvox) 
  
  Wtilde <- rep(0, nvox)
  
  for (ii in 1:3) {
    for (jj in 1:3) {
      for (kk in 1:3) {
        for (ll in 1:3) {
          
          Wtilde <- Wtilde + evec[ii, i, ] * evec[jj, j, ] * evec[kk, k, ] * evec[ll, l, ] * W[ii, jj, kk, ll, ]   
          
        }
      }
    }
  }
  
  invisible(Wtilde)
  
}


defineKurtosisTensor <- function(DK) {
  
  W <- array(0, dim = c(3, 3, 3, 3, dim(DK)[2]))
  
  W[1, 1, 1, 1, ] <- DK[1, ]
  W[2, 2, 2, 2, ] <- DK[2, ]
  W[3, 3, 3, 3, ] <- DK[3, ]
  
  W[1, 1, 1, 2, ] <- DK[4, ]
  W[1, 1, 2, 1, ] <- DK[4, ]
  W[1, 2, 1, 1, ] <- DK[4, ]
  W[2, 1, 1, 1, ] <- DK[4, ]
  
  W[1, 1, 1, 3, ] <- DK[5, ]
  W[1, 1, 3, 1, ] <- DK[5, ]
  W[1, 3, 1, 1, ] <- DK[5, ]
  W[3, 1, 1, 1, ] <- DK[5, ]
  
  W[2, 2, 2, 1, ] <- DK[6, ]
  W[2, 2, 1, 2, ] <- DK[6, ]
  W[2, 1, 2, 2, ] <- DK[6, ]
  W[1, 2, 2, 2, ] <- DK[6, ]
  
  W[2, 2, 2, 3, ] <- DK[7, ]
  W[2, 2, 3, 2, ] <- DK[7, ]
  W[2, 3, 2, 2, ] <- DK[7, ]
  W[3, 2, 2, 2, ] <- DK[7, ]
  
  W[3, 3, 3, 1, ] <- DK[8, ]
  W[3, 3, 1, 3, ] <- DK[8, ]
  W[3, 1, 3, 3, ] <- DK[8, ]
  W[1, 3, 3, 3, ] <- DK[8, ]
  
  W[3, 3, 3, 2, ] <- DK[9, ]
  W[3, 3, 2, 3, ] <- DK[9, ]
  W[3, 2, 3, 3, ] <- DK[9, ]
  W[2, 3, 3, 3, ] <- DK[9, ]
  
  W[1, 1, 2, 2, ] <- DK[10, ]
  W[1, 2, 1, 2, ] <- DK[10, ]
  W[1, 2, 2, 1, ] <- DK[10, ]
  W[2, 1, 2, 1, ] <- DK[10, ]
  W[2, 1, 1, 2, ] <- DK[10, ]
  W[2, 2, 1, 1, ] <- DK[10, ]
  
  W[1, 1, 3, 3, ] <- DK[11, ]
  W[1, 3, 1, 3, ] <- DK[11, ]
  W[1, 3, 3, 1, ] <- DK[11, ]
  W[3, 1, 3, 1, ] <- DK[11, ]
  W[3, 1, 1, 3, ] <- DK[11, ]
  W[3, 3, 1, 1, ] <- DK[11, ]
  
  W[2, 2, 3, 3, ] <- DK[12, ]
  W[2, 3, 2, 3, ] <- DK[12, ]
  W[2, 3, 3, 2, ] <- DK[12, ]
  W[3, 2, 3, 2, ] <- DK[12, ]
  W[3, 2, 2, 3, ] <- DK[12, ]
  W[3, 3, 2, 2, ] <- DK[12, ]
  
  W[1, 1, 2, 3, ] <- DK[13, ]
  W[1, 1, 3, 2, ] <- DK[13, ]
  W[3, 1, 1, 2, ] <- DK[13, ]
  W[2, 1, 1, 3, ] <- DK[13, ]
  W[2, 3, 1, 1, ] <- DK[13, ]
  W[3, 2, 1, 1, ] <- DK[13, ]
  W[1, 3, 1, 2, ] <- DK[13, ]
  W[1, 2, 1, 3, ] <- DK[13, ]
  W[3, 1, 2, 1, ] <- DK[13, ]
  W[2, 1, 3, 1, ] <- DK[13, ]
  W[1, 2, 3, 1, ] <- DK[13, ]
  W[1, 3, 2, 1, ] <- DK[13, ]
  
  W[2, 2, 1, 3, ] <- DK[14, ]
  W[2, 2, 3, 1, ] <- DK[14, ]
  W[3, 2, 2, 1, ] <- DK[14, ]
  W[1, 2, 2, 3, ] <- DK[14, ]
  W[1, 3, 2, 2, ] <- DK[14, ]
  W[3, 1, 2, 2, ] <- DK[14, ]
  W[3, 2, 1, 2, ] <- DK[14, ]
  W[2, 3, 2, 1, ] <- DK[14, ]
  W[1, 2, 3, 2, ] <- DK[14, ]
  W[2, 1, 2, 3, ] <- DK[14, ]
  W[2, 1, 3, 2, ] <- DK[14, ]
  W[2, 3, 1, 2, ] <- DK[14, ]
  
  W[3, 3, 1, 2, ] <- DK[15, ]
  W[3, 3, 2, 1, ] <- DK[15, ]
  W[2, 3, 3, 1, ] <- DK[15, ]
  W[1, 3, 3, 2, ] <- DK[15, ]
  W[1, 2, 3, 3, ] <- DK[15, ]
  W[2, 1, 3, 3, ] <- DK[15, ]
  W[1, 3, 2, 3, ] <- DK[15, ]
  W[2, 3, 1, 3, ] <- DK[15, ]
  W[3, 2, 3, 1, ] <- DK[15, ]
  W[3, 1, 3, 2, ] <- DK[15, ]
  W[3, 1, 2, 3, ] <- DK[15, ]
  W[3, 2, 1, 3, ] <- DK[15, ]
  
  invisible( W)
}

## TODO: For efficiency we might combine both functions
##       Many doubled calculations


kurtosisFunctionF1<- function(l1, l2, l3) {
  
  #require(gsl)
  ## Tabesh et al. Eq. [28], [A9, A12]
  
  ## consider removable singularities!! 
  F1 <- numeric(length(l1))
  
  ind12 <- abs(l1 - l2) > 1e-10 ## l1 != l2
  ind13 <- abs(l1 - l3) > 1e-10 ## l1 != l3
  
  ind <- (ind12) & (ind13)
  l1i <- l1[ind]
  l2i <- l2[ind]
  l3i <- l3[ind]
  if (any(ind)) F1[ind] <- 
    (l1i+l2i+l3i)^2 /18 /(l1i-l2i)/(l1i-l3i) *
    (ellint_RF(l1i/l2i, l1i/l3i, 1) * sqrt(l2i*l3i)/l1i + 
        ellint_RD(l1i/l2i, l1i/l3i, 1) * (3*l1i^2 - l1i*l2i - l1i*l3i - l2i*l3i)/(3*l1i*sqrt(l2i*l3i)) 
      - 1) 
  ind1 <- (!ind12) & (ind13)
  if (any(ind1)) F1[ind1] <- kurtosisFunctionF2(l3[ind1], l1[ind1], l1[ind1]) / 2 
  
  ind2 <- (ind12) & (!ind13)
  if (any(ind2)) F1[ind2] <- kurtosisFunctionF2(l2[ind2], l1[ind2], l1[ind2]) / 2 
  
  ##  at singularities ind3
  ##  we have ellint_RF(1, 1, 1) = 1, ellint_RD(1, 1, 1) 
  ##  Result should be Inf ???
  ind3 <- (!ind12) & (!ind13)
  if (any(ind3)) F1[ind3] <- 1/5
  
  F1
}

kurtosisFunctionF2 <- function(l1, l2, l3) {
  
  #require(gsl)
  ## Tabesh et al. Eq. [28], [A10, A11, A12]
  
  alpha <- function(x) {
    z <- rep(1, length(x))
    z[x > 0] <- 1/sqrt(x[x > 0]) * atanh(sqrt(x[x > 0]))
    z[x < 0] <- 1/sqrt(-x[x < 0]) * atan(sqrt(-x[x < 0]))
    z
  }
  
  ## consider removable singularities!!
  F2 <- numeric(length(l1))
  
  ind23 <- abs(l2 - l3) > 1e-10 ## l2 != l3
  ind12 <- abs(l1 - l2) > 1e-10 ## l1 != l2
  if (any(ind23)){
    l1i <- l1[ind23]
    l2i <- l2[ind23]
    l3i <- l3[ind23]  
    F2[ind23] <- 
      (l1i + l2i + l3i)^2 / (l2i - l3i)^2 / 3 * 
      (ellint_RF(l1i/l2i, l1i/l3i, 1) * (l2i + l3i) / sqrt(l2i*l3i) + 
          ellint_RD(l1i/l2i, l1i/l3i, 1) * (2*l1i - l2i - l3i) / 3 / sqrt(l2i*l3i) 
        - 2 ) 
  }
  ind1 <- (!ind23) & (ind12)
  if (any(ind1)) {
    l1i <- l1[ind1]
    l2i <- l2[ind1]
    l3i <- l3[ind1]  
    F2[ind1] <- 
      6 * (l1i + 2*l3i)^2 / 144 / l3i^2 / (l1i - l3i)^2 *
      (l3i * (l1i + 2*l3i) + l1i * (l1i - 4 * l3i) * alpha(1 - l1i / l3i))
  }
  ind2 <- (!ind23) & (!ind12)
  if (any(ind2)) F2[ind2] <- 6/15
  
  F2
}

pseudoinverseSVD <- function(xxx, eps=1e-8) {
  svdresult <- svd(xxx)
  d <-  svdresult$d
  dinv <- 1/d
  dinv[abs(d) < eps] <- 0
  svdresult$v %*% diag(dinv) %*% t(svdresult$u)
}
