
hyavar <- function(yuima, bw, nonneg = TRUE, psd = TRUE){
  
  x <- get.zoo.data(yuima)
  n.series <- length(x)
  
  # allocate memory
  ser.X <- vector(n.series, mode="list")     # data in 'x'
  ser.times <- vector(n.series, mode="list") # time index in 'x'
  ser.diffX <- vector(n.series, mode="list")
  ser.num <- integer(n.series) # number of observations
  avar.cov <- diag(n.series) # asymptotic variances of cce(yuima)$covmat
  avar.cor <- diag(0, n.series) # asymptotic variances of cce(yuima)$cormat
  
  for(i in 1:n.series){
    
    # NA data must be skipped
    idt <- which(is.na(x[[i]]))
    if(length(idt>0)){
      x[[i]] <- x[[i]][-idt]
    }
    if(length(x[[i]])<2) {
      stop("length of data (w/o NA) must be more than 1")
    }
    
    # set data and time index
    ser.X[[i]] <- as.numeric(x[[i]]) # we need to transform data into numeric to avoid problem with duplicated indexes below
    ser.times[[i]] <- as.numeric(time(x[[i]]))
    ser.diffX[[i]] <- diff(ser.X[[i]])
    # set number of observations
    ser.num[i] <- length(ser.X[[i]])
    
  }
  
  # Computing the Hayashi-Yoshida estimator
  result <- cce(yuima, psd = psd)
  cmat <- result$covmat
  
  if(n.series > 1){
    
    if(missing(bw)){
      
      bw <- matrix(0, n.series, n.series)
      
      for(i in 1:(n.series - 1)){
        
        for(j in (i + 1):n.series){
          
          Init <- min(ser.times[[i]][1], ser.times[[j]][1])
          Term <- max(ser.times[[i]][ser.num[i]], ser.times[[j]][ser.num[j]])
          bw[i, j] <- (Term - Init) * min(ser.num[c(i,j)])^(-0.45)
          bw[j, i] <- bw[i, j]
          
        }
        
      }
      
    }
    
  }
  
  bw <- matrix(bw, n.series, n.series)
  
  for(i in 1:n.series){
    
    avar.cov[i, i] <- (2/3) * sum(ser.diffX[[i]]^4)
    
    for(j in 1:i){
      
      if(i > j){
        
        N.max <- max(ser.num[i],ser.num[j])
        
        avar <- .C("hyavar",
                   as.double(ser.times[[i]]),
                   as.double(ser.times[[j]]),
                   as.integer(ser.num[i]),
                   as.integer(ser.num[j]),
                   as.double(ser.X[[i]]),
                   as.double(ser.X[[j]]),
                   as.integer(N.max),
                   as.double(bw[i, j]),
                   avar = double(4))$avar
        ## avar[1] = var(HYij), avar[2] = cov(HYii, HYij), avar[3] = cov(HYij, HYjj), avar[4] = cov(HYii, HYjj)
        
        avar.cov[i, j] <- avar[1]
        avar.cov[j, i] <- avar[1]
        
        Pi <- matrix(c(avar.cov[i,i], avar[2], avar[4], 
                       avar[2], avar[1], avar[3], 
                       avar[4], avar[3], avar.cov[j,j]), 3, 3)
        
        if(nonneg){ # Taking the spectral absolute value of Pi
          tmp <- eigen(Pi, symmetric = TRUE)
          Pi <- tmp$vectors %*% diag(abs(tmp$values)) %*% t(tmp$vectors)
        }
        
        d <- c(-0.5 * cmat[i,j]/cmat[i,i], 1, -0.5 * cmat[i,j]/cmat[j,j])
        
        avar.cor[i, j] <- d %*% Pi %*% d / (cmat[i,i] * cmat[j,j])
        avar.cor[j, i] <- avar.cor[i, j]
        
      }
      
    }
    
  }
  
  return(list(covmat = cmat, cormat = result$cormat, avar.cov = avar.cov, avar.cor = avar.cor))
}

