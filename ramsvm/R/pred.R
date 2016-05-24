pred_linear <- function(x.test, 
                        beta, 
                        beta0, 
                        k) {

  XI <- XI.gen(k = k, kd = as.double(k))

  beta0 <- matrix(beta0, 
                  nrow = nrow(x.test),  
                  ncol = ncol(beta), 
                  byrow=TRUE)

  f.matrix <- t(x.test %*% beta + beta0)

  inner.matrix <- matrix(data = 0.0, nrow = nrow(x.test), ncol = k)

  for( ii in 1L:k ) inner.matrix[,ii] <- colSums(f.matrix * XI[,ii])

  z <- apply(X = inner.matrix, MARGIN = 1L, FUN = pred)

  return(z)

}

pred_kernel <- function(x.test, 
                        x,
                        kernel,
                        kparam,
                        beta,  
                        beta0,  
                        k) {

  XI <- XI.gen(k = k, kd = as.double(k))

  K.test <- Kmat(x = x.test, 
                 y = x, 
                 kernel = kernel, 
                 kparam = kparam)

  beta0 <- matrix(beta0, 
                  nrow = nrow(x.test),  
                  ncol = ncol(beta), 
                  byrow=TRUE)

  f.matrix <- t(K.test %*% beta + beta0)

  inner.matrix <- matrix(data = 0.0, nrow = nrow(x.test), ncol = k)

  for( ii in 1L:k ) inner.matrix[,ii] <- colSums(f.matrix * XI[,ii])

  z <- apply(X = inner.matrix, MARGIN = 1L, FUN = pred)

  return(z)

}

pred <- function(f) {
  tst <- sapply(f, function(i){isTRUE(all.equal(i,max(f)))})
  y <- min(which(tst))
  return(y)
}
