if (!file.exists("test-dataset.Rdata")) {

  library(robustvarComp)
  if (require(mvtnorm) & require(Matrix)) {
    ones <- function (n, p = n)
      matrix(rep(1, n * p), nrow = n, ncol = p)

    size <- 20
    i <- 2
    j <- 2
    n <- 3
    p <- i*j*n ## number of columns of Y  
    k <- 5 ## number of beta parameters
    R <- 3 ## number of gamma parameters
    beta <- c(0,rep(2,k)) ## values beta parameters k+1 parameters,
                      ## the first the intercept
    gamma <- c(1,1,2)/4 ## values gamma parameters
    eta0 <- 1/4 ## value of eta_0
            ## recall that the vector eta is eta_0 * gamma  
    V <- array(0, dim=c(p,p,R)) ## V_r matrix/ces 
                             ## (just one in this example) 
                             ## with dimension pXp
    V[,,1] <- V1 <- diag(i)%x%ones(j)%x%ones(n)
    V[,,2] <- V2 <- ones(i)%x%diag(j)%x%ones(n)
    V[,,3] <- V3 <- diag(i)%x%diag(j)%x%ones(n)

    S <- eta0*(robustvarComp:::Vprod(V, gamma)) ## Variance and Covariance matrix of Y
                 ## i.e. \eta_0 (V_0+\sum_{r=1}^{R}\gamma_{0r} V_{r})
                 ## where V_0 = I_0
                 ## V0 is added in Vprod function automatically

    VV1 <- VV2 <- VV3 <- list()
    for (N in 1:size) {
      VV1[[N]] <- V1
      VV2[[N]] <- V2
      VV3[[N]] <- V3              
    }
    v1 <- as.matrix(bdiag(VV1))
    v2 <- as.matrix(bdiag(VV2))
    v3 <- as.matrix(bdiag(VV3))

    set.seed(1234)
    x <- array(0, dim=c(p,size,k+1))
    for (N in 1:size)
      x[,N,] <- cbind(rep(1,p), matrix(rnorm(p*k, mean=0, sd=1), p, k))
    e <- rmvnorm(n=size, mean=rep(0, p), sigma=S)
    y <- t(e)+robustvarComp:::xprod(x, beta) ## the Y observations
  
    Y <- as.vector(y)
    X <- matrix(nrow=0, ncol=k+1)
    for (N in 1:size) {
      X <- rbind(X, x[,N,])
    }
            
    dati <- data.frame(Y=c(Y), X=X)

    test.control <- varComprob.control()
    test.init.control <- varComprob.control(cov.init="covOGK")
  
    save(beta, gamma, eta0, y, x, V, Y, X, v1, v2, v3, n, p, k, R, dati, test.control, test.init.control, file="test-dataset.Rdata")
  }
}
