
context("sparseMVN")

test_that("sparseMVN", {

    require(mvtnorm)
    set.seed(123)
    N <- 10 ## number of samples
    m <- 15  ## number of blocks in sparse covariance matrix
    p <- 3 ## size of each block
    k <- 7  ## 
    
    ## build block-arrow covariance/precision matrix for test

    mu <- seq(-3,3,length=p*m+k)
    Q1 <- tril(kronecker(Matrix(seq(0.1,p,length=p*p),p,p),diag(m)))
    Q2 <- cBind(Q1,Matrix(0,m*p,k))
    Q3 <- rBind(Q2,cBind(Matrix(rnorm(k*m*p),k,m*p),Diagonal(k)))
    CV <- Matrix::tcrossprod(Q3)
    chol.CV <- Matrix::Cholesky(CV)  ## creates a dCHMsimpl object
    
    ## convert to base R for comparison

    CV.base <- as(CV, "matrix")

    ## sample N MVN's

    for (prec in c(FALSE, TRUE)) {
        x.sp <- rmvn.sparse(N, mu, chol.CV,prec=prec)
        
        ## check dimensions
        
        ## each row is a draw
        ## each col is a variable
        expect_equal(NROW(x.sp), N)  
        expect_equal(NCOL(x.sp), length(mu))  
        
        ## computing log densities using dmvn.sparse
        d.sp <- dmvn.sparse(x.sp, mu, chol.CV, prec=prec)

        if (prec) {
            mat <- solve(CV.base)
        } else {
            mat <- CV.base
        }
         
        ## computing log densities using dmvnorm
        d.dns <- dmvnorm(x.sp, mu, mat, log=TRUE)
        expect_equal(d.sp, d.dns)
        }
})
          
