context("Tests using Rosenbrock function")

test_that("Rosenbrock", {

    f.rosen <- function(V) {
        
        N <- length(V)/2
        x <- V[seq(1,2*N-1,by=2)]
        y <- V[seq(2,2*N,by=2)]
        return(sum(100*(x^2-y)^2+(x-1)^2))
    }
    
    df.rosen <- function(V) {
        N <- length(V)/2
        x <- V[seq(1,2*N-1,by=2)]
        y <- V[seq(2,2*N,by=2)]
        
        t <- x^2-y
        dxi <- 400*t*x+2*(x-1)
        dyi <- -200*t
        return(as.vector(rbind(dxi,dyi)))
    }
    
    hess.rosen <- function(V) {
        
        N <- length(V)/2
        x <- V[seq(1,2*N-1,by=2)]
        y <- V[seq(2,2*N,by=2)]
        d0 <- rep(200,N*2)
        d0[seq(1,(2*N-1),by=2)] <- 1200*x^2-400*y+2
        d1 <- rep(0,2*N-1)
        d1[seq(1,(2*N-1),by=2)] <- -400*x
        
        H <- bandSparse(2*N,
                        k=c(-1,0,1),
                        diagonals=list(d1,d0,d1),
                        symmetric=FALSE,
                        giveCsparse=TRUE)
        return(drop0(H))
    }

    set.seed(123)
    N <- 3
    start <- as.vector(rnorm(2*N,-1,3))

    m <- list(list(hs=hess.rosen, method="Sparse", precond=0),
              list(hs=NULL, method="BFGS", precond=0),
              list(hs=NULL, method="SR1", precond=0),
              list(hs=hess.rosen, method="Sparse", precond=1)
              )

    for (meth in m) {        
        opt0 <- trust.optim(start,
                            fn=f.rosen,
                            gr=df.rosen,
                            hs=meth$hs,
                            method=meth$method,
                            control=list(
                                preconditioner=meth$precond,
                                report.freq=0L,
                                maxit=5000L
                            )
                            )
        
        expect_equal(opt0$fval, 0, tolerance=1e-8)
        expect_equal(opt0$solution, rep(1,2*N), tolerance=1e-8)
        expect_match(opt0$status, "Success")
        expect_match(opt0$method, meth$method)
    }
})
          
          
    
