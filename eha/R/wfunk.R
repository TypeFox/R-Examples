wfunk <- function(beta = NULL, lambda, p, X = NULL, Y,
                  offset = rep(0, length(Y)),
                  ord = 2, pfixed = FALSE){
    
    ## Returns loglik, score, and information (=-fpp)
    ## For one stratum (only)!!
    
    if (ord < 0) return(NULL)
    nn <- NROW(Y)
    if (NCOL(Y) == 2) Y <- cbind(rep(0, nn), Y)
    
    if (is.null(X)){
        if (pfixed){
            bdim <- 1
            b <- -log(lambda)
        }else{
            bdim <- 2
            b <- c(-log(lambda), log(p))
        }
        
        mb <- 0
        
        fit <- .Fortran("wfuncnull",
                        as.integer(ord),
                        as.integer(pfixed),
                        as.double(p),
                        as.integer(bdim),
                        as.double(b),
                        as.integer(nn),
                        as.double(Y[, 1]),
                        as.double(Y[, 2]),
                        as.integer(Y[, 3]),
                                        #
                        f = double(1),
                        fp = double(bdim),
                        fpp = double(bdim * bdim),
                        ok = integer(1),
                        ## DUP = FALSE,
                        PACKAGE = "eha"
                        )
    }else{
        mb <- NCOL(X)
        if (length(beta) != mb) stop("beta mis-specified!")
        if (pfixed){
            bdim <- mb + 1
            b <- c(beta, -log(lambda))
        }else{
            bdim <- mb + 2
            b <- c(beta, -log(lambda), log(p))
        }
        
        fit <- .Fortran("wfunc", ## Returns -loglik, -score, +information
                        as.integer(ord),
                        as.integer(pfixed),
                        as.double(p),
                        as.integer(bdim),
                        as.integer(mb),
                        as.double(b),
                                        #
                        as.integer(nn),
                        as.double(t(X)),
                        as.double(Y[, 1]), ## enter
                        as.double(Y[, 2]), ## exit
                        as.integer(Y[, 3]), ## event
                        as.double(offset),
                                        #
                        f = double(1),
                        fp = double(bdim),
                        fpp = double(bdim * bdim),
                        ok = integer(1),
                        ## DUP = FALSE,
                        PACKAGE = "eha"
                        )
    }
    
    ret <- list(f = -fit$f)
    if (ord >= 1){
        ## The delta method: (log(lambda) <--> -log(lambda):
        xx <- rep(1, bdim)
        xx[mb + 1] <- -1
        ret$fp <- -xx * fit$fp
        if (ord >= 2){
            xx <- diag(xx)
            ret$fpp <- xx %*% matrix(fit$fpp, ncol = bdim) %*% t(xx)
        }
    }    
    ret
}
