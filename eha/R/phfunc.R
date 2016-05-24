phfunc <- function(beta = NULL, lambda, p, X = NULL, Y,
                  offset = rep(0, length(Y)),
                  ord = 2, pfixed = FALSE, dist = "weibull"){

## Returns loglik, score, and information (=-fpp)
## For one stratum (only)!!

  if (ord < 0) return(NULL)

  if (dist == "weibull"){
    dis <- 0
  }else if(dist == "loglogistic"){
    dis <- 1
  }else if (dist == "lognormal"){
    dis <- 2
  }


    nn <- NROW(Y)
    if (NCOL(Y) == 2) Y <- cbind(rep(0, nn), Y)

    if (is.null(X)){
        stop("No X not implemented yet")
        if (FALSE){    #### NOT yet!!
            if (pfixed){
                bdim <- 1
                b <- log(lambda)
            }else{
                bdim <- 2
                b <- c(log(lambda), log(p))
            }

            fit <- .Fortran("phfuncnull",
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
        }

    }else{
        mb <- NCOL(X)
        if (length(beta) != mb) stop("beta mis-specified!")
        if (pfixed){
          stop("p fixed not implemented yet")
            bdim <- mb + 1
            b <- c(beta, log(lambda))
        }else{
            bdim <- mb + 2
            b <- c(beta, log(lambda), log(p))
        }
        b <- beta
        alpha <- log(lambda)
        gamma <- log(p)

        d0 <- .C("loglik_ph",
                 as.integer(dis),
                 as.integer(mb),
                 as.double(b),
                 as.double(alpha),
                 as.double(gamma),
                 as.integer(nn),
                 as.double(t(X)),
                 as.double(Y[, 1]),
                 as.double(Y[, 2]),
                 as.integer(Y[, 3]),
                 as.double(offset),
                 f = double(1), ## Return value
                 ## DUP = FALSE,
                 PACKAGE = "eha")
        ret <- list(f = - d0$f)
        if (ord >= 1){
          d1 <- .C("d_loglik_ph",
                   as.integer(dis),
                   as.integer(mb),
                   as.double(b),
                   as.double(alpha),
                   as.double(gamma),
                   as.integer(nn),
                   as.double(t(X)),
                   as.double(Y[, 1]),
                   as.double(Y[, 2]),
                   as.integer(Y[, 3]),
                   as.double(offset),
                   fp = double(bdim), ## Return value
                   ## DUP = FALSE,
                   PACKAGE = "eha")
          ret$fp <- d1$fp
          if (ord >= 2){
            d2 <- .C("d2_loglik_ph",
                     as.integer(dis),
                     as.integer(mb),
                     as.double(b),
                     as.double(alpha),
                     as.double(gamma),
                     as.integer(nn),
                     as.double(t(X)),
                     as.double(Y[, 1]),
                     as.double(Y[, 2]),
                     as.integer(Y[, 3]),
                     as.double(offset),
                     fpp = double(bdim * bdim), ## Return value
                     ## DUP = FALSE,
                     PACKAGE = "eha")
            ret$fpp <- matrix(d2$fpp, ncol = bdim)
          }
        }
    }

  return(ret)
}

