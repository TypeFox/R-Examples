"lmoms.cov" <-
function(x, nmom=5, as.pwm=FALSE, showC=FALSE,
            se=c("NA", "lamse", "lmrse", "pwmse"), ...) {
    se <- match.arg(se)
    if(se == "pwmse") as.pwm <- TRUE

    n <- length(x)
    if(n <= nmom) {
       warning("insufficient length of data")
       return(NA)
    }
    x <- sort(x)

    # falling factorial function
    ff <- function(a, b) gamma(b+1)*choose(a,b)
    b <- pwm(x, nmom=nmom)$betas

    C <- matrix(ncol=nmom, nrow=nmom)
    for(r in 1:nmom) {
       rr <- r-1
       row <- sapply(0:rr, function(k) (-1)^(rr-k)*choose(rr,k)*choose(rr+k,k) )
       if(length(row) > 0) row <- c(row, rep(0, nmom-length(row)))
       C[r,] <- row
    }
    names <- sapply(1:nmom, function(k) paste0("beta2lam",k))
    colnames(C) <- rownames(C) <- names
    if(showC) print(C)

    XX <- x %*% t(x) # matrix to precompute the x[i]*x[j] term for speed enhancement
    # without cost of readability of the code

    ThetaHAT <- matrix(ncol=nmom, nrow=nmom)
    # Slow algorithm but only requires eight lines of code and mimics the symbology
    # used in eq. 22 of Elamir, E.A.H., and Seheult, A.H., 2004, Exact variance
    # structure of sample L-moments: Journal of Statistical Planning and Inference,
    # v. 124, pp. 337--359.
    for(k in 0:(nmom-1)) { # the var-covar structure of the PWMs
       ThetaHAT[k+1,(k+1):nmom] <-
          sapply(k:(nmom-1), function(l) {
             b[k+1]*b[l+1] - sum(sapply(1:(n-1), function(i) {
                                sum(sapply((i+1):n, function(j) {
         return((ff(i-1, k)*ff(j-k-2,l) + ff(i-1, l)*ff(j-l-2,k)) * XX[i,j])
                                }))   }))/ff(n, k+l+2)   })
       ThetaHAT[,k+1] <- ThetaHAT[k+1,]
    }

    if(as.pwm) {
       names <- sapply(0:(nmom-1), function(k) paste0("beta",k))
       colnames(ThetaHAT) <- rownames(ThetaHAT) <- names
       if(se == "pwmse") {
          return(diag(ThetaHAT))
       } else {
          return(ThetaHAT)
       }
    } else {
       LamHAT <- C %*% ThetaHAT %*% t(C)
       names <- sapply(1:nmom, function(k) paste0("lam",k))
       colnames(LamHAT) <- rownames(LamHAT) <- names
       if(se == "lamse") { # standard errors of the lambdas desired
          zz <- sqrt(diag(LamHAT))
          names(zz) <- sapply(1:nmom, function(r) paste0("se.lam",r))
          return(zz)
       } else if(se == "lmrse") { # standard errors of the lambdas and ratios desired
          if(nmom <= 2) { # just the lambdas available
             zz <- sqrt(diag(LamHAT))
             names(zz) <- sapply(1:nmom, function(r) paste0("se.lam",r))
             expectV <- lmr$lambdas[2]; varV <- LamHAT[2,2] # Elamir and Seheult notation
             expectU <- lmr$lambdas[1]; varU <- LamHAT[1,1] # Elamir and Seheult notation
             LCVse <- sqrt((varV/expectV^2 + varU/expectU^2 -
                            2*LamHAT[1,2]/(expectU*expectV))*(expectU/expectV)^2)
             attr(zz, "se.lcv") <- LCVse
             return(zz)
          }
          lmr <- lmoms(x, nmom=nmom) # need the lambda expectations
          expectV <- lmr$lambdas[2]; varV <- LamHAT[2,2] # Elamir and Seheult notation
          B <- varV/expectV^2 # a constant
          lmom.ratio.variances <- sapply(3:nmom, function(r) {
              expectU <- lmr$lambdas[r]; varU <- LamHAT[r,r]
              A <- varU/expectU^2; C <- 2*LamHAT[2,r]/(expectU*expectV)
              D <- (expectU/expectV)^2; return((A+B-C)*D)  })
          zz <- sqrt(c(diag(LamHAT)[1:2], lmom.ratio.variances)) # standard errors
          expectU <- lmr$lambdas[1]; varU <- LamHAT[1,1] # Elamir and Seheult notation
          LCVse <- sqrt((varV/expectV^2 + varU/expectU^2 -
                         2*LamHAT[1,2]/(expectU*expectV))*(expectU/expectV)^2)
          names(zz) <- c(sapply(1:2,    function(r) paste0("se.lam",r)),
                         sapply(3:nmom, function(r) paste0("se.lmr",r)))
          attr(zz, "se.lcv") <- sqrt(LCVse)
          return(zz)
       } else {
          return(LamHAT)
       }
    }
}

