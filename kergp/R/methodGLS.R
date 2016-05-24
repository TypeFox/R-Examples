##============================================================================
##
## 'gls' is for a covAll with known params. This function returns the
## estimated 'beta.hat' for beta in the model
## 
##  y = F beta + eta (+ epsilon)
##
## where 'eta' is Norm(0, K) with variance K := [ k(x_i, x_j) ] and
## k ins the kernel function in 'object'. The optional term + epsilon
## is for an extra withe noise the variance of which can be ggiven by
## 'varNoise'.
##
## Let 'L' be the root in the "LLprime" Cholesky decomposition
## C = L %*% t(L). With Linv := inv(L) we get an OLS regression by
## transforming the model into 
##
##  (Linv %*% y) = (Linv %*% F) beta + (Linv %*% eta)
##
## and then using classical QR tips in regression, see 
## Kenneth Lange (2010), "Numerical Analysis for Statisticians" 2nd ed.
## pp. 102-103. Springer-Verlag,
##
## The function also returns z := Linv %*%eta.hat. Note that the name 'z'
## is for compatibility with DiceKriging but could be changed.
##
##============================================================================
setMethod("gls",
          signature = signature(object = "covAll"),
          definition = function(object, y, X, F = NULL, varNoise = NULL, 
              beta = NULL, checkNames = TRUE, ...) {
            
              ##cat("gls estimation with \"covTS\" covariance\n")
              n <- length(y)
              if (nrow(X) != n) stop("'X' must be a matrix with length(y) rows")
              
              C <- covMat(object, X, compGrad = FALSE, checkNames = checkNames)
              
              if (!is.null(varNoise)) {
                  noise <- TRUE
                  varNoise <- as.numeric(varNoise)
                  if (is.na(varNoise) || varNoise < 0) {
                      stop("'varNoise' must be a positive numeric")
                  }
                  diag(C) <- diag(C) + varNoise
              } else {
                  noise <- FALSE
              }
              L <- t(chol(C))
              
              if (!is.null(F)) {
                  p <- NCOL(F)
                  p1 <- p + 1L
                  FStarPlus <- forwardsolve(L, cbind(F, y))
                  qrFStarPlus <- qr(FStarPlus)
                  RStarPlus <- qr.R(qrFStarPlus)
                  if (length(beta)>0){
                      betaHat <- beta
                      eStar <- FStarPlus %*% (rbind(-as.matrix(betaHat), 1))
                      sseStar <- crossprod(eStar, eStar)
                  } else {
                      betaHat <- backsolve(RStarPlus[1L:p, 1L:p], RStarPlus[1L:p,  p1])  
                      eStar  <- qr.Q(qrFStarPlus)[ , p1] * RStarPlus[p1, p1]
                      sseStar <- RStarPlus[p1, p1]^2
                  }
                  FStar <- FStarPlus[ , 1L:p, drop = FALSE]
                  RStar <- RStarPlus[1L:p, 1L:p]  
              } else {
                  betaHat <-  NULL
                  eStar <- forwardsolve(L, y)
                  sseStar <- crossprod(eStar, eStar)
                  FStar <- NULL
                  RStar <- NULL
              }
              names(betaHat) <- colnames(F)
              
              return(list(betaHat = betaHat,    ## estimate
                   L = L,                ## Cholesky Lower root of C 
                   eStar = eStar,        ## inv(L) %*% error
                   sseStar = sseStar,    ## t(eStar) %*% eStar
                   FStar = FStar,        ## inv(L) %*% F
                   RStar = RStar,        ## qr.r(FStar)
                   noise = noise,
                   varNoise = varNoise))  ## copy
          }
          )
