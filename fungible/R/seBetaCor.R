###################################################################### 
# This function computes the asymptotic covariance matrix of         # 
# standardized regression coefficients using correlations.           # 
#                                                                    #
# Arguments:                                                         # 
#                                                                    # 
#  R       - A p x p predictor correlation matrix.                   # 
#  rxy     - A p x 1 vector of predictor-criterion correlations      #  
#  Nobs    - number of observations.                                 # 
#  alpha   - desired Type I error rate; default = .05.               #  
#  digits  - number of significant digits to print; default = 3.     # 
#  covmat  - A (p+1)p/2 x (p+1)p/2 covariance matrix of correlations.#
#            default = 'normal'. The default option computes an      #
#            asymptotic covariance matrix under the assumption of    #
#            multivariate normal data. Users can supply a covariance #
#            matrix under asymptotic distribution free (ADF) or      #
#            elliptical distributions when available.                #
#                                                                    #  
# Output                                                             # 
#                                                                    # 
#   cov.mat - covariance matrix of standardized regression           # 
#             coefficients.                                          # 
#   SEs     - vector of standard errors for the standardized         # 
#             regression coefficients.                               # 
#   alpha   - desired Type I error rate.                             #   
#   CIs     - (1-alpha)% confidence intervals for standardized       #
#             regression coefficients.                               #
######################################################################

seBetaCor <- function(R, rxy, Nobs, alpha = .05, digits = 3, 
                       covmat = 'normal') { 

#~~~~~~~~~~~~~~~~~~~~~~~~ Internal Functions ~~~~~~~~~~~~~~~~~~~~~~~~# 	
# Dp: Duplicator Matrix	 
  Dp <- function(p) {
    M <- matrix(nrow = p, ncol = p)
    M[ lower.tri(M, diag = T) ] <- seq( p*(p + 1)/2 )
    M[ upper.tri(M, diag = F) ] <- t(M)[ upper.tri(M, diag = F) ]
    D <- outer(c(M), unique(c(M)),
                FUN = function(x, y) as.numeric(x == y) )
    D
  }

# row.remove: Removes rows from a symmetric transition matrix 
# to create a correlation transition matrix  
# (see Browne & Shapiro, (1986); Nel, 1985).

  row.remove <- function(p) {
    p1 <- p2 <- p
    rows <- rep(1,p)
    for(i in 2:p) {
      rows[i] <- rows[i] + p1
      p1 <- p1 + (p2-1)
      p2 <- p2 - 1
    }
    rows
  }

# cor.covariance: Create a covariance matrix among predictor- 
# and predictor/criterion-correlations assuming multivariate
# normality (see Nel, 1985).  

  cor.covariance <- function(R, Nobs) {
  
  # Symmetric patterned matrix (Nel, p. 142)
    Ms <- function(p) {
    M <- matrix(c( rep( c( rep( c(1, rep(0, times = p*p +
               (p - 1) )), times = p - 1), 1,
                rep(0, times = p) ), times = p - 1),
                rep( c( 1, rep(0, times = p*p + (p - 1)) ),
                times = p - 1 ), 1 ), nrow = p^2)
    (M + diag(p^2))/2   
  }

  # Diagonal patterned matrix (Nel, p. 142).
    Md <- function(p) {
      pl <- seq(1,(p^2),by=(p+1))
      dg <- rep(0,p^2)
      dg[pl] <- 1
      diag(dg)
    }

  # Nel's (1985, p. 143) Psi matrix.
    Psi <- function(R) {
      p <- ncol(R)
      id <- diag(p)
      .5*(4*Ms(p) %*% (R %x% R) %*% Ms(p) - 2*(R %x% R)
          %*% Md(p) %*% (id %x% R + R %x% id) -
          2*(id %x% R + R %x% id) %*% Md(p) %*% (R %x% R) +
          (id %x% R + R %x% id) %*% Md(p) %*% (R %x% R) %*%
          Md(p) %*% (id %x% R + R %x% id))
    }
  
    p <- ncol(R)
  # Create symmetric transition matrix
    Kp <- solve(t(Dp(p)) %*% Dp(p)) %*% t(Dp(p))
  
  # Create correlation transition matrix
    Kpc <- Kp[-row.remove(p),] 
  
    (Kpc %*% Psi(R) %*% t(Kpc))/Nobs # The desired cov matrix
  }  # End cor.covariance
#~~~~~~~~~~~~~~~~~~~~~~ End Internal Functions ~~~~~~~~~~~~~~~~~~~~~~#
  
  R <- as.matrix(R)
  
  p <- ncol(R)
  Rinv <- solve(R)
  
  if(p == 1) {
  	beta.cov <- ((1 - rxy^2)^2)/(Nobs - 3)
  	ses <- sqrt(beta.cov)
  } else {

  # Covarianc matrix of predictor and predictor-criterion correlations
    sR <- rbind(cbind(R, rxy),c(rxy, 1))
  
    if(is(covmat)[1] == 'matrix') Sigma <- covmat
    else Sigma <- cor.covariance(sR, Nobs)

  # Create symmetric transition matrix (see Nel, 1985)
    Kp <- solve(t(Dp(p)) %*% Dp(p)) %*% t(Dp(p))

  # Create correlation transition matrix (see Browne & Shapiro, 1986).
    Kpc <- as.matrix(Kp[-row.remove(p),] )
    if(ncol(Kpc) == 1) Kpc <- t(Kpc)

  # Derivatives of beta wrt predictor correlations (Rxx)  
    db.drxx <- -2 * ( ( t( rxy ) %*% Rinv) %x% Rinv ) %*% t(Kpc)

  # Derivatives of beta wrt predictor-criterion correlations (rxy)
    db.drxy <- Rinv

  # Concatenate derivatives
    jacob <- cbind(db.drxx,db.drxy)

  # Reorder derivatives to match the order of covariances and 
  # variances in Sigma

    rxx.nms <- matrix(0,p,p)
    rxy.nms <- c(rep(0,p+1))
    for(i in 1:p) for(j in 1:p) rxx.nms[i,j] <- paste("rx",i,"rx",j,sep='')
    for(i in 1:p+1) rxy.nms[i] <- paste("rx",i,"y",sep='')

    nm.mat <- rbind(cbind(rxx.nms,rxy.nms[-(p+1)]),rxy.nms)
    old.ord <- nm.mat[lower.tri(nm.mat)]
    new.ord <- c(rxx.nms[lower.tri(rxx.nms)],rxy.nms)

    jacobian <- jacob[,match(old.ord,new.ord)]  

  # Create covariance matrix of standardized regression coefficients 
  # using the (Nobs-3) correction suggested by Yuan and Chan (2011)  

    beta.cov <- jacobian %*% Sigma %*% t(jacobian) * Nobs/(Nobs-3)
    beta.nms <- NULL   
    for(i in 1:p) beta.nms[i] <- paste("beta",i,sep='')
    rownames(beta.cov) <- colnames(beta.cov) <- beta.nms

    ses <- sqrt(diag(beta.cov))
  
  }
  
  CIs <- as.data.frame(matrix(0, p, 3))
  colnames(CIs) <- c("lbound", "estimate", "ubound")
  for(i in 1:p) rownames(CIs)[i] <- paste("beta_", i, sep='')

  tc <- qt(alpha / 2, Nobs - p - 1, lower.tail = FALSE)
  beta <- Rinv %*% rxy
 
  for(i in 1:p) {
    CIs[i,] <- c(beta[i] - tc * ses[i], beta[i], beta[i] + tc * ses[i])
  }

  cat("\n", 100 * (1 - alpha),
      "% CIs for Standardized Regression Coefficients: \n\n",sep='')

  print(round(CIs,digits))
  invisible(list(cov.Beta=beta.cov,se.Beta=ses,alpha=alpha,CI.beta=CIs)) 
}	
##############    END OF FUNCTION ###################################








                 