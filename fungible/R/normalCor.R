#---------------------------------------------------------------#
# Covariance matrix of correlations                             #
#                                                               #
# Nel, D.G. (1985). A matrix derivation of the                  #
# asymptotic covariance matrix of sample correlation            #
# coefficients. Linear algebra and its applications,            #
# 67:137-145.                                                   #
#                                                               #
# Arguments:                                                    #
#   R      - A p x p correlationmatrix                          #
#   Nobs   - number of observations                             #
#                                                               #      
# Output                                                        #
#          - a covariance matrix of correlations                #
#                                                               #
#---------------------------------------------------------------#

normalCor <- function(R, Nobs) {
	
# Duplicator Matrix (Nel, p. 138).	
	
  Dp <- function(p) {
    M <- matrix(nrow = p, ncol = p)
  	M[ lower.tri(M, diag = T) ] <- seq( p*(p + 1)/2 )
  	M[ upper.tri(M, diag = F) ] <- t(M)[ upper.tri(M, diag = F) ]
  	D <- outer(c(M), unique(c(M)), 
               FUN = function(x, y) as.numeric(x == y) )
  	D
  }

# Symmetric patterned matrix (Nel, p. 142)

  Ms <- function(p) {
   M <- matrix(c( rep( c( rep( c(1, rep(0, times = p*p + 
                  (p - 1) )),
                  times = p - 1), 1, rep(0, times = p) ),
                  times = p - 1),
                  rep( c( 1, rep(0, times = p*p + (p - 1)) ), 
                  times = p - 1 ), 1 ),
                  nrow = p^2)
   (M + diag(p^2))/2
   }


# Diagonal patterned matrix (Nel, p. 142).

  Md <- function(p) {
    pl <- seq(1,(p^2),by=(p+1))
    dg <- rep(0,p^2)
    dg[pl] <- 1
    diag(dg)
  }

# Nel's (p. 143) Psi matrix. 

  Psi <- function(R) {
    p <- ncol(R)
    id <- diag(p)
    .5*(4*Ms(p) %*% (R %x% R) %*% Ms(p) - 2*(R %x% R) 
        %*% Md(p) %*% (id %x% R + R %x% id) -
        2*(id %x% R + R %x% id) %*% Md(p) %*% (R %x% R) +
        (id %x% R + R %x% id) %*% Md(p) %*% (R %x% R) %*% 
        Md(p) %*% (id %x% R + R %x% id))
   }
   
# Remove rows from Kp to convert a symmetric 
# transition matrix into a correlation transition 
# matrix (Nel, p. 143).

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
  
  p <- ncol(R)

# Create symmetric transition matrix  
  Kp <- solve(t(Dp(p)) %*% Dp(p)) %*% t(Dp(p))

# Create correlation transition matrix (see Shapiro & Browne, 1986).
  
  Kpc <- Kp[-row.remove(p),]
  
  normalCovMat<-(Kpc %*% Psi(R) %*% t(Kpc))/Nobs # The desired cov matrix

  normalCovMat

 }  # End cor.covariance










