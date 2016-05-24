#################################################
# Generating correlation matrices               #
# with prespecified Eigen Values                #
# Based on Marsaglia and Olkin (1984) Algorithm #
#                                               #
# Jeff Jones                                    #
#                                               #
# Input:                                        #
#                                               #
#   eigenval - vector of desired eigen values   #
#   seed     - set.seed(seed)                   #
#                                               #
# Output:                                       #
#                                               #
#	 R        - correlation matrix               #
#################################################
#################################################


genCorr <- function(eigenval,seed="rand"){
	
  if(!isTRUE(all.equal(sum(eigenval),length(eigenval)))) stop("Sum of eigenvalues not equal to Number of Variables\n")
  if(seed!="rand") set.seed(seed)

  norm <- function(x) x/as.numeric( sqrt(x %*%t(x)))

## step (i) in Marsaglia & Olkin

  Nvar <- length(eigenval)
  L <- diag(eigenval)
  E <- diag(Nvar)
  P <- matrix(0,Nvar,Nvar)

  if(identical(eigenval,rep(1,Nvar))) return(diag(Nvar))
  else {
  for(i in 1:(Nvar-1)) {

      xi <- as.vector(rnorm(Nvar))                  # step (ii) 
      xi <- xi%*%E                                  # create vector in rowspace of E

      a <- (xi %*% (diag(Nvar)-L) %*% t(xi))

      d_sq <- -1

      while(d_sq <= 0) {                             # if b^2 - ac < 0, recompute eta

        eta <- as.vector(rnorm(Nvar))                # step (iii)
        eta <- eta%*%E                               # create another vector in rowspace of E

        b <- (xi %*% (diag(Nvar)-L) %*% t(eta))      # step (iv)

        cc <- (eta %*% (diag(Nvar)-L) %*% t(eta))

        d_sq <- ((b^2) - (a*cc))

      } # end while

    ## step (v)

    r <- as.numeric(((b + sign(runif(1,-1,1))*sqrt(d_sq)))/(a))
    zeta.tmp <- (r*xi - eta)
    zeta <- sign(runif(1,-1,1))*norm(zeta.tmp)

    P[i,] <- zeta
    E <- E - t(zeta)%*%zeta
  
  } # end for

  P[Nvar,] <- norm(as.vector(rnorm(Nvar))%*%E)      # normalized random vector in rowspace of E

  R <- P%*%L%*%t(P)                                 # correlation matrix

  R     
  }

}

