info.ordinal.cat <- function( model, link, theta, x, icat ){

### Returns the contribution to the information matrix from an observation
### of icat.

### model is 1 for linear model
###          2 for quadratic model

### link is 1 for logistic link
###         2 for complementary log-log

### theta is vector of parameter values
###         theta[1:(ncat-1)] are intercept terms
###         theta[ncat] is coefficient of x
###         theta[ncat+1] is coefficient of x^2 for quadratic

### x scalar value of covariate

### icat is outcome value for which contribution to into evaluated

  ntheta <- length(theta)

  if ( model == 1 ) ncat <- ntheta   else 
          ncat <- ntheta - 1

  Du.i <- Du.im1 <- rep(0, ntheta)
  Dg.i <- Dg.im1 <- rep(0, 3)

  not.first <- icat > 1
  not.last <- icat < ncat
  if (ncat == 1) not.last <- 1

### Set up value, first and second derivative of link with respect to u
### for u_i and u_im1

  if (not.last) {
    u.i <- theta[icat] + theta[ncat]*x
    if (model == 2) u.i <- u.i + theta[ncat+1]*x^2 
    Dg.i <- derivs.link( u.i, link )  
    Du.i[icat] <- 1
    Du.i[ncat] <- x
    if (model ==2) Du.i[ntheta] <- x^2  }

  if (not.first) {
    u.im1 <- theta[icat-1] + theta[ncat]*x
    if (model == 2) u.im1 <- u.im1 + theta[ncat+1]*x^2 
    Dg.im1 <- derivs.link( u.im1, link )  
    Du.im1[icat-1] <- 1
    Du.im1[ncat] <- x
    if (model ==2) Du.im1[ntheta] <- x^2  }

### Compute p.i, probability of landing in icat

  if (ncat == 1 || icat == 1) p.i <- Dg.i[1] else {
      if (icat < ncat) p.i <- Dg.i[1] - Dg.im1[1]  else 
             p.i <- 1 - Dg.im1[1]   }

### Compute second derivatives multiplied by p.i

  ans <- outer(Du.i,Du.i) * ( Dg.i[3] - Dg.i[2]^2 / p.i ) -
    outer(Du.im1,Du.im1) * ( Dg.im1[3] + Dg.im1[2]^2 / p.i ) + 
      ( outer(Du.i,Du.im1) + outer(Du.im1,Du.i) ) * 
         ( Dg.i[2] * Dg.im1[2] / p.i )

### Return negative of second derivatives multiplied by p.i

  return( -ans )    }
