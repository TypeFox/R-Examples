mstep <- function( x, theta, phi )
{
  ## Purpose: M step for theta
  ## ----------------------------------------------------------------------
  ## Arguments: x - table of 0/1s
  ##            phi - prob of length
  ##            theta - intensity parm
  ## ----------------------------------------------------------------------
  ## Author: Charles Berry, Date:  7 May 2011, 11:34

  ## f <- function(theta) sum( -(1-x) * theta * phi + x * log( 1 - exp(-theta*phi)) )
  
  exptp <- exp( -phi %o% theta )
  denom <- 1 - exptp
  ## deprecated  denom[ x == 0 ] <- 1.0
  grad <- colSums( - (1-x) * phi  + x * phi * exptp / denom )
  curv <- - colSums( phi^2 * exptp/ denom )
  theta <- theta - grad / curv
  list(theta=theta,grad=grad,info=curv )
}
