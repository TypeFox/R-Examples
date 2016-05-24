"poi.beta.laplace" <-
function( data, alpha=get.alpha(data), phi=1, nu=0.1,
                           beta=seq(from=0,to=7,length=1000) )
{
  h <- function(x, P)
    P$phi*x - (P$n*P$alpha+P$nu-1)*log(x) + sum( (P$alpha+P$y)*log(P$X+x) )
  h2 <- function( u, phi, n, nu, y, X, alpha )
    (n*alpha+nu-1)/u^2 - sum( (alpha+y)/(u+X)^2 )
  y <- data$y
  X <- data$x
  n <- length(data$y)
  P <- list( phi=phi, nu=nu, n=n, y=y, X=X, alpha=alpha )
  h.min <- optim( par=1, fn=h, method="L-BFGS-B", hessian=TRUE, lower=0, upper=Inf, P=P )
  u <- h.min$parameters
  conv <- h.min$message
  H2 <- det(h.min$hessian)
  out <- sum( lgamma(alpha+y)-lgamma(y+1)-lgamma(alpha)+y*log(X)) - 
       lgamma(nu) + nu*log(phi) -
       h.min$value+0.5*log( 2*pi/H2 )
  PDF <- beta
  for (i in 1:length(beta))
  PDF[i] <- - h(beta[i],P) - (- h.min$value + 0.5*log( 2*pi/H2 ) )
  int <- sum(exp(PDF))*(beta[2]-beta[1])

  list( int=int, conv=conv, x=beta, y=exp(PDF) )
}

