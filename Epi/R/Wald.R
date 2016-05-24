Wald <-
function( obj, H0=0, ... )
{
rl <- ci.lin( obj, ..., vcov=TRUE )
beta <- rl$est
vcov <- rl$vcov
if( missing( H0 ) ) H0 <- beta*0
if( length(H0) != length(beta) ) stop( "H0 has length ", length(H0),
          " but the set of selected paramteters has length ",
          length(beta), ":\n",
          paste(round(beta,options()[["digits"]]),collapse=" ") )
chi <- t( beta-H0 ) %*% solve( vcov, beta-H0 )
 df <- length(beta)
  p <- 1 - pchisq( chi, df )
c( "Chisq"=chi, "d.f."=df, "P"=p )
}

