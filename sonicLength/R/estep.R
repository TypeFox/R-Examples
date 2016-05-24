 Ey.given.x <- function(x,theta,phi )
  {
    ## Purpose: E step
    ## ----------------------------------------------------------------------
    ## Arguments: x - table of 0/1s
    ##            phi - prob of length
    ##            theta - intensity parm
    ## ----------------------------------------------------------------------
    ## Author: Charles Berry, Date: 26 Mar 2011, 16:45
    lambda <- matrix(phi,ncol=1) %*% matrix(theta,nrow=1)
    expnl <- exp(-lambda)
    denom <- ( 1 - expnl ) #dpois(0, lambda ))
    res <- x * lambda / denom
    if (length(wd <- which(denom==0))) res[ wd ] <- x [ wd ] # just in case lambda is really small
    attr(res,"logLik") <- sum( -lambda[ x==0 ],  log1p( -expnl[x==1] ))
    res
  }


pr.y.given.x <- function(x, theta,phi,kmax=20)
  {
    ## Purpose: distribution of sonicants
    ## ----------------------------------------------------------------------
    ## Arguments: x -table of 0/1s a matrix whose rows may
    ##                    correspond to unique lengths (phi) with rownames
    ##            theta -  like maxEM$theta
    ##            phi -   like maxEM$phi
    ##            kmax - highest count to bother with (all higher values
    ##                        are globbed together in the result)
    ## ----------------------------------------------------------------------
    ## Author: Charles Berry, Date: 19 Jun 2011, 12:44
    lambda <- phi %o% theta
    lambda.list <- apply( lambda*x, 2 , function(x) x[ x > 0 ] )
    sapply( lambda.list,
           function(x) {
             lx <- length(x)
             if (lx == 1 ) {
               res <- dpois( 1:kmax, x ) / (1 - exp( -x ) )
               c(res, pmax( 0.0, 1-sum(res)))
             } else {
               if ( lx > kmax ) {
                 rep(0:1,c(kmax,1))
               } else {
                 ## do a convolution
                 parts <- outer(x,1:kmax, function(x,y) dpois(y,x)/(1-exp(-x)))
                 parts <- cbind( parts, 1 - rowSums( parts ))
                 cres <- parts[1,]
                 for ( i in 2:nrow(parts)) cres <- convolve(cres,rev(parts[i,]),type='open')
                 res <- c( rep( 0 , lx-1 ), cres[ 1:(kmax-lx+1)])
                 res[kmax+1] <- pmax(0.0, 1 - sum(res) )
                 res
               }
             }
           })
  }

