PostSummaries <- function(post.means, post.sds, hypers = NULL, qtheta = NULL,
                           alpha.grid = NULL, ngrid = NULL, smooth = 0) {
   #############################################################
   ###  Last edited: 7/13/14
   ###  Input
   ###     post.means -  a vector of posterior means
   ###     post.sds - a vector posterior standard deviations
   ###        Note: posteriors are assumed to be Normal
   ###     hypers - a list of hyperparameters for a Normal prior
   ###     qtheta - if hypers is not entered, a function for computing
   ###               the upper tails of theta
   ###     alpha.grid -  a grid of alpha-values in (0,1)
   ###     ngrid  - length of alpha.grid
   ###     smooth - smoothing parameter for the Vmarginal curve   
   #############################################################
  
   if(is.null(qtheta) & is.null(hypers)) {
       stop("You must enter a list of hyperparameter or a quantile 
             function for theta")
   }
   nunits <- length(post.means)
   rvalues <- rep(1,nunits)
   if(is.null(alpha.grid)) {
       ### initialize alpha.grid if not entered by user
       alpha.grid <- MakeGrid(nunits, type = "log", ngrid = ngrid)
   }
   if(!is.null(alpha.grid)) {
    
       # check that alpha.grid values are strictly between 0 and 1  
       alpha.check <- all((alpha.grid > 0) & (alpha.grid < 1))
       if(!alpha.check) {
           stop("All the alpha.grid values must be strictly between zero and one")
       }
   }
   
   if(is.null(qtheta)) {
        if(class(hypers)!="list") {
             stop("hypers must be a list with two components")
        }
        theta.alpha <- qnorm(alpha.grid,mean=hypers$mean,sd=hypers$sd,lower.tail=F)
        Thetfun <- function(alpha) { qnorm(alpha,mean=hypers$mean,sd=hypers$sd,lower.tail=F) }  
   }
   else {
       if(!is.null(hypers)) {
           if(class(hypers)!="list") {
               stop("hypers must be a list with two components")
           }
       
           ### This is the case where both qtheta and hypers are non-NULL
           warning("The prior will be determined by the hypers argument")
            
           theta.alpha <- qnorm(alpha.grid,mean=hypers$mean,sd=hypers$sd,lower.tail=F)
           Thetfun <- function(alpha) { qnorm(alpha,mean=hypers$mean,sd=hypers$sd,lower.tail=F) } 
       }
       else {
           ### qtheta should be a function since the r-values are determined
           ### through root-finding
       
           if(class(qtheta)!="function") {
               stop("qtheta must be a function")
           }
           theta.alpha <- qtheta(alpha.grid)
           Thetfun <- qtheta
       }
   }

   ngrid <- length(alpha.grid)
   lam <- numeric(ngrid)
   for( j in 1:ngrid)
   {
       V <- pnorm(theta.alpha[j],mean=post.means,sd=post.sds,lower.tail=F)
       lam[j] <- quantile(V, probs = 1 - alpha.grid[j], names = FALSE, type = 1)
   }
   
   if(smooth=="none") {
        cc2 <- approxfun(alpha.grid, lam, yleft = 1, yright = 0)
        lam.smooth <- approxfun( c(0,cc2$x,1), c(1,cc2$y,0))
    }
    else {
        cc2 <- supsmu( alpha.grid, lam, bass = smooth)
        lam.smooth <- approxfun( c(0,cc2$x,1), c(1,cc2$y,0))
    }

    dfun <- function(alpha, pm, psd)
    {
       dd <- lam.smooth(alpha) - pnorm(Thetfun(alpha),mean=pm,sd=psd,lower.tail=F)
       dd
    }
   
    rvalues <- rep(0,nunits)
    for(i in 1:nunits) {
       rvalues[i] <- uniroot(dfun, lower = 0, upper = 1,pm = post.means[i], psd = post.sds[i])$root
    }
   
    #rvalues <- Mroot(dfun_multi,lower=rep(0,nunits),upper=rep(1,nunits))

    ans <- list()
    class(ans) <- "rvals"
    ord <- order( rvalues )
    bar <- data.frame( RValue=rvalues, RV.rank=rank(rvalues),
                     PostMean=post.means,PM.rank=rank(-post.means))
    ans$main <- bar[ord,]
    ans$aux <- list(V=V,alpha.grid=alpha.grid,unsorted=bar,Vmarginals=lam,
                   Vmarginals.smooth=lam.smooth)
    ans$rvalues <- rvalues
    return(ans)
}

