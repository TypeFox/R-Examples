
rvaluesMCMC <- function(output, qtheta, alpha.grid = NULL, ngrid = NULL,
                        smooth = "none")  {
  
  ##########################################################################
  ##  Last edited: 7/10/14
  ##  Input
  ##    output - an array of mcmc output (num_units x num_posterior_samples)
  ##    qtheta - either a function computing theta quantiles (upper tail) or a 
  ##             vector of theta quantiles
  ##    alpha.grid - grid of values over (0,1); used for the discrete approx.
  ##                 approach for computing r-values.
  ######
  ##   Output
  ##     object of class rvalues 
  ###########################################################################
  
    nunits <- nrow(output)
    rvalues <- rep(1,nunits)
    if(is.null(alpha.grid)) {
       ### initialize alpha_grid if not entered by user
       alpha.grid <- MakeGrid(nunits,type="log",ngrid=ngrid)
    }
    if(!is.null(alpha.grid)) {
        ngrid <- length(alpha.grid)
    
        # check that alpha.grid values are strictly between 0 and 1  
        alpha.check <- all((alpha.grid > 0) & (alpha.grid < 1))
        if(!alpha.check) {
            stop("All the alpha.grid values must be strictly between zero and one")
        }
    }
    
    ### check qtheta
    if(class(qtheta)!="function" & class(qtheta) != "numeric") {
      stop("qtheta must either be a function or a vector")
    }
    if(class(qtheta)=="function") {
      ### if qtheta is a function
      theta.alpha <- qtheta(alpha.grid)
    }
    else {
      ### if qtheta is a vector
      if(length(theta.alpha) != length(alpha.grid)) {
        stop("qtheta and alpha.grid must have the same length") 
      }
      theta.alpha <- qtheta
    }
    ### theta_alpha represents the quantiles of theta corresponding to 
    ### the alpha.grid

    ngrid <- length(theta.alpha)
       
    ## Valpha[i,j] = P(\theta_{i} \geq \theta_{\alpha_{j}}|X_{i}), where,
    ## marginally: P(\theta_{i} \geq \theta_{\alpha_{j}}) = 1 - \alpha_{j}
 
    V <- matrix(nrow = nunits, ncol = ngrid)
    for(i in 1:nunits) {
        V[i,] <- ecdff(output[i,], theta.alpha)
    }
    
    ## Compute "raw" \lambda_{\alpha} function
    lam <- numeric(ngrid)
    for( j in 1:ngrid )  {
       lam[j] <- quantile(V[,j], probs = 1 - alpha.grid[j], names = FALSE, type = 1)
    }
    ## Smooth \lambda_{\alpha}
    if(smooth=="none") {
        cc2 <- approxfun(alpha.grid, lam, yleft = 1, yright = 0)
        lam.smooth.eval <- cc2(alpha.grid)
        lam.smooth <- approxfun( c(0,alpha.grid,1), c(1,lam.smooth.eval,0))
    }
    else {
        cc2 <- supsmu( alpha.grid, lam, bass= smooth )
        lam.smooth.eval <- cc2$y
        lam.smooth <- approxfun( c(0,cc2$x,1), c(1,cc2$y,0))
    }
    ### For each row of Valpha, determine the index at which
    ### Valpha[i,] intersects lambda_{\alpha}
    rvalues <- VVcut(V, lam.smooth.eval, nrow(V), ngrid, alpha.grid)
    
    #rvalues <- alpha.grid[cut.ind]
 
    ### add stuff when returning ...
    thetaPM <- rowMeans(output)
    ans <- list()
    class(ans) <- "rvals"
    ord <- order( rvalues, -thetaPM )
    bar <- data.frame( RValue=rvalues, RV.rank=rank(rvalues),
                     PostMean=thetaPM,PM.rank=rank(-thetaPM))
    ans$main <- bar[ord,]
    ans$aux <- list(alpha.grid=alpha.grid,unsorted=bar,V=V,
                    Vmarginals=lam, Vmarginals.smooth = lam.smooth)
    ans$rvalues <- rvalues
    return(ans)
}

ecdff <- function(x,theta.alpha) 
{
  x <- sort(x,decreasing=T)
  n <- length(x)
  vals <- unique(x)
  if(length(vals)==1) {
     ## In this case the empirical cdf is a point mass at vals
     return(ifelse(theta.alpha >= vals,0,1))
  }
  else {
     rval <- approxfun(vals, cumsum(tabulate(match(x, vals)))/n, yleft = 1, yright = 0)
     return(rval(theta.alpha))
  }
}
