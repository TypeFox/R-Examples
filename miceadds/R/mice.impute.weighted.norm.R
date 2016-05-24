mice.impute.weighted.norm <- function(y, ry, x, ridge = .00001 , pls.facs = NULL ,
                    imputationWeights = NULL ,
                    interactions = NULL , quadratics = NULL , 
                     ...){
    x <- cbind(1, as.matrix(x))
    # .weighted.norm.draw
    xobs <- x[ry,]
    yobs <- y[ry]
    if ( is.null( imputationWeights ) ){ imputationWeights <- rep(1 , length(y) ) }
    weights.obs <- imputationWeights[ ry   ]
    # standardize all weights to one
    weights.obs <- length(weights.obs) * weights.obs / sum( weights.obs )
    
    #.+.+.+.+.+.+.+.+.+.+.+.+
    # PLS interactions and quadratics
    newstate <- get( "newstate" , pos = parent.frame() )  
    vname <- get("vname", pos = parent.frame()) # get variable name         
    plsout <- .aux.pls.imputation( newstate = newstate , vname = vname , pls.impMethod = "norm" , 
                    x = x[,-1] , y = y , ry=ry , imputationWeights = imputationWeights , 
                    interactions = interactions, quadratics = quadratics ,  pls.facs = pls.facs ,  ... )
    pls.facs <- plsout$pls.facs
    yimp <- plsout$yimp
   
      if ( is.null(pls.facs) ){ 
          parm <- .weighted.norm.draw( yobs = yobs , xobs = xobs , ry = ry , y = y , x = x ,
                              weights.obs = weights.obs , ... )   
          yimp <- x[!ry,  ] %*% parm$beta + stats::rnorm(sum(!ry)) * parm$sigma
                                }
      return(yimp)
    }

	
	
	
#***** .weighted.norm.draw *******************
.weighted.norm.draw <- function( yobs , xobs , ry , y , x , weights.obs , ridge = .00001 , ... ){
    WW <- diag( weights.obs )
    # X'*W*X
    xtx <- t(xobs) %*% WW %*% xobs
    pen <- ridge * diag(xtx)
    if (length(pen)==1) pen <- matrix(pen)
    v <- solve(xtx + diag(pen))
    # V * ( X'*W*Y)
    coef <- t(yobs %*% WW %*% xobs %*% v)
    residuals <- yobs - xobs %*% coef
    # calculate weighted residuals
    residuals2 <- weights.obs * residuals
    sigma.star <- sqrt(sum((residuals2)^2)/ stats::rchisq(1, sum(ry) - ncol(x)))  
    beta.star <- coef + (t(chol((v + t(v))/2)) %*% stats::rnorm(ncol(x))) * sigma.star
    parm <- list(coef, beta.star, sigma.star)     
    names(parm) <- c("coef","beta", "sigma") 
    return(parm)
    }
#----------------------------------------------------------------------
