mice.impute.weighted.pmm <- function (y, ry, x,  imputationWeights = NULL , 
                                    pls.facs = NULL ,  interactions = NULL , quadratics = NULL ,  ...){
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
    plsout <- .aux.pls.imputation( newstate = newstate , vname = vname , pls.impMethod = "pmm" , 
                    x = x[,-1] , y = y , ry=ry , imputationWeights = imputationWeights , 
                    interactions = interactions , quadratics = quadratics ,  pls.facs = pls.facs ,  ... )
    # save PLS result
    pls.facs <- plsout$pls.facs
    yimp <- plsout$yimp        
    if (is.null(pls.facs)){
            parm <- .weighted.norm.draw( yobs = yobs , xobs = xobs , ry = ry , y = y , x = x ,
                                weights.obs = weights.obs , ... )   
            yhatobs <- x[ry,] %*% parm$coef
            yhatmis <- x[!ry,] %*% parm$beta
            yimp <- apply(as.array(yhatmis), 1, .pmm.match, yhat = yhatobs,
                                y = y[ry], ... )
                            }    
    return(yimp)
    
    }
