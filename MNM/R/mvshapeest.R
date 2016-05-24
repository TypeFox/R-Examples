mv.shape.est <- function(X,score="identity", estimate="outer", location=NULL, na.action=na.fail, ...)
    {
    score <- match.arg(score,c("identity", "sign", "rank", "symmsign"))
    estimate <- match.arg(estimate,c("inner", "outer"))
    
    X<-na.action(X)
    
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    X<-as.matrix(X)
    n<-dim(X)[1]
    p<-dim(X)[2]
    
    X.names <- colnames(X)
    
    if (!is.null(location)) if (any(c(!is.vector(location), length(location)!=p))) stop("'location must be NULL or a vector of length p'")
    
    covMat <- switch( score,
        
        "identity" = if (is.null(location)) cov(X) else covOrigin(X,location=location)
        ,
        "sign" = switch(estimate,
                        "outer" =  if (is.null(location)) SCov(X, location=spatial.median(X), ...) else SCov(X, location=location)
                        ,
                        "inner" = {
                                    if (is.null(location)) SHAPE <- HR.Mest(X,...)$scatter else SHAPE <- tyler.shape(X, location=location, ...)
                                    p* SHAPE / sum(diag(SHAPE))
                                  }  
                        ) 
        ,
        "symmsign" = switch(estimate,
                        "outer" =  {
                                    if (!is.null(location)) warning("when 'score = symmsign' then location will be ignored")
                                    SSCov(X)
                                    }
                        ,
                        "inner" = {
                                    if (!is.null(location)) warning("when 'score = symmsign' then location will be ignored")
                                    SHAPE <- duembgen.shape(X,...)
                                    p* SHAPE / sum(diag(SHAPE))
                                    }
                        ) 
        ,
        "rank" = switch(estimate,
                        "outer" =  {
                                    if (!is.null(location)) warning("when 'score = rank' then location will be ignored")
                                    RCov(X)
                                    }
                        ,
                        "inner" = {
                                    if (!is.null(location)) warning("when 'score = rank' then location will be ignored")
                                    SHAPE <- rank.shape(X,...)
                                    p* SHAPE / sum(diag(SHAPE))
                                    }
                        )
        )
    
    colnames(covMat) <- X.names
    rownames(covMat) <- X.names
    return(covMat)
    }
