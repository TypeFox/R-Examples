mice.impute.2l.latentgroupmean.ML <- function (y, ry, x, type , 
                    pls.facs = NULL , imputationWeights = NULL ,
                    interactions = NULL , quadratics = NULL , 
					EAP = FALSE , ...){  
    # latent group mean
    cluster <- as.numeric( x[ , type == - 2] )
    covariates <- as.matrix( x[ , type == 1 ] )
    colnames(covariates) <- colnames(x)[ type == 1 ]
    y <- x[ , type == 2 ]
    # distinguish cases with and without covariates
    if ( sum( type == 1 ) > 0){ 
        cov2 <- cbind( cluster , covariates )
        # aggregate covariates
        covaggr <- mice.impute.2l.groupmean( y , ry , x = cov2 , type = c(-2 , rep(1,ncol(covariates) ) ) ,
                       grmeanwarning=FALSE )
        colnames(covaggr) <- colnames(covariates)
        # aggregation at level 2
        covaggr.l2 <- stats::aggregate( covaggr , list( cluster ) , mean )
        colnames(covaggr.l2)[-1] <- colnames(covaggr)
        y.l2 <- stats::aggregate( y , list( cluster ) , mean , na.rm=T )
        h1 <- as.matrix(covaggr.l2[,-1] )
        colnames(h1) <- colnames(covaggr)
        # PLS interactions and quadratics
        newstate <- get( "newstate" , pos = parent.frame() )  
        vname <- get("vname", pos = parent.frame()) # get variable name         
        plsout <- .aux.pls.imputation( newstate = newstate , vname = vname , pls.impMethod = "xplsfacs" , 
                        x = h1 , y = y.l2[,2] , ry= ( ! is.na(y.l2[,2] )) , 
                        imputationWeights = rep( 1 , nrow(covaggr.l2)) , 
                        interactions = interactions, quadratics = quadratics ,  pls.facs = pls.facs ,  ... )
        # imputation PLS
        if( ! is.null( plsout$yimp  ) ){ 
            covaggr.l2r <- as.matrix(plsout$yimp[,-1])
            covaggr <- as.matrix( covaggr.l2r[ match( cluster , covaggr.l2[,1] ) , ] )
                }		
        mod <- lme4::lmer( y ~ as.matrix( covaggr  ) + ( 1 | cluster ) )
                    } else {
        mod <- lme4::lmer( y ~ ( 1 | cluster ) )
                    }							
#    modr <- ranef( mod , postVar = TRUE )
    modr <- lme4::ranef( mod , condVar = TRUE )
    modr <- modr$cluster	
#    modf <- fixef( mod )	
	modf <- mod@resp$mu		
    # extract cluster indices
    a1 <- stats::aggregate( modf , list( cluster) , mean )
    a1[,3] <-  sqrt( attr( modr , "postVar" )[1,1,] )
    # match cluster indices
    ind <- match( cluster , a1[,1] )
    ximp <- stats::rnorm( nrow(a1) , mean = a1[,2] , 
				           sd = (1-EAP)*a1[,3] )[ ind ]
    return(ximp)
}
