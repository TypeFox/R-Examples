mice.impute.2l.latentgroupmean.MCMC <- function (y, ry, x, type , 
                    pls.facs = NULL , imputationWeights = NULL ,
                    interactions = NULL , quadratics = NULL , 
                    mcmc.burnin=100 , mcmc.iter=1000 , 
					EAP = FALSE , ...){  

    # latent group mean
    cluster <- as.numeric( x[ , type == - 2] )
    covariates <- as.matrix( x[ , type == 1 ] )
    colnames(covariates) <- colnames(x)[ type == 1 ]
    y <- x[ , type == 2 ]

    # distinguish cases with and without covariates
    if ( sum( type == 1 ) > 0){ 

        # aggregate covariates
        cov2 <- cbind( cluster , covariates )
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

        rownames(covaggr) <- NULL # prevent warning
        mcmcdf <- data.frame(y,cluster,covaggr)
        covnames <- paste0("x",1:ncol(covaggr))
        colnames(mcmcdf) <- c("y","cluster",covnames)
        mcmcfml <- stats::as.formula(paste0(c("y~1",covnames),collapse="+"))

        #yaggr <- mice.impute.2l.groupmean( y , ry , x=cbind(cluster,y) , type=c(-2,1) ,
        #         grmeanwarning=FALSE )
        #prior <- list(R=list(V=max(var(y-yaggr),.01), nu=1),
        #              G=list(G1=list(V=max(var(yaggr),.01), nu=1)))
        prior <- list(R=list(V=1, nu=1), G=list(G1=list(V=1, nu=1)))

        mod <- MCMCglmm::MCMCglmm(mcmcfml, random=~cluster, data=mcmcdf,
                        thin=1, burnin=mcmc.burnin, nitt=mcmc.burnin+mcmc.iter, 
                        verbose=FALSE)

    } else {

        mcmcdf <- data.frame(y=y,cluster=cluster)
        mcmcfml <- y~1

        prior <- list(R=list(V=1, nu=1), G=list(G1=list(V=1, nu=1)))

        mod <- MCMCglmm::MCMCglmm(mcmcfml, random=~cluster, data=mcmcdf,
                        thin=1, burnin=mcmc.burnin, nitt=mcmc.burnin+mcmc.iter, 
                        verbose=FALSE)

    }

    modf <- stats::model.matrix(mcmcfml, data=mcmcdf) %*% apply(mod$Sol,2,mean)
    psi2 <- mean(mod$VCV[,1])
    sig2 <- mean(mod$VCV[,2])

    # aggregate: cluster size, fixed prediction, residual from fixed
    a1 <- stats::aggregate( cbind( 1 , modf , y-modf ) , list( cluster) , sum )
    a1[,3] <- a1[,3]/a1[,2] # average from fixed
    a1[,4] <- a1[,4]/a1[,2] # average deviation from fixed
    a1[,5] <- ( psi2 / (psi2+sig2/a1[,2]) ) * a1[,4] # EAPs of random effects
    a1[,6] <- psi2*sig2/(sig2 + psi2*a1[,2])

    # match cluster indices
    ind <- match( cluster , a1[,1] )
    ximp <- stats::rnorm( nrow(a1) , mean = a1[,3]+a1[,5] , 
				sd = (1 - EAP )*sqrt(a1[,6]) )[ ind ]
    return(ximp)

}
