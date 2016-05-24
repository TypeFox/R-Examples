########################################## 
# AUXILIAR FUNCTIONS AND PACKAGES REQUIRED
########################################## 

# KernSmooth package.
#source("distpred_aux2.R")# AUXILIARES1 file.
#source("distpred_aux.R")# AUXILIARES file.
#require(KernSmooth)
################################ 
# INPUT PARAMETERS
################################ 

# load series dataset as a list (thus permitting to deal with series of different length)


diss.PRED = function( x, y, h, B=500, logarithm.x=FALSE, logarithm.y=FALSE, differences.x=0, differences.y=0, plot=FALSE) {
    .ts.sanity.check(x, y)
    series <- list(x=x,y=y)
    logarithms = c(logarithm.x, logarithm.y)
    differences = c(differences.x, differences.y)
    # GENERAL INPUT PARAMETER VALUES
    #B <- 1000  		# number of bootstrap resamples 
    #k <- 5			# length of the forecast horizon
    # d <- 1			# autoregressive order
    
    
    # PARAMETERS REQUIRED TO GENERATE PREDICTIONS
    number.x <- 401	  	# number of equispaced points where estimators will be computed 	
    length.grid <- 400	# size of bands-grid to compute the cross validation selector
    nw.band <- 2		# nw.band determines the kind of automatic bandwidth selector 
    # considered to construct the nadaraya-watson estimator:
    # 1 when the cross-validation bandwidth selector is used
    # 2 when the plug-in bandwidth by Ruppert,Sheater and Wand is used
    l <- 4			# neighborhood of X_i whose observations are left out to estimate 
    # m(X_i) in the cross-validation algorithm. So: X_i is left out 
    # when l=0, X_{i-1},X_i, X_{i+1} when l=1, X_{i-2},X_{i-1},X_i, 
    # X_{i+1}, X_{i+2} when l=2 and so on.
    innov.sobrantes <- 100	# number of innovations generated to be burnt 
    
    # PARAMETERS REQUIRED TO CONSTRUCT THE BOOTSTRAP FORECAST DENSITIES 
    type.bw.forecast.dens <- "SJ-dpi" 
    # bandwidth selector to estimate the bootstrap forecast densities. 
    # Possible choice belongs to the set ("nrd0", "nrd", "ucv", "bcv", 
    # "SJ-ste", "SJ-dpi")I (see "bw.nrd" and "density" functions in R.
    
    
    k <- h
    
    ####################################################################################### 
    # MAIN PROGRAM
    ####################################################################################### 
    
    N <- length(series)						# number of series
    nombres <- names(series) 					# names of series
    
    ########################################################################################
    # 0. TRANSFORMATION ( logarithm and differences)
    ########################################################################################
    
    orig.series = cbind(x,y) #store untransformed series
    logarithms = 1 - logarithms #the transformation functions take 0 as do logarithm, counter-intuitive
    trans.series = transf.TRAMO( cbind(x,y), logarithms, differences )
    series = list( x= trans.series$T.Ser[,1], y= trans.series$T.Ser[,2] )
    
    
    
    
    
    #######################################################################################
    # 1. COMPUTING BOOTSTRAP VECTORS OF PREDICTIONS OF LENGTH k
    #######################################################################################
    
    # The bootstrap prediction paths will be stored in 
    path.k.prediction <- array(0,dim=c(k,B,N))
    dimnames(path.k.prediction)[[3]] <- nombres
    
    for ( punt in 1:N )
    {
        # Threshold to compute the nonparametric estimator 
        CM <- 3*sd(series[[punt]])
        # Nonparametric approach to generate bootstrap predictions: 
        aux <- pred.AB.CB(series[[punt]], nw.band, length.grid, l, CM, B, innov.sobrantes, 2, k, 0, 1) #only conditional BOOTSTRAP				
        # conditional bootstrap 
        path.k.prediction[,,punt] <- aux$CB						
        
    }
    
    #######################################################################################
    # 1.1  BACKTRANSFORM THE PREDICTIONS
    #######################################################################################
    
    back.pred = backtransf.TRAMO( path.k.prediction, trans.series$M.L.Ser, logarithms, differences, trans.series$Medias.log.series )
    path.k.prediction = back.pred
    
    
    ##########################################################################################
    # 2. COMPUTING BOOTSTRAP FORECAST DENSITIES
    ##########################################################################################
    
    # Bootstrap predictions at horizon k 
    k.prediction <- path.k.prediction[k,,]
    
    # The bandwidths to compute the bootstrap forecast densities are stored in: 
    bw.k.prediction <- array(0,dim=c(N))
    dimnames(bw.k.prediction)[[1]] <- nombres
    
    # The bootstrap forecast densities are stored in: 
    density.k.prediction <- array(0,dim=c(512,2,N))
    
    for (ser in 1:N)
        for (met.boot in 1:1) # only autoregressive
        {
            auxiliary <- density( k.prediction[,ser], bw=type.bw.forecast.dens) 
            bw.k.prediction[ser] <- auxiliary$bw
            density.k.prediction[,,ser] <- cbind(auxiliary$x,auxiliary$y)
        } 
    
    ##########################################################################################
    # 3. COMPUTING PAIRWISE DISTANCES MATRIX BASED ON PAIRS OF BOOTSTRAP FORECAST DENSITIES
    ##########################################################################################
    
    # The pairwise L1 distances are stored in: 
    DL1 <- array(0,dim=c(N,N))
    dimnames(DL1)[[1]] <- nombres;  dimnames(DL1)[[2]] <- nombres
    
    
    for (ser1 in 1:(N-1))
        for (ser2 in (ser1+1):N)
            for (met.boot in 1:1)
            {
                r <- range(density.k.prediction[,1,c(ser1,ser2)])
                a <- r[1] - 0.025*(r[2]-r[1]) ; b <- r[2] + 0.025*(r[2]-r[1])
                integrand_L1 <- function(u) 
                {
                    abs( estimator.density( density.k.prediction[,,ser1], u, bw.k.prediction[ser1] ) - estimator.density( density.k.prediction[,,ser2], u, bw.k.prediction[ser2] ) )
                }
                tryCatch( {
                    DL1[ser1,ser2] <- integrate( integrand_L1, lower=a, upper=b)$value }, error = function(e) {
                        plot.default( density.k.prediction[,,1], type="l", col="red", lwd=2, main="Error, problem intergrating the L1 distance of these densities, approximating...", xlab="", ylab="Density",
                                      xlim=c(min(density.k.prediction[,1,1], density.k.prediction[,1,2]), max(density.k.prediction[,1,1], density.k.prediction[,1,2])), ylim=c(0, max(density.k.prediction[,2,1], density.k.prediction[,2,2]) ) )
                        lines( density.k.prediction[,,2], col="blue", lwd=2 )
                        legend("topright", pch=16, col=c("red", "blue"), legend= c("x", "y") )
                        DL1[ser1,ser2] <- integrate( integrand_L1, lower=a, upper=b, stop.on.error=FALSE)$value
                    })
            }
    
    if (plot) {
        
        plot.default( density.k.prediction[,,1], type="l", col="red", lwd=2, main=paste("Prediction Density distance \nFor horizon:", k), xlab="", ylab="Density",
                      xlim=c(min(density.k.prediction[,1,1], density.k.prediction[,1,2]), max(density.k.prediction[,1,1], density.k.prediction[,1,2])), ylim=c(0, max(density.k.prediction[,2,1], density.k.prediction[,2,2]) ) )
        lines( density.k.prediction[,,2], col="blue", lwd=2 )
        legend("topright", pch=16, col=c("red", "blue"), legend= c("x", "y") )
    }
    
    list( L1dist = DL1[1,2], dens.x = density.k.prediction[,,1], dens.y = density.k.prediction[,,2] ,
          bx.x = bw.k.prediction[1], bw.y = bw.k.prediction[2])
}


 L1dist <- function( density.x, density.y, bw.x, bw.y) {
 r <- range(density.x[,1], density.y[,1])
 a <- r[1] - 0.025*(r[2]-r[1]) ; b <- r[2] + 0.025*(r[2]-r[1])
 integrand_L1 <- function(u) 
 {
     abs( estimator.density( density.x, u, bw.x ) - estimator.density( density.y, u, bw.y ) )
 }
 DL1 <- 2
 tryCatch( {
     DL1 <- integrate( integrand_L1, lower=a, upper=b)$value }, error = function(e) {
         plot.default( density.x, type="l", col="red", lwd=2, main="Error, problem intergrating the L1 distance of these densities, approximating...", xlab="", ylab="Density",
                       xlim=c(min(density.x[,1], density.y[,1]), max(density.x[,1], density.y[,1])), ylim=c(0, max(density.x[,2], density.y[,2]) ) )
         lines( density.y, col="blue", lwd=2 )
         legend("topright", pch=16, col=c("red", "blue"), legend= c("x", "y") )
         DL1 <- integrate( integrand_L1, lower=a, upper=b, stop.on.error=FALSE)$value
     })
    DL1
}

#prototyping multiple parameter per series distance
multidiss.PRED = function( series, h, B=500, logarithms=NULL, differences=NULL, plot=FALSE) {
    n <- length(series)
    distances <- matrix(0, n, n)
    rownames(distances) <- names(series)
    if (is.null(logarithms)) {
        logarithms <- rep(FALSE, n)
    }
    if (is.null(differences)) {
        differences <- rep(0, n)
    }
    densities <- list()
    #calc first all the individual densities by applying diss.PRED with the first
    individual.dens <- list()
    ii = 1
    while (ii < n) {
        dP <- diss.PRED(series[[ii]], series[[ii+1]], h , B, logarithms[ii], logarithms[ii+1] , differences[ii], differences[ii+1], FALSE )
        individual.dens[[ii]] <- list(dens=dP$dens.x, bw=dP$bw.x)
        individual.dens[[ii+1]] <- list(dens=dP$dens.y, bw=dP$bw.y)
        ii = ii + 2
        if (ii == n) ii = ii -1
    }
    for ( i in 1:(n-1) ) {
        for (j in (i+1):n ) {
            distance <- L1dist(individual.dens[[i]]$dens, individual.dens[[j]]$dens, individual.dens[[i]]$bw, individual.dens[[j]]$bw )
            distances[ i, j] <- distance
            distances[ j, i] <- distance
            densities[[i]] <- individual.dens[[i]]$dens
            densities[[j]] <- individual.dens[[j]]$dens
        }
    }
    #old way, less efficient
#     for ( i in 1:(n-1) ) {
#         for (j in (i+1):n ) {
#             distance <- diss.PRED(series[[i]], series[[j]], h , B, logarithms[i], logarithms[j] , differences[i], differences[j], FALSE )
#             distances[ i, j] <- distance$L1dist
#             distances[ j, i] <- distance$L1dist
#             densities[[i]] <- distance$dens.x
#             densities[[j]] <- distance$dens.y
#         }
#     }
    
    if (plot) {
        colors <- rainbow(length(densities))
        linetypes <- rep(1:6, ceiling(length(densities)/6))
        #linetypes = sample(linetypes)
        predxlim <- range(lapply( densities, function(x) { range(x[,1]) }))
        predylim <- range(lapply( densities, function(x) { range(x[,2]) }))
        plot( densities[[1]], type="l", col=colors[1], ylim=predylim, xlim=predxlim , xlab="", ylab="", main=paste("Prediction Densities \nAt horizon:", h), lty=linetypes[1], lwd=1)
        for ( i in 2:length(densities) ) {
            lines(densities[[i]], col=colors[i], lty=linetypes[i], lwd=1)
        }
        legend("topright", col=colors, legend= rownames(distances) , lty=linetypes, lwd=2, ncol=2 )
    }
    
    
    
    list( dist = as.dist((distances)), densities = densities )
}

