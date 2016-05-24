 

#-----------------------------------------------------
# nonparametric estimation of conditional covariance
ccov.np <- function( data , score , bwscale = 1.1 , thetagrid = seq( -3,3,len=200) , 
		progress=TRUE ){
            # number of Items I
            I <- ncol(data)
            # z-standardization of score
            score <- scale( score )[,1]
            # matrix of item response functions
            if (progress ){
                cat("Pairwise Estimation of Conditional Covariances\n" )            
                cat("...........................................................\n" )
                cat("Nonparametric ICC estimation \n " ) }
            icc.items <- matrix( 0 , length(thetagrid) , I )
            if ( I >= 20 ){ 
                    display <- seq( 1 , I , floor( I/20  ) )[ 2:20 ]
                        } else { display <- 20 }
            i <- 1
            for ( ii in 1:I ){
                x <- score[ ! is.na( data[,ii] )  ]
                y <- data[ ! is.na( data[,ii] ) , ii ]
                icc.items[,ii] <- stats::ksmooth( x  , y , bandwidth = bwscale * length(x)^(-1/5)  , 
							x.points = thetagrid , kernel="normal")$y
                if ( i < 20 ){ if ( ii == display[i] & progress ){ 
							cat( paste( 5*i  , "% " , sep="" ) ) ; i <- i + 1 ; 
                                                        if (i == 11){ cat("\n" ) }
                                                        utils::flush.console()} }
                }
            if ( progress){ cat("\n") }
            # weights thetagrid
            wgt.thetagrid <- stats::dnorm(thetagrid)
            wgt.thetagrid <- wgt.thetagrid 
            if (progress ){
                cat("...........................................................\n" )
                cat("Nonparametric Estimation of conditional covariances \n " ) 
                utils::flush.console()
                    }
            # calculation of conditional covariance
            ccov.table <- data.frame( "item1ID" = rep( 1:I , I ) , "item2ID" = rep( 1:I , each = I ) )
            ccov.table <- ccov.table[ ccov.table$item1ID < ccov.table$item2ID , ]
            ccov.table$N <- apply( ccov.table , 1 , FUN = function(ll){ 
                        sum( rowSums( is.na( data[ , c( ll[1] , ll[2] ) ] ) ) == 0 ) } )
            ccov.table <- ccov.table[ ccov.table$N > 0 , ]
            ccov.table$item1 <- colnames(data)[ ccov.table$item1ID ]
            ccov.table$item2 <- colnames(data)[ ccov.table$item2ID ]
            ccov.table$itempair <- paste( ccov.table$item1 , ccov.table$item2 , sep="-" )
            # smoothing all item pairs
            # calculate conditional covariances
            FF <- nrow( ccov.table )
            ccor.matrix <- ccov.matrix <- prod.matrix <- matrix( 0 , nrow= length(thetagrid ) , ncol = FF )
            ii <- 1
            for (ff in 1:FF){
                display <- seq( 1 , FF , floor( FF/20  ) )[ 2:20 ]
                data.ff <- data[ , c( ccov.table[ff,1] , ccov.table[ff,2] ) ]
                which.ff <- which( rowSums( is.na( data.ff ) ) == 0  )
                data.ff <- data.ff[ which.ff , ]
#                y <- data[ ! is.na( data[,ii] ) , ii ]     #       Bug: 2011-07-20
                prod.matrix[,ff] <- stats::ksmooth( x = score[ which.ff]  , y = data.ff[,1]*data.ff[,2] , 
                        bandwidth = bwscale * length(which.ff)^(-1/5)  , x.points = thetagrid , kernel="normal")$y
                ccov.matrix[ , ff ] <- prod.matrix[,ff] -  icc.items[, ccov.table[ff,1] ] * icc.items[, ccov.table[ff,2] ]
                if ( ii < 20 ){ if ( ff == display[ii] & progress==T ){ cat( paste( 5*ii  , "% " , sep="" ) ) ; ii <- ii + 1 ; 
				utils::flush.console()
                                                            if (ii == 11){ cat("\n" ) }
                            } }
                }
            # remove NAs from ccov.matrix
            ccov.matrix[ is.na( ccov.matrix) ] <- 0
            if ( progress ){ cat("\n") }
            # calculate (weighted) conditional covariance
            ccov.table$ccov <- apply( ccov.matrix , 2 , FUN = function(sp){ 
						stats::weighted.mean( sp  , wgt.thetagrid ) } )
            res <- list( "ccov.table" = ccov.table , "ccov.matrix" = ccov.matrix ,
                            "data" = data , "score" = score , "icc.items" = icc.items )
            return( res ) 
            }
#------------------------------------------------------------------------------------------------------------------------------

##NS export(detect.index)
#-----------------------------------------------------------------------------------------------------
detect.index <- function( ccovtable , itemcluster ){
    # INPUT:
    # result from ccov.np 
    # itemcluster ... identifies an item cluster for each item
    #.............................
    # calculate delta
    ccovtable <- ccovtable$ccov.table
    ccovtable$delta <- ifelse( itemcluster[ ccovtable$item1ID ] == itemcluster[ ccovtable$item2ID ] , 1 , -1 )
    # calculate indizes
    indizes <- c( 100*mean( ccovtable$ccov * ccovtable$delta ) , 
                mean( sign( ccovtable$ccov ) * ccovtable$delta ) , 
                    sum( ccovtable$ccov  * ccovtable$delta ) / sum( abs( ccovtable$ccov ) ) )
	# calculate weighted indizes 
    weighted.indizes <- c( 100* stats::weighted.mean( ccovtable$ccov * ccovtable$delta , sqrt(ccovtable$N) ) ,
        stats::weighted.mean( sign( ccovtable$ccov ) * ccovtable$delta , sqrt(ccovtable$N) ) , 
        sum( ccovtable$ccov  * ccovtable$delta * sqrt(ccovtable$N) ) / sum( abs( ccovtable$ccov ) * sqrt(ccovtable$N) ) )
    res <- data.frame( "unweighted" = indizes , "weighted" = weighted.indizes )
    rownames(res) <- c("DETECT" , "ASSI" , "RATIO" )
    res
    }
#-----------------------------------------------------------------------------------------------------




##NS export(conf.detect)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Confirmatory DETECT analysis
conf.detect <- function( data , score , itemcluster , bwscale = 1.1 , progress = TRUE ,
                            thetagrid = seq( -3,3,len=200)  ){
    cat("-----------------------------------------------------------\n" )
    cat("Confirmatory DETECT Analysis \n" ) ; flush.console()
    h1 <- is.matrix( score )
    if (h1 ){ PP <- ncol(score) }
    if (! h1  ){  cat("Conditioning on 1 Score\n" )  } else {
            cat(paste("Conditioning on ",PP, " Scores\n" , sep="") ) }
    cat(paste("Bandwidth Scale:" , bwscale , "\n" ) ) 
    utils::flush.console()
    if ( ! h1 ){
        ccovtable <- ccov.np( data , score = score, bwscale = bwscale , 
                                progress= progress , thetagrid = thetagrid )
        res <- detect.index( ccovtable , itemcluster = itemcluster )
                    } else {
            ccovtable.list <- list()
            for (pp in 1:PP){
                cat( paste( "DETECT Calculation Score " , pp , "\n" , sep="") ) ; 
				utils::flush.console()
                ccovtable.list[[pp]] <- ccov.np( data , score = score[,pp], 
                                    bwscale = bwscale , progress= FALSE )
                    }  
        detect.list <- lapply( ccovtable.list , FUN = function( ccovtable ){ 
                    detect.index( ccovtable , itemcluster=itemcluster ) } )
        detect.matrix <- matrix( unlist( lapply( detect.list , FUN = function( ll){ c( ll[1,] , ll[2,] , ll[3,] ) } ) ) , nrow=PP , byrow=T)
        detect.summary <- data.frame( "NScores" = PP , "Mean" = colMeans( detect.matrix ) , 
                    "SD" = apply( detect.matrix , 2 , stats::sd ) , 
                    "Min" = apply( detect.matrix , 2 , min ) , 
                    "Max" = apply( detect.matrix , 2 , max ) 
                    )
        rownames(detect.summary) <- c("DETECT Unweighted" , "DETECT Weighted" , "ASSI Unweighted" , "ASSI Weighted" ,
                        "RATIO Unweighted" , "RATIO Weighted" )
            }
    cat("-----------------------------------------------------------\n" )            
    if ( ! h1){   res <- list(  "detect" = res , "ccovtable" = ccovtable , "detect.summary" = res ) } else
            {     res <- list(  "detect" = detect.list , "ccovtable" = ccovtable.list , "detect.summary" = detect.summary ) }
    print(round(res$detect.summary,3))
    return( res )
    }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##NS export(expl.detect)
###############################################################################
# Exploratory DETECT analysis
expl.detect <- function( data , score , nclusters , N.est = NULL , seed=897 , bwscale = 1.1 ){
	set.seed( seed )
    # number of items
    I <- ncol(data)
    # sample for estimation
    N <- nrow(data)
    if ( is.null( N.est ) ){ N.est <- floor(N/2) }
    estsample <- sort( sample( 1:N , floor( N.est ) ) )
    # validation sample
    valsample <- setdiff( 1:N , estsample )
    #**********************************
    # Maximizing DETECT index
    #**********************************
    # nonparametric estimation of conditional covariance
    cc <- ccov.np( data = data[ estsample,] , score = score[estsample] , bwscale = bwscale )
    ccov.matrix <- .create.ccov( cc , data = data[ estsample,]  )
    # create distance matrix
    cc1 <- max(ccov.matrix) - ccov.matrix
    # Ward Hierarchical Clustering
    d <- stats::as.dist(cc1)
    fit <- stats::hclust(d, method="ward")         # hierarchical cluster analysis
    clusterfit <- fit
    itemcluster <- data.frame( matrix( 0 , I , nclusters ) )
    itemcluster[,1] <- colnames(data)
    colnames(itemcluster) <- c( "item" , paste( "cluster" , 2:nclusters , sep="") )
    detect.unweighted <- detect.weighted <- NULL
    for (k in 2:nclusters){ 
        itemcluster[,k] <- stats::cutree( fit , k=k ) 
        h1 <- detect.index( ccovtable=cc , itemcluster = itemcluster[,k] )
        detect.unweighted <- rbind( detect.unweighted , h1$unweighted )
        detect.weighted <- rbind( detect.weighted , h1$weighted )    
                }
    colnames(detect.unweighted) <- paste( c( "DETECT" , "ASSI" , "RATIO" ) , ".est" , sep="")
    colnames(detect.weighted) <- paste( c( "DETECT" , "ASSI" , "RATIO" ) , ".est" , sep="")
    dfr1 <- data.frame( "N.Cluster" = 2:nclusters )
    dfr1$N.items <- I
    dfr1$N.est <- N.est
    dfr1$N.val <- length(valsample)
    dfr1$size.cluster <- sapply( 2:nclusters  , FUN = function(tt){ paste( table( itemcluster[,tt] ) , collapse="-" ) } )
    detu <- data.frame( dfr1 , detect.unweighted )
    detw <- data.frame( dfr1 , detect.weighted )
    #************************************
    # Validating DETECT index
    #************************************
	if ( length(valsample) > 0 ){
		cc <- ccov.np( data = data[ valsample,] , score = score[valsample] , bwscale = bwscale )
		detect.unweighted <- detect.weighted <- NULL
		for (k in 2:nclusters){ 
			h1 <- detect.index( ccovtable=cc , itemcluster = itemcluster[,k] )
			detect.unweighted <- rbind( detect.unweighted , h1$unweighted )
			detect.weighted <- rbind( detect.weighted , h1$weighted )    
					}
		colnames(detect.unweighted) <- paste( c( "DETECT" , "ASSI" , "RATIO" ) , ".val" , sep="")
		colnames(detect.weighted) <- paste( c( "DETECT" , "ASSI" , "RATIO" ) , ".val" , sep="")
		detu <- data.frame( detu , detect.unweighted )
		detw <- data.frame( detw , detect.weighted )
			}
	cat("\n\nDETECT (unweighted)\n\n")
	clopt <- which.max( detu$DETECT.est ) + 1 
	cat("Optimal Cluster Size is " , clopt , " (Maximum of DETECT Index)\n\n" )
	detu1 <- detu
	for (vv in 6:ncol(detu)){ detu1[,vv] <- round( detu1[,vv] , 3) }
    print(detu1)
    res <- list( "detect.unweighted" = detect.unweighted , "detect.weighted" = detect.weighted ,
                    "clusterfit" = clusterfit , "itemcluster" = itemcluster )
	# plot cluster solution
	graphics::plot( res$clusterfit , 
			main = paste( "Cluster Dendogram with " , clopt , " Clusters" , sep="")
					)
    stats::rect.hclust(res$clusterfit, k=clopt, border="red")
	return(res)
    }
###############################################################################

##NS # export(.create.cov)
#*********************************************************
# auxiliary function for creating a covariance matrix
.create.ccov <- function( cc , data ){
    ccc <- cc$ccov.table
    I <- max( ccc$item1ID , ccc$item2ID )
    ccov.matrix <- matrix( 0 , I , I)
    rownames(ccov.matrix) <- colnames(ccov.matrix) <- colnames(data)
    LL <- nrow(ccc)
    for (ll in 1:LL){ 
            ccov.matrix[ ccc$item1ID[ll] , ccc$item2ID[ll] ] <- ccc$ccov[ll]
            ccov.matrix[ ccc$item2ID[ll] , ccc$item1ID[ll] ] <- ccov.matrix[ ccc$item1ID[ll] , ccc$item2ID[ll] ]
                        }
    return( ccov.matrix) 
        }
#*********************************************************
