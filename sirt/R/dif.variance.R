 

#-------------------------------------------------------------------------#
# Routine for calculating DIF variance   (Camilli & Penfield, 1997)      #
dif.variance <- function( dif , se.dif , items = paste("item",1:length(dif),sep="") ){
        # calculating the weights
        wi <-  1/se.dif^2
        # mean DIF
        md <- mean(dif)
        # weighted tau2 estimator
        # formula (5) of PENFIELD & ALGINA (2006)
        weighted.tauq <- sum( wi^2 * ( dif - md )^2 - wi  ) / sum( wi^2 )
        # unweighted tau2 estimator
        # formula (5) of PENFIELD & ALGINA (2006)
        unweighted.tauq <- sum( ( dif - md )^2 - se.dif^2 ) / ( length(items) - 1 )
        # calculation of variances v_i
        vi <- weighted.tauq + se.dif^2
		weighted.tauq[ weighted.tauq < 0 ] <- 0			
		unweighted.tauq[ unweighted.tauq < 0 ] <- 0
        # Empirical Bayes DIF estimate
        lambda.i <- weighted.tauq / vi
        eb.dif <- lambda.i * ( dif - md ) + ( 1 - lambda.i) *md 
        list( weighted.DIFSD = sqrt(weighted.tauq) , 
				unweighted.DIFSD = sqrt(unweighted.tauq) ,
				mean.se.dif = sqrt( mean( se.dif^2 ) ) , eb.dif = eb.dif  )
        }
#-------------------------------------------------------------------------#


#################################################################################
dif.strata.variance <- function( dif , se.dif , itemcluster ){ 
    # stratified dif variance
    # means in differential item functioning
    # itemcluster is a vector of strata corresponding to items
    stratadif <- stats::aggregate( 1+0*dif , list( itemcluster ) , sum , na.rm = T )
    colnames(stratadif) <- c("strata" , "N.Items" )
    stratadif <- data.frame(stratadif)
    stratadif$M.DIF <- stats::aggregate( dif , list( itemcluster ) , mean , na.rm = T )[,2]
    # DIF in strata
    SS <- nrow(stratadif)
    for (ss in 1:SS){
        # ss <- 1
        items.ss <- which( itemcluster == stratadif[ss,"strata"]  )
        dif.ss <- dif[ items.ss ]
        difv.ss <- dif.variance( dif = dif.ss , se.dif = se.dif[ items.ss ] )
        stratadif$weighted.tau[ss] <- difv.ss$weighted.DIFSD
        stratadif$unweighted.tau[ss] <- difv.ss$unweighted.DIFSD
                }
        stratadif[ is.na(stratadif ) ] <- 0
    res <- list( "stratadif" = stratadif ,
                "weighted.DIFSD" = sum( stratadif$N.Items /  sum( stratadif$N.Items ) * 
							stratadif$weighted.tau ) ,
                "unweighted.DIFSD" = sum( 
					( stratadif$N.Items -1 )/  ( sum( stratadif$N.Items ) - 1) *  stratadif$unweighted.tau )
                        )
    return(res)
    }
#################################################################################

