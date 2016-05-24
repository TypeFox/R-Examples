personfit.stat <-
function( dat , abil , b ){
    dfr <- data.frame( "case" = seq( 1 , nrow(dat)) )
	dfr$NItems <- rowSums( 1 - is.na(dat)) 
    dfr$mean <- rowMeans( dat , na.rm=TRUE )	
    dfr$abil <- abil    
    # person fit statistics
    pval <- colMeans( dat , na.rm=TRUE )	
    dfr$caution <- pf.caution( data=dat , pval = pval )
    dfr$depend <- pf.dependability( data=dat , pval = pval )
    h1 <- pf.eci(data=dat, pval=pval, wle=abil, itemdiff=b)
    dfr <- cbind( dfr , h1$ECI )
    h1 <- pf.l0(data=dat, wle=abil, itemdiff=b)
    dfr$l0 <- h1$l0
    dfr$lz <- h1$lz
    h1 <- pf.outfit.infit(data=dat, wle=abil, itemdiff=b)
    dfr$outfit <- h1$Outfit
    dfr$infit <- h1$Infit
    dfr$rpbis <- pf.rpbis(data=dat, pval=pval)
    dfr$rpbis.itemdiff <- pf.rpbis.itemdiff(data=dat, itemdiff=b)
    dfr$U3 <- pf.U3( data=dat , pval = pval )
    return( dfr )
            }
