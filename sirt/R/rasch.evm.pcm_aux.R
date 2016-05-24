
rasch.evm.pcm.dif <- function( b_evm , item , PARS_vcov , I , G , group.unique ,
		dat , dat.resp ){
   #********************************************************************
    difstats <- data.frame( "item" = colnames(dat) , "X2" = NA , "df" = NA , "p" = NA , "p.holm"=NA )  
   for (jj in 1:I){          
#        jj <- 4
        ii <- jj
        ind.ii <- which( item$item==ii )
        nj <- length(ind.ii)
        delta.jj <- rep(NA , nj*G )
        varmat.jj <- matrix(0 , nj*G , nj*G )
        Rdesign.jj <- matrix(0,nj*(G-1) , nj*G )
        for (gg in 1:G){
            delta.jj[  1:nj + nj*(gg-1)  ] <- b_evm[[gg]][ ind.ii ]           
            varmat.jj[ 1:nj + nj*(gg-1) , 1:nj + nj*(gg-1) ] <- PARS_vcov[[gg]][ ind.ii , ind.ii ]
            if (gg <G){
                for (vv in 1:nj){
                        Rdesign.jj[ vv + nj*(gg-1) , vv + nj*(gg-1) ] <- 1
                        Rdesign.jj[ vv + nj*(gg-1) , vv + nj*(gg) ] <- -1
                                }
                        }
                        }
                                               
        # calculate test value
        d0 <- Rdesign.jj %*% delta.jj
        # calculate variance matrix    
        ivm <- solve( Rdesign.jj %*% varmat.jj %*% t(Rdesign.jj ) )
        difstats[jj,"X2"] <- ( t(d0) %*% ivm %*% d0 )[1,1]
        difstats[jj,"df"] <- nrow(Rdesign.jj)
	
    # compute item-wise DIF statistics
    difjj <- stats::aggregate( delta.jj , list( rep(1:G , each=nj) ) ,  mean )[,2]
    for (gg1 in 1:(G-1)){
    for (gg2 in (gg1+1):(G)){ # gg1 <- 1 ; gg2 <- 2
        difstats[jj,paste0("DIF.",group.unique[gg1],".",group.unique[gg2]) ] <- difjj[gg1] - difjj[gg2]
                            }
        }
}
	difstats$p <- 1 - stats::pchisq( difstats$X2 , df=difstats$df )
	difstats$p.holm <- stats::p.adjust( difstats$p )               
#	difstats$V <- sqrt( difstats$X2 / ( colSums( dat.resp ) * difstats$df ) )	
	return(difstats)
}