
############################################################
# computation of Q3 statistic from residuals
.tam.q3.q3sub <- function( residM ){
    corM <- stats::cor(residM  , use="pairwise.complete.obs")
    I <- ncol(residM)
    dfr <- matrix( corM , ncol=1 )
    dfr <- data.frame( "index1" = rep(1:I,each=I) , "index2"= rep(1:I,I)  , "Q3" = dfr )
    dfr <- dfr[ dfr$index1 < dfr$index2 , ]
    dfr$aQ3 <- dfr$Q3 - mean( dfr$Q3)
    return(dfr)
        }
#############################################################

#############################################################
# Jackknife estimation with bias correction
.tam.q3.jackknife2 <- function( ms1 , ms.jack ){
        # pseudo values		
#		res0 <- jackknife_calc( ms1 , as.matrix(ms.jack ) )
		ms <- data.frame("val" = ms1 )
		JJ <- ncol(ms.jack)
        psx <- ms1 + ( JJ-1 ) * ( ms1 - ms.jack )		
        # jackknife estimate
        ms$jkunits <- JJ    
        ms$jk_est <- rowMeans( psx )

        ms$jk_se <- sqrt( rowSums( ( psx - ms$jk_est )^2 ) / (JJ-1 ) / JJ  )
#        ms$jk_se <- sqrt( apply( ms.jack , 1 , FUN = function(ll){ 
#                sum( ( ll - mean(ll) )^2 ) } ) * (JJ-1) / JJ )
        ms$est_low <- ms$jk_est - 1.96 * ms$jk_se
        ms$est_upp <- ms$jk_est + 1.96 * ms$jk_se
		return(ms)
		}
##############################################################


