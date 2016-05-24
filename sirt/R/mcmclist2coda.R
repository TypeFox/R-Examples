
##############################################################
# write elements from mcmcmlist into code file
mcmclist2coda <- function( mcmclist , name , coda.digits=5 ){
    m1 <- mcmclist[[1]]
    vars <- colnames(m1)
    # create codaIndex file
    BB <- nrow(m1)
    VV <- length(vars)
    c1 <- paste( vars , seq( 1 , BB*VV , BB ) , seq( BB , BB*VV , BB ) )
    base::writeLines( c1 , paste0( name , "_codaIndex.txt" ) )
    # create coda file
    m2 <- matrix( m1 , ncol=1 )
    m2 <- paste( rep(1:BB , VV ) , round( m2[,1] , coda.digits )  )
    base::writeLines( m2 , paste0( name , "_coda1.txt" ) )
            }
#########################################################
