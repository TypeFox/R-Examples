
#--------------------------------------------------------------------------
# -----------------------------------------------------------
# data preparations for rasch.jml and rasch.mml
.data.prep <- function( dat , weights = NULL , use.freqpatt = TRUE ){
        #-------------------
        # freq.patt ... should frequency pattern be taken into account?
        #               default = TRUE
        #-------------------
        # should items being excluded?
# a0 <- Sys.time()		
        item.means <- base::colMeans( dat , na.rm=T )
        item.elim <- which( item.means %in% c(0,1))
        if ( length( item.elim ) > 0 ){
                stop( cat( paste( "There are" , length(item.elim) , "Items with no variance!") ) ) 
                    }
        if ( any( is.na( item.means )) ){ stop( "There are items which contains only missings!") }
         n <- nrow(dat)
        I <- ncol(dat)  
		if( is.null(weights) ){  weights <- rep( 1 , n ) }
        # indicator for nonmissing response
        dat.9 <- dat
        dat.9[ is.na(dat) ] <- 9
         # pattern                           
        if ( use.freqpatt ){ 		#		
			freq.patt <- apply( dat.9 , 1 , FUN = function(ll){ paste(ll , collapse="" ) } )  #
#		freq.patt <- dat.9[,1]
#		for (ii in 2:I){ freq.patt <- paste0( freq.patt , dat.9[,ii] ) }	
        # frequency pattern with frequencies
					dat1 <- data.frame( table( freq.patt ) )
# cat("a220") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1										
                        } else { 
					freq.patt <- paste("FP" , 1000000 + 1:n , sep="")
					dat1 <- data.frame( freq.patt )
					colnames(dat1)[1] <- "freq.patt"
							}
         # weighting the frequencies if survey weights are supplied
        if( !is.null(weights) ){        
            # standardize weights
            weights <- weights / sum( weights ) * n 
            if ( use.freqpatt ){ 
                    dat1[,2] <- stats::aggregate( weights , list( freq.patt) , sum )[,2]
                        } else { 
					dat1[,"Freq"] <- weights 
						}
                    }	
         # item pattern corresponding to frequency pattern
		if ( use.freqpatt){ 
			dat2 <- matrix( as.numeric( unlist( strsplit( paste(dat1[,1]) , "" ) ) ) , ncol= ncol(dat) , byrow=T)
					} else { 
					dat2 <- dat.9 }
        dat2.resp <- 1 * ( dat2 !=9 )
		dat2[ dat2 == 9 ] <- 0
         # mean right
        dat1$mean <- rowSums( dat2 * dat2.resp )  / rowSums( dat2.resp )     
        freq.patt <- data.frame(  freq.patt , rowMeans( dat , na.rm=TRUE ) , 1:n )
        colnames(freq.patt)[2:3] <- c("mean" , "index" )
        list( "dat" = dat , "dat2" = dat2 , "dat2.resp" = dat2.resp , "dat1" = dat1 , 
				"freq.patt" = freq.patt  , "I" = I , "n" = n ,
                "dat9" = dat.9  )       
        }
    #*******************
    #   OUTPUT:
    #   dat         ... original data
    #   dat2        ... reduced original data. Each different item resonse is represented by one row.
    #   dat2.resp   ... indicator response matrix for each item response pattern
    #   dat1        ...     --- ADD AN EXPLANATION HERE     ----
    #   freq.patt   ... absolute frequency of each item response pattern
    #   I           ... number of items
    #   n           ... number of subjects
    #   dat9        ... This is the original data exact from the fact that missings are recoded by 9.
#-----------------------------------------------------------------
# -------------------------------------------------------------------

#--------------------------------------------------------
# Small function which helps for printing purposes
.prnum <- function( matr , digits ){
        VV <- ncol(matr)
        for (vv in 1:VV){
        # vv <- 1
        if ( is.numeric( matr[,vv]) ){ matr[,vv] <- round( matr[,vv] , digits ) }
                        }
        print(matr)
                   }
#--------------------------------------------------------



############################################################
# Function for calculation of a response pattern 
#   for dichotomous responses
resp.pattern2 <- function( x ){
    n <- nrow(x)
    p <- ncol(x)
    mdp <- (x %*% (2^((1:ncol(x)) - 1))) + 1
    misspattern <- mdp[,1]
    misspattern <- list( "miss.pattern" = mdp[,1] , 
                "mp.index" = match( mdp[,1] , sort( unique(mdp[,1] ) ) ) )
    return( misspattern )
        }   
########################################################
