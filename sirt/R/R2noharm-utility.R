



#..............................................................................
# utility functions for R2noharm

#---------------------------------------------------------------
# Read Tanaka Index                                             
.noharm.tanaka <- function( noharmout1 ){
    i1 <- grep( "Tanaka" , noharmout1 )
    tanaka <- as.numeric( strsplit( paste(noharmout1[i1]) , split = "=" )[[1]][2] )
    return(tanaka)
    }
#---------------------------------------------------------------

#-----------------------------------------------------------------
# Root mean square of residuals
.noharm.rmsr <- function( noharmout1 ){
    i1 <- grep( "Root mean square of residuals" , noharmout1 )
    RMSR <- as.numeric( strsplit( paste(noharmout1[i1]) , split = "=" )[[1]][2] )
    return(RMSR)
    }
#-----------------------------------------------------------------


#------------------------------------------------------------------
# UNIQUE VARIANCES, thresholds, 1-dimensional parametrizations    #
.noharm.itemlevel <- function( noharmout1 , itemlevel.type , I=I,dat=dat){
        i1 <- grep( itemlevel.type , noharmout1 )
        VV <- floor(I/9)+1 * ( I/9 != floor( I/9) )
        uniquevar <- noharmout1[ seq( i1+5 , i1 + 5 + 4*(VV-1) , 4 ) ]
		#*** addition ARb 2014-01-29
			uniquevar <- uniquevar[ substring( uniquevar , 1 , 1 ) == " " ]
		#***
        uniquevar <- paste( uniquevar , collapse= "  " )
        uniquevar <- strsplit( uniquevar , split = " " )[[1]]
        uniquevar <- as.numeric( paste( uniquevar[ uniquevar != "" ] ) )
		uniquevar <- uniquevar[1:I]
        names(uniquevar) <- colnames(dat)
        return( uniquevar )
        }
#------------------------------------------------------------------


#-------------------------------------------------------------------
# Function for extracting factor loadings                           
.noharm.loadings <- function( noharmout1 , type1 , type2 , dimensions=dimensions , I=I , dat = dat){
    i1 <- grep( type1 , noharmout1 )[1]
    i2 <- grep( type2 , noharmout1 )[1]
    rows <- seq( i1+2 , i2-4  )
    if ( ( dimensions == 1 ) & ( type1 == "Factor Loadings" ) ){
		rows <- seq(i1+2 , i2-6 )
			}
    r1 <- intersect( which( noharmout1 == "" ) , rows )
    VV <- length(r1)/2
    r2 <- stats::aggregate( r1 , list(rep( 1:VV , each=2 )) , mean )
    colnames(r2) <- c("index" , "index.row" )
    index.m <- data.frame( r2 , "begin" = r2[,2] + 2 , 
                    "end" = c( r2[,2][-1] - 2 , i2-4 ) )
		if ( dimensions == 1 & type1 == "Factor Loadings"){
				index.m$end <- i2 - 6
					}
    loading.matrix <- matrix(0,I,dimensions)
    rownames(loading.matrix) <- colnames(dat)
    colnames(loading.matrix) <- paste( "F" , 1:dimensions, sep="")
    for (vv in 1:VV){
        ind.vv.col <- noharmout1[ index.m$index.row[vv] ]
        ind.vv.col <- strsplit( ind.vv.col , split=" " )[[1]]
        ind.vv.col <- as.numeric( ind.vv.col[ ind.vv.col != "" ] )
        ind.vv <- seq( index.m$begin[vv]  , index.m$end[vv] )
        loading.vv <- noharmout1[ ind.vv ]
        loading.vv <- sapply( loading.vv , FUN = function(ll){
                            ll1 <- strsplit( paste(ll) , split=" " )[[1]] 
                            as.numeric( ll1[ ll1 != ll1[1] ] )
                                                } )
		loading.vv <- as.matrix( loading.vv )
        ind.vv.row <- as.vector( loading.vv[1,] )
        loading.vv <-  loading.vv[-1,]
        loading.matrix[ ind.vv.row , ind.vv.col ] <- t(loading.vv) 
        }  
        return( loading.matrix )
        }
#--------------------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Function for extracting Factor Correlations                              
.noharm.correlations <- function( noharmout1 , 
						type1 = "Factor Correlations" , 
						type2 = "Residual Matrix" ,
                                    dimensions=dimensions , dat = dat ){
    i1 <- grep( type1, noharmout1 )[1]
    i2 <- grep( type2 , noharmout1 )[1]
    rows <- seq( i1+2 , i2-4 )
    r1 <- intersect( which( noharmout1 == "" ) , rows )
    VV <- length(r1)/2
    r2 <- stats::aggregate( r1 , list(rep( 1:VV , each=2 )) , mean )
    colnames(r2) <- c("index" , "index.row" )
    index.m <- data.frame( r2 , "begin" = r2[,2] + 2 , 
                    "end" = c( r2[,2][-1] - 2 , i2-4 ) )
    correlation.matrix <- matrix(NA,dimensions,dimensions)
    rownames(correlation.matrix) <- colnames(correlation.matrix) <- paste( "F" , 1:dimensions, sep="")
    for (vv in 1:VV){
        ind.vv.col <- noharmout1[ index.m$index.row[vv] ]
        ind.vv.col <- strsplit( ind.vv.col , split=" " )[[1]]
        ind.vv.col <- as.numeric( ind.vv.col[ ind.vv.col != "" ] )
        ind.vv <- seq( index.m$begin[vv]  , index.m$end[vv] )
        correlation.vv <- noharmout1[ ind.vv ]
        correlation.vv <- sapply( correlation.vv , FUN = function(ll){
                            ll1 <- strsplit( paste(ll) , split=" " )[[1]] 
                            as.vector( as.numeric( ll1[ ll1 != "" ] ))
                                                } )
        if ( is.list( correlation.vv ) ){ 
            ind.vv.row <- as.vector(unlist( lapply( correlation.vv , FUN = function(ll){ ll[1] } ) ))
            correlation.vv <-  lapply( correlation.vv , FUN = function(ll){ ll[-1] } ) 
            } else {
            correlation.vv <- as.vector( correlation.vv )
            ind.vv.row <- correlation.vv[1]
            correlation.vv <- as.list( correlation.vv[-1] )
            }
        RVV <- length(correlation.vv)
        for (rvv in 1:RVV){
            correlation.vv.rvv <- correlation.vv[[rvv]]
            i1.rvv <- ind.vv.row[rvv]
            i2.rvv <- ind.vv.col[ seq( 1 , length(correlation.vv.rvv)) ]              
            correlation.matrix[ i2.rvv , i1.rvv ] <- correlation.matrix[ i1.rvv , i2.rvv ] <- correlation.vv.rvv
                }
        }  
    correlation.matrix 
        }
#---------------------------------------------------------------------------------------------
        


#--------------------------------------------------------------------------#
# Extracting Residual Matrix                                               #
.noharm.residuals <- function(noharmout1, I = I , dat = dat ){
    i1 <- grep( "Residual Matrix" , noharmout1 )
    i2 <- grep( "Sum of squares of residuals" , noharmout1 )
    rows <- seq( i1+2 , i2-4 )
    r1 <- intersect( which( noharmout1 == "" ) , rows )
    VV <- length(r1)/2
    r2 <- stats::aggregate( r1 , list(rep( 1:VV , each=2 )) , mean )
    colnames(r2) <- c("index" , "index.row" )
    index.m <- data.frame( r2 , "begin" = r2[,2] + 2 , 
#                    "end" = c( r2[,2][-1] - 2 , i2-4 ) )
                    "end" = c( r2[,2][-1] - 2 , i2-3 ) )
    resid.matrix <- matrix(0,I,I)
    rownames(resid.matrix) <- colnames(resid.matrix) <- colnames(dat)
    for (vv in 1:VV){
        ind.vv.col <- noharmout1[ index.m$index.row[vv] ]
        ind.vv.col <- strsplit( ind.vv.col , split=" " )[[1]]
        ind.vv.col <- as.numeric( ind.vv.col[ ind.vv.col != "" ] )
        ind.vv <- seq( index.m$begin[vv]  , index.m$end[vv] )
        resid.vv <- noharmout1[ ind.vv ]
        resid.vv <- sapply( resid.vv , FUN = function(ll){
                            ll1 <- strsplit( paste(ll) , split=" " )[[1]] 
                            as.numeric( ll1[ ll1 != "" ] )
                                                } )								
    if (is.list( resid.vv)){ 
			gh1 <- lapply( resid.vv , FUN = function(ll){ length(ll) } )
			ind1 <- which (gh1 == 0)
			if ( length(ind1) > 0 ){
				for (ii in ind1 ){	resid.vv[[ii]] <- NULL }
						}
            ind.vv.row <- as.vector(unlist( lapply( resid.vv , FUN = function(ll){ ll[1] } ) ))
            resid.vv <-  lapply( resid.vv , FUN = function(ll){ ll[-1] } ) 
                        } else {
                        resid.vv <- as.vector(resid.vv)
                        ind.vv.row <- resid.vv[1]
                        resid.vv <- as.list( resid.vv[-1] )
                        }
        RVV <- length(resid.vv)
        for (rvv in 1:RVV){
            resid.vv.rvv <- resid.vv[[rvv]]
            i1.rvv <- ind.vv.row[rvv]
            i2.rvv <- ind.vv.col[ seq( 1 , length(resid.vv.rvv)) ]              
            resid.matrix[ i2.rvv , i1.rvv ] <- resid.matrix[ i1.rvv , i2.rvv ] <- resid.vv.rvv
                }
        }  
        return( resid.matrix )
        }
#-------------------------------------------------------------------------------------
