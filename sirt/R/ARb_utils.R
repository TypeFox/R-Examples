


############################################################################
# This function reads fwf files                                            #
read.fwf2 <- function( file , format , variables = NULL){
    ff <- base::readLines( file )
    ind.ff1 <- c( 1, base::cumsum(format)[- length(format) ] + 1 )
    ind.ff2 <- base::cumsum(format)
    I <- length(format)
    n <- length( ff )
    dfr <- data.frame( matrix(0 , nrow= n , ncol=I ) )
    for (ii in 1:I){  
		dfr[,ii ] <- base::as.numeric( substring( ff , ind.ff1[ii] , ind.ff2[ii] )  ) 
					}
    if (!is.null(variables)){ 
		colnames(dfr) <- variables 
			}
    return(dfr)
    } 
############################################################################

#----------------------------------------------------------------
# utility function for formatting output in write.fwf2
.write.format2 <- function( vec1 , ff , fr ){
    if (fr == 0){
        vec2 <- round( vec1 , fr )
        blank.vv <- paste( rep( " " , ff ) , collapse="" )
        vec2 <- paste( substring( blank.vv , 1 , ff - nchar(vec2) ) , vec2 , sep="")    
            } else {
        d.vv <- round( vec1  , fr ) + 10^(-(fr+1))
        # generate blank
        blank.vv <- paste( rep( " " , ff+1 ) , collapse="" )
        d.vv <- paste( substring( blank.vv , 1 , ff+1 - nchar(d.vv) ) , d.vv , sep="")
        g.vv <- grep("NA" ,d.vv)
        d.vv[ g.vv  ] <- ifelse( ff > 1 ,  gsub( "NA" , " ." , d.vv[g.vv] ) , gsub( "NA" , "." , d.vv[g.vv] ) )
        vec2 <- substring( d.vv , 1 , ff )
        vec2
            }
    return(vec2 )
    }
#---------------------------------------------------------------

##############################################################################
write.fwf2 <- function( dat  , format.full , format.round , savename ){
        if (is.null( colnames(dat) ) ){ 
			colnames(dat) <- paste( "V" , 1:( ncol(dat) ) , sep="") 
							}
        matr <- matrix( " " , nrow= nrow(dat) , ncol = ncol(dat) )
        ind1 <- which( format.round <= 0  )
        format.full[ ind1 ] <- format.full[ind1] 
        format.round[ ind1 ] <- format.round[ind1]        
        for (vv in 1:( ncol(matr) ) ){
            fvv <- format.round[vv]
            fff <- format.full[vv]
            matr[,vv] <- .write.format2( vec1 = dat[,vv] , ff = fff , fr = fvv )
                }
        matr <- base::apply( matr , 1 , FUN = function(ll){ paste( ll , 
						collapse="" ) } )			
		if ( is.vector(matr) ){
			base::writeLines( matr , paste( savename , ".dat" , sep="") )
					} else {
			utils::write.table( matr , paste( savename , ".dat" , sep="") , 
						row.names=F , col.names=F)
						}
        dfr <- data.frame( "variable" = colnames(dat) , 
                    "begin" = c( 1 , cumsum( format.full )[ - ncol(dat) ] 
								+ 1 ) , 
                    "end" = cumsum( format.full )  ,
                    "length" = format.full
                            )
        utils::write.table( dfr , paste( savename , "__LEGEND.txt",sep="") , 
				quote=FALSE , row.names=FALSE , col.names=TRUE)
				
        return(dfr)
        }
##############################################################################
	
