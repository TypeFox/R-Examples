micombine.cor <- function( mi.res , variables = 1:( ncol(mi.list[[1]]) ) ,  conf.level = .95 ,
		method="pearson"  ){
    # INPUT:
    # mi.res    ... MICE object
    # variables ... which variables are selected?
	if ( class(mi.res) == "mids.1chain"){
		mi.res <- mi.res$midsobj
			}	
	mi.list <- .milist( mi.res )		
    N <- nrow( mi.list[[1]] )
    VV <- length(variables)    
    # check if variables are given in character form
    if (is.character(variables)){ variables <- which( colnames(mi.list[[1]]) %in%  variables ) }
    dfr <- NULL
    for ( i in 1:(VV-1) ){
        for (j in (i+1):VV){
            if (i != j ){
            ii <- variables[i]
            jj <- variables[j]
                if ( i != j){ 
                    # calculate correlation coefficients
                    cor.ii.jj <- unlist( lapply( mi.list , FUN = function(dat){  
								stats::cor( dat[ , ii] , dat[,jj] , method=method) 
									} ) )
                    res.ii.jj <- .sub.micombine.cor( cor.list = cor.ii.jj , N = N , conf.level = conf.level )
                    dfr <- rbind( dfr , c( ii , jj , res.ii.jj ) )
            }   }}}
    vars <- colnames( mi.list[[1]] )
	dfr1 <- dfr
	dfr <- rbind( dfr , dfr1[ , c(2,1,seq(3,ncol(dfr) )) ] )
    if (VV == 2 ){  dfr <- dfr[ , -c(1:2) ] 
                    cat( paste( "Correlation of " , vars[ii] , " with " , vars[jj] )  , "\n" ) 
                    print( round( dfr[1,] , 4 ) )
        } else {
    dfr <- data.frame( "variable1" = vars[ dfr[,1] ] , "variable2" = vars[ dfr[,2] ] , dfr[ , -c(1:2) ] )
    print( data.frame( dfr[,1:2] , round( dfr[,-c(1:2)] , 4 ) ) )
            }
    invisible(dfr)
    }

	
	
	



#----------------------------------------------------------------------------------------------------#
# subroutine for combining correlations for multiply imputed data                                    #
.sub.micombine.cor <- function( cor.list , N , conf.level ){
        # convert correlations to Fisher transformed values
        fisher.cor.list <- as.list(1/2*log( ( 1 + cor.list) / ( 1 - cor.list ) ))
        var.fisher <- as.list( rep( 1/(N-3) , length(cor.list) ) )		
        # combination of point estimators according Rubin's formula
        fisher.cor.combine <- MIcombine( fisher.cor.list , var.fisher)		
        zr <- coef(fisher.cor.combine)
        zr.se <- sqrt( fisher.cor.combine$variance )[1,1]
        t.zr <- zr / zr.se
        fisher2cor <- function(z){ ( exp(2*z) - 1 )/ ( exp(2*z) + 1 ) }
        res <- c( "r" = fisher2cor(zr)  ,  
            "fisher_r" = zr ,
            "fisher_rse" = zr.se ,
			"fmi" = fisher.cor.combine$missinfo ,			
            "t" = t.zr  , 
            "p" = 2 * stats::pnorm( abs(t.zr) , lower.tail = FALSE ) ,
             fisher2cor( zr + stats::qnorm( ( 1 - conf.level ) / 2 ) * zr.se ) , 
             fisher2cor( zr - stats::qnorm( ( 1 - conf.level ) / 2 ) * zr.se ) )
            names(res)[7] <- paste( "lower" , round(100*conf.level,2),sep="")
            names(res)[8] <- paste( "upper" , round(100*conf.level,2),sep="")
        res <- c( res , - ( res[8] - res[7] ) / ( 2* stats::qnorm( ( 1 - conf.level )/2 ) ) )
        names(res)[9] <- "rse"
        res <- res[ c(1,9,2:8) ]
        return(res)
            }
#----------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------
# extract multiply imputed data sets from mice into a list of data frames
.milist <- function( mi.res ){ 
    mi.list <- NULL
    M <- mi.res$m   # extrahiere Anzahl der Imputationen
    for (ii in 1:M){ 
        mi.list[[ii]] <- mice::complete( mi.res , action= ii )
        }
    return(mi.list)
     }
#-------------------------------------------------------------------------


