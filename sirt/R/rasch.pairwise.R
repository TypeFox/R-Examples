

#---------------------------------------------------------------#
# Rasch estimation (Approximate Method)                         #
# also called MINCHI method                                     #
# Handbook of Statistics Vol. 26                                #
# Chapter of G. Fischer: p. 544 ff.                             #
# pairwise likelihood method
##NS export(rasch.pairwise)
rasch.pairwise <- function( dat , conv = 0.0001 , maxiter = 3000 , 
	progress = TRUE , b.init = NULL , zerosum = FALSE ){
	#******************************************************#
    # INPUT:                                               #
    # dat      ... data frame                              #
    # conv     ... convergence in epsilon parameters       #
    # maxiter  ... maximal number of iterations            #
    # progress ... display progress in MINCHI method       #
    # beta.init ... initial beta values                    #
    s1 <- Sys.time()
	CALL <- match.call()
        # should items being excluded?
        item.elim <- which( colMeans( dat , na.rm=T ) %in% c(0,1))
        if (length(item.elim)>0){ dat <- dat[ , - item.elim ] }

        dat <- as.matrix(dat)
        # data preparation
        dat.resp <- 1 - is.na( dat )
        dat.9 <- dat
        dat.9[ is.na(dat) ] <- 9
        # calculate n_{ij}
        n.ij <- t( dat.9 * dat.resp ) %*%  ( ( 1 - dat.9 ) * dat.resp  )
        # which item pairs occur in estimation procedure
        delta.ij <- 1 * ( n.ij + t( n.ij ) > 0 ) 

        # initial values for beta
        if( is.null( b.init) ){ beta <- - stats::qlogis( colMeans( dat , na.rm=T ) ) } else { beta <- b.init } 
        # calculate y_{ij} values
        y.ij <- n.ij / ( n.ij + t( n.ij) )
        y.ij[ delta.ij == 0 ] <- 0
        y.ji <- t( y.ij )

        eps <- exp( - beta )
        change <- 1
        iter <- 0
        while( change > conv & iter < maxiter ){
                eps0 <- eps
                eps <- sqrt( rowSums( y.ij * eps * delta.ij ) / colSums( y.ij * 1 / eps ) )
				if (zerosum){
					b1 <- - log(eps)
					b2 <- b1 - mean(b1)
					eps <- exp( - b2 )	
								}
                change <- max( abs( eps0 - eps ) ) 
                iter <- iter + 1
                if ( progress ){
                    cat( "PL Iter." , iter , ": max. parm. change = " , 
                                round( max(abs( log(eps0) - log(eps))) , 6 ) , "\n")
                    utils::flush.console()
                        }                
                }
        item <- data.frame( "N" = colSums(1 -is.na(dat)) , "p" = colMeans( dat , na.rm=T ) , 
#                                "itemdiff" = scale( - log(eps) , scale=F ) 
						"b" = - log(eps) ,
						"itemcluster"= rep(0,ncol(dat))
									)
        s2 <- Sys.time()									
        res <- list( "b" = - log( eps ) , "eps" = eps , "iter" = iter , "conv" = conv , "dat" = dat ,
                    "freq.ij" = n.ij  , "item" = item , "fct" = "rasch.pairwise",
					"s1"=s1 , "s2"=s2 , CALL = CALL ) 
        class(res) <- "rasch.pairwise"
        return(res)
    }    
#-------------------------------------------------------------------




#**********************************************************
# Summary for rasch.minchi object                         *
##NS S3method(summary,rasch.pairwise)
summary.rasch.pairwise <- function( object , ...){
    cat("------------------------------------------- \n")
    d1 <- utils::packageDescription("sirt")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
	cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
	cat("Computation Time:" , print(object$s2 - object$s1), "\n\n")
    cat(paste0("Function '" , object$fct , "'") , "\n") 
	
	cat("Call:\n", paste(deparse(object$CALL), sep = "\n", collapse = "\n"), 
				"\n\n", sep = "")				
    cat("------------------------------------------- \n")
    cat("Pairwise likelihood estimation \n")
    cat("Rasch Model \n")
    cat("------------------------------------------- \n")
    cat("Item Parameters \n")
    print( round( object$item , 3 ))                
                }
#*******************************************************



