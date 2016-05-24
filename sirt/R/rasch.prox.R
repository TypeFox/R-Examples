


#------------------------------------------------------------------------------------------------------------------------------
# PROX routine for Rasch models
rasch.prox <- function( dat , dat.resp = 1 - is.na(dat) , 
				freq=rep(1,nrow(dat)) , 
				conv = 0.001 , maxiter = 30 , progress = FALSE ){
        # INPUT:
        # dat       ... data frame
        # dat.resp  ... missing response pattern
        # freq      ... frequencies of item response pattern
        # conv      ... convergence criterion (maximal change for item parameter estimates)
        # maxiter   ... maximal number of PROX iterations
        # print     ... should estimation process being displayed       
        I <- ncol( dat )
        IP <- nrow(dat )
        # calculate logit frequency for every item
        s.i <- colSums( dat * dat.resp * freq )
        n.i <- colSums( dat.resp * freq )
        logit.freq.i <- stats::qlogis( ( s.i / n.i + .01)/1.02 )
        # calculate logit frequency for every person (response pattern)
        r.n <- rowSums( dat * dat.resp * freq )
        n.n <- rowSums( dat.resp * freq )
        logit.freq.n <- stats::qlogis( ( r.n / n.n + .01 )/1.02 )
        # initial mean and SD of logit abilities for persons which answers item i 
        d.i <- mu.i <- rep(0,I)
        sigma.i <- rep(1,I)       
        # initial distribution of logit item diffuculties encountered by person n (respective pattern n)
        b.n <- mu.n <- rep( 0 , IP )
        sigma.n <- rep( 1, IP )
        # Number of subjects per item
        N.i <- colSums( dat.resp * freq )
        # Number of items per subject pattern
        N.n <- rowSums( dat.resp  )
        iter <- 0
        par.change <- 1
		Mdf <- dat.resp * freq
        # PROX algorithm                            #
        while ( ( par.change > conv ) & ( iter < maxiter )) { 
                d.i0 <- d.i
                b.n0 <- b.n
                # update item difficulty
                d.i <- mu.i - sqrt( ( 1 + sigma.i^2 / 2.9 ) ) * logit.freq.i        
                # update ability estimate for pattern
                b.n <- mu.n + sqrt( ( 1 + sigma.n^2 / 2.9 ) ) * logit.freq.n 
                # center persons
                b.n <- b.n - stats::weighted.mean( b.n , freq )
                # update mu.i and sigma.i 
                # (mean and standard deviation of the logit abilities of the person encountering item i)
                mu.i <- colSums( Mdf * b.n ) / N.i
				g1 <- ( colSums( b.n^2 * Mdf ) - N.i * mu.i^2  ) / ( N.i - 1 )
				sigma.i <- sqrt( ifelse( g1 < 0 , .0001 , g1 ) )
                # update mu.n and sigma.n
                # (mean and standard deviation of the logit difficulties encountered by person n)
#                d.i.m <- outer( rep(1 , IP ) , d.i )
				d.i.m <- matrix( d.i , nrow=IP , ncol=I , byrow=T )
                mu.n <- rowSums( dat.resp *  d.i.m ) / N.n
				g1 <- ( rowSums( d.i.m^2 * dat.resp ) - N.n * mu.n^2  ) / ( N.n - 1 )
                sigma.n <- sqrt( ifelse( g1 < 0 , .0001 , g1 ) )
                iter <- iter + 1
                par.change <- max( abs(d.i - d.i0) )
                if ( progress){
                    cat( "PROX Iter." , iter , ": max. parm. change = " , round( par.change , 6 ) , "\n")
                    utils::flush.console()
                        }
        }    
    return( list( "b" = d.i , "theta" = b.n , "iter" = iter ,
					"sigma.i" = sigma.i , "sigma.n" = sigma.n ) )
    } 
#------------------------------------------------------------------------------------------------------------------------------

