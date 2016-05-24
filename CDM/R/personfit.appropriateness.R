
########################################################
# personfit.appropriateness statistics
########################################################

########################################################################################
# function personfit appropriateness statistic
personfit.appropriateness <- function( data , probs , skillclassprobs , h = .001 ,
        eps = 1E-10 , maxiter = 30 , conv = 1E-5 , max.increment = .1 , progress=TRUE ){
    # data input
    data.resp <- 1 - is.na( data )
    data[ is.na(data) ] <- 0
    N <- nrow(data)
    L <- dim(probs)[3]    
    skillclassprobsM <- matrix( skillclassprobs , nrow=N , ncol=L , byrow=TRUE )
    I <- ncol(data)
    # algorithm appr = 1 
	rho <- rep(.6 , N)
    appr.type <- 1
	if (progress){
		cat("***************************************\n")
		cat("Appropriateness Type " , appr.type )
					}
    res1 <- .calc.personfit.appr.algorithm( probs , data , data.resp , N , L , I , 
            appr.type  , rho , skillclassprobsM , eps=eps , maxiter=maxiter , conv=conv ,
            max.increment=max.increment , h=h , progress=progress )
    
    # algorithm appr = 0
    appr.type <- 0
	if (progress){
		cat("***************************************\n")	
		cat("Appropriateness Type " , appr.type )
					}	
	rho <- rep(.5,N)
    res0 <- .calc.personfit.appr.algorithm( probs , data , data.resp , N , L , I , 
            appr.type  , rho , skillclassprobsM , eps=eps , maxiter=maxiter , conv=conv ,
            max.increment=max.increment , h=h , progress=progress )
    # summaries
	dfr <- data.frame( "appr.type" = c(1,0) , "M.rho" = c( mean(res1$rho) , mean(res0$rho)) ,
					"SD.rho" = c( stats::sd( res1$rho) , stats::sd(res0$rho) )
						)
	dfr$median.SE.rho <- c( stats::median( res1$se.rho) , stats::median(res0$se.rho) )					
	dfr$prop.sign.T2 <- c( mean(res1$p<.05) , mean(res0$p <.05) )
	rownames(dfr) <- c("Spuriously High Scorers" ,
					    "Spuriously Low Scorers" )
    res2 <- list( "summary" = dfr , 
		"personfit.appr.type1"=res1 , "personfit.appr.type0" = res0 )
	class(res2) <- "personfit.appropriateness"	
    return(res2)
            }
####################################################################################
# S3 methods
# summary
summary.personfit.appropriateness <- function( object , digits=3 , ... ){
    print( round( object$summary , digits=3) )
                }
#***********************************
# plot method
plot.personfit.appropriateness <- function( x , cexpch=.65 , ... ){
	graphics::par(mfrow=c(2,2)) 
	# type=1	
	x1 <- x$personfit.appr.type1     
	N1 <- nrow(x1)	
	graphics::hist( x1$rho , main = "Appropriateness Type 1",
		freq=TRUE, breaks=seq(0,1 , length=20) , xlab=expression(rho) ,
		ylim=c(0,N1) )
	graphics::plot( c(0,1) , c(0,.5) , type="n" , xlab=expression(rho) , ylab= expression( p( T[2] ) ) ,
				main="Spuriously High Scorer"  )
	x1a <- x1[ x1$p >= .05 , ]
	graphics::points( x1a$rho , x1a$p , pch=1 , cex = cexpch )
	x1a <- x1[ x1$p < .05 , ]
	graphics::points( x1a$rho , x1a$p , pch=17 , cex = cexpch , col=2)
	# type=0
	x1 <- x$personfit.appr.type0     
	graphics::hist( x1$rho , main = "Appropriateness Type 0",
		freq=TRUE, breaks=seq(0,1 , length=20) , xlab=expression(rho) ,
		ylim=c(0,N1) )

	graphics::plot( c(0,1) , c(0,.5) , type="n" , xlab=expression(rho) , ylab= expression( p( T[2] ) ) ,
				main="Spuriously Low Scorer")
	x1a <- x1[ x1$p >= .05 , ]
	graphics::points( x1a$rho , x1a$p , pch=1 , cex = cexpch )
	x1a <- x1[ x1$p < .05 , ]
	graphics::points( x1a$rho , x1a$p , pch=17 , cex = cexpch , col=2)

	graphics::par( mfrow=c(1,1))
		}
####################################################################################
		
#########################################################################################
# algorithm calculation appropriateness statistics
.calc.personfit.appr.algorithm <- function( probs , data , data.resp , N , L , I , 
           appr.type  , rho , skillclassprobsM , eps=1E-15 , maxiter=30 , conv=1E-8 ,
           max.increment=.1 , h=.0001 , progress=TRUE ){
    rho <- rep(.4 , N )
    abs.incr <- 1
    iter <- 0
    while( ( abs.incr > conv) & ( iter < maxiter ) ){
        rho0 <- rho
        rho[ rho > 1 - h ] <- 1 - 2*h
        rho[ rho < 2*h ] <- 2*h
        ll0 <- .calc.ll.personfit.appropriateness( probs , data , data.resp , N , L , I , 
                appr.type=appr.type , rho=rho , skillclassprobsM , eps = eps )
        rho1 <- rho + h
        ll1 <- .calc.ll.personfit.appropriateness( probs , data , data.resp , N , L , I , 
                appr.type=appr.type , rho=rho1 , skillclassprobsM , eps = eps )
        rho2 <- rho - h
        ll2 <- .calc.ll.personfit.appropriateness( probs , data , data.resp , N , L , I , 
                appr.type=appr.type , rho=rho2 , skillclassprobsM , eps = eps )
        # first derivative
        deriv1 <- ( ll1 - ll0 ) / h
        # second derivative
        deriv2 <- ( ll1 - 2*ll0 + ll2 ) / h^2
        # update rho        
        increment <- deriv1 / abs( deriv2 )
        increment <- ifelse( abs(increment) > max.increment , max.increment*sign(increment) , increment )           
        rho <- rho + increment
        abs.incr <- max( abs( rho - rho0) )
		
		max.increment <- abs.incr * .9^(iter-1 )		
        iter <- iter+1
        if (progress){
          cat("\nIteration" , iter , "| Max. rho parameter change =" , round( abs.incr , 8 ) ) ;
          utils::flush.console()
                    }
                   }
		cat("\n")
        #******************* end iterations

        rho[ rho > 1 - h ] <- 1 - 2*h
        rho[ rho < 2*h ] <- 2*h
        
		# log-likelihood evaluated at rho=0		
		ll.rho0 <- .calc.ll.personfit.appropriateness( probs , data , data.resp , N , L , I , 
                  appr.type=appr.type , rho=0*rho+h/2 , skillclassprobsM , eps = eps )                    

		res <- data.frame( "rho" = rho , "se.rho"=sqrt( 1 / abs( deriv2 ) ) , 
				"ll.rho" = ll0 , "ll.0" = ll.rho0 , "T2" = 2*(ll0-ll.rho0)  )  
		res[ res$T2 < 0 , "T2" ] <- 0				
		res$p <- ( 1 - stats::pchisq( abs( res$T2 ) , df=1) ) / 2				 
		return(res)
        }
#########################################################################################

#########################################################################
# calculation of individual likelihood
.calc.ll.personfit.appropriateness <- function( probs , data , data.resp , N , L , I , 
        appr.type , rho , skillclassprobsM , eps = 1E-10 ){
    ll <- matrix( 1 , nrow=N , ncol=L)
    for (ii in 1:I){
        # ii <- 1
        prob.ii <- probs[ ii , 2 , ]
        prob.ii[ prob.ii < eps ] <- eps    
        prob.iiM <- matrix( prob.ii , nrow=N , ncol=L , byrow=TRUE )    
        prob.iiM <- (1-rho)*prob.iiM + rho * appr.type    
        ll.ii <- data[,ii] * prob.iiM + (1-data[,ii]) * (1-prob.iiM)
        ll.ii <- data.resp[,ii] * ll.ii + ( 1 - data.resp[,ii] )
        ll <- ll * ll.ii
                }			
    # calculate individual total log-likelihood
    ll0 <- log( rowSums( ll*skillclassprobsM ) )
    return(ll0)
            }    
#########################################################################