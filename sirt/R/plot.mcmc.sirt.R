
######################################################
# plot results of objects of class mcmc.sirt
plot.mcmc.sirt <- function( x , layout=1 , conflevel=.90 , 
	round.summ=3 , lag.max= .1 , col.smooth="red" , lwd.smooth=2 , 
	col.ci="orange" , cex.summ=1 , ask=FALSE , ... ){
	
	object <- x	# rename x into object
	mcmcobj <- (object$mcmcobj)[[1]]
	lag.max <- round( nrow(mcmcobj) * lag.max )
	
	# layout type
	# layout=1 : standard output from coda package	
    if (layout==1){  plot(object$mcmcobj , ask=ask , ...) }         

	#***************
	# layout =2
	if (layout==2){
		mcmcobj <- (object$mcmcobj)[[1]]
		lag.max <- min( nrow(mcmcobj) , lag.max )
		# index vector
		a1 <- attr(mcmcobj,"mcpar")
		iterindex <- seq(a1[1] , a1[2] , a1[3] )
		smcmcobj <- object$summary.mcmcobj
		VV <- ncol(mcmcobj)		
		ci.quant <- - stats::qnorm( (1-conflevel)/2 )
		graphics::par( mfrow=c(2,2))
		for (vv in 1:VV){
#	vv <- 15
			x.vv <- as.vector( mcmcobj[,vv] )
			parm.vv <- colnames(mcmcobj)[vv]
			sparm.vv <- smcmcobj[ smcmcobj$parameter == parm.vv , ]
			#***
			# traceplot
			graphics::plot( iterindex , x.vv , type="l" ,  main= paste0( "Traceplot of " , parm.vv ) ,
				xlab="Iterations" , ylab="" , ... )
			x1 <- as.numeric( x.vv )
			xmin <- min(x1)
			xmax <- max(x1)
			# l1 <- loess( x1 ~ iterindex )$fitted 
			# include moving average here!!
			l1 <- .movingAverage(x1 , B = round( lag.max / 2 ) , fill=FALSE)
			graphics::lines( iterindex ,l1  , col= col.smooth , lwd= lwd.smooth )
			#***
			# density estimate
			graphics::plot( stats::density( x.vv ) , main= paste0( "Density of " , parm.vv ) )
			
			c1 <- stats::quantile( x1 , ( 1 - conflevel  ) / 2 )			
			c2 <- stats::quantile( x1 , 1 - ( 1 - conflevel  ) / 2 )			
#			lines( sparm.vv$Mean + c(-1,1)*ci.quant * sparm.vv$SD , c(0,0) , col=col.ci , lwd=3 )
			graphics::lines( c(c1,c2) , c(0,0) , col=col.ci , lwd=3 )
			graphics::points( sparm.vv$Mean , 0 , pch=17 , col= col.ci , cex=1.5)
			#***
			# autocorrelation function
			stats::acf( x.vv , lag.max=lag.max ,
				main= paste0( "Autocorrelation of " , parm.vv ) )
			#***
			# numerical summary
			graphics::plot( c(0,1) , c(0,1) , axes=FALSE , xlab="" , ylab="", 
					main= paste0( "Summary of " , parm.vv ) , type="n" , ...)
			x0 <- 0 ; y0 <- 0
			# heights.summ = c( .05 ,  .20 , .35 ,  .5 , .65 , .8 , .95)
			heights.summ = c( .05 ,  .15 , .25 ,  .35 , .45 , .55 , .65 , .75)
			graphics::text( x0 + .0015 , y0 + heights.summ[8] , "Posterior Mean =" , cex= cex.summ  , pos=4)
			graphics::text( x0 + .5 , y0 + heights.summ[8] , 
				paste0( format.numb( x = mean( x1 ) , digits = round.summ)  )  , pos=4 )
			hvv <- heights.summ[7]	
			graphics::text( x0 + .0015 , y0 + hvv , "Posterior Mode =" , cex= cex.summ  , pos=4)				
			graphics::text( x0 + .5 , y0 + hvv , 
				paste0( format.numb( x = sparm.vv$MAP , digits = round.summ)  )  , pos=4 )			
				
			graphics::text( x0 + .0015 , y0 + heights.summ[6] , "Posterior SD   =" , cex= cex.summ  , pos=4)
			graphics::text( x0 + .5 , y0 + heights.summ[6] , 
				paste0( format.numb( x = stats::sd( x1 ) , digits = round.summ)  )  , pos=4 )

			hvv <- heights.summ[5]
			graphics::text( x0 + .0015 , y0 + hvv , 
							paste( round(100*conflevel ) , "% Credibility Interval = " ,sep="") ,
							cex= cex.summ , pos=4 )

			hvv <- heights.summ[4]
				ci.lower <- format.numb( stats::quantile( x1 , ( 1 - conflevel  ) / 2 ) , digits = round.summ )
				ci.upper <- format.numb( stats::quantile( x1 , 1- ( 1 - conflevel  ) / 2 ) , digits = round.summ )            
			graphics::text( x0 + .25 , y0 + hvv , 
							paste( "[" , ci.lower ,    "," , ci.upper , "]" ,  sep="") ,
							cex= cex.summ  , pos=4)
			hvv <- heights.summ[3]
			graphics::text( x0 + .0015 , y0 + hvv , "Rhat =" , cex= cex.summ  , pos=4)
			graphics::text( x0 + .5 , y0 + hvv , 
				paste0( format.numb( x = sparm.vv$Rhat , digits = 2)  )  , pos=4 )
			hvv <- heights.summ[2]
			graphics::text( x0 + .0015 , y0 + hvv , "PercSERatio =" , cex= cex.summ  , pos=4)
			graphics::text( x0 + .5 , y0 + hvv , 
				paste0( format.numb( x = sparm.vv$PercSERatio , digits = 1)  )  , pos=4 )

			hvv <- heights.summ[1]
			graphics::text( x0 + .0015 , y0 + hvv , "Effective Sample Size =" , cex= cex.summ  , pos=4)
			graphics::text( x0 + .705 , y0 + hvv , 
				paste0( format.numb( x = sparm.vv$effSize , digits = 1)  )  , pos=4 )
				
				
			graphics::par(ask=ask)				
					}				
			graphics::par(mfrow=c(1,1))
			 }

                    }
#######################################################
# format numbers
format.numb <- function( x , digits ){ 
    a1 <- round( x , digits ) + 10^{-(digits +1 ) } 
    a1 <- substring( a1 , 1 , nchar(a1) - 1 )
    return(a1)
    }
#######################################################



# moving window average for time series
.movingAverage <- function(x, B, fill=TRUE){

    x1 <- cumsum(x)
    N <- length(x)
    y <- rep(NA,N)
    i <- seq(B+1 , N-B)
    xdiff <- x1[ -seq(1,B) ] - x1[ -seq(N-B+1,N) ] 
    xdiff <- xdiff[ - seq(1,B) ] 

    y[i]  <- ( x1[i] + xdiff - c(0,x1[ -seq(N-2*B,N) ]) ) / (2*B+1)

  # fill NAs at beginning and end of time series
  if(fill){
    j <- seq(0,B-1)
    ybeg <- sapply(j, function(z) sum( x[ seq(1,(2*z+1)) ]) / (2*z+1) )
    yend <- sapply(rev(j), function(z) sum( x[ seq(N-2*z,N) ] ) / (2*z+1) )
    y[j+1] <- ybeg
    y[ base::rev(N-j) ] <- yend
  }

  y
}	