fitgraph.default <-
function( fitEst, fitLB, fitUB, itemLabels
	                 , mainTitle   = 'Fit Plot'
					 , pch    = 18
	                 , fitColours = c('gray70','gray60','gray50','gray40','gray0')
	                 , xlab = "Items"
	                 , cex = 1.25
	                 ,... )
{

	par(mar=c(4.1, 4.1, 3.1, 1.1))

	nI <- length(fitEst)
	item <- c(1:nI)
	#dev.new()
	plot( fitEst ~ item
		, type = "n"
		, axes = FALSE
		, ylab = "Fit"
		, xlab = ""
		, ylim = c( max( min( fitEst ) - .3, 0) , max( fitEst ) + .3)
		, xlim = c( min( item ), max( item )))
		

	title( main = mainTitle )
	
	xcorrs <-c( 0, rep( 1:( nI - 1), each = 2), nI)
	xcorrs <- xcorrs + .5
	yUBs = rep(fitUB,each = 2)
	yLBs = rep(fitLB,each = 2)
	
	
	polygon( c( min( xcorrs ), xcorrs, max(xcorrs)), c( 0, yUBs, 0 )
		   , col    = fitColours[1], border = NA) 
		   
		
	polygon( c( min( xcorrs ), xcorrs, max(xcorrs)), c( 0, yLBs, 0 )
		   , col    = "white", border = NA)
		   
		   
		
	 lines( c(head(yLBs, n = 1),yUBs) ~ c(head(xcorrs, n = 1),xcorrs)
	      , type = "l", col  = fitColours[3]
	      , lwd  = 2, lend = 2, ljoin = 2, lty  = "solid") 
    
	 lines( c(yLBs, tail(yUBs, n = 1))  ~ c(xcorrs,tail(xcorrs, n = 1))
	      , type = "l", col  = fitColours[3]
	      , lwd  = 2, lend = 2, ljoin = 2, lty  = "solid") 

	axis( 2, las = 1, tcl = -0.5
		, at = 1
		, cex.axis  = 1, font.axis = 2, col.axis  = fitColours[5])

	axis( 2, las = 1, tcl = -0.2
		, at = seq( from = 1.2, to = ( max( fitEst ) + .3 ), by = .2)
		, cex.axis  = .8, font.axis = 2, col.axis  = fitColours[5])

	axis( 2, las = 1, tcl = -0.2
		, at = seq( from = 0.8, to = ( min( fitEst ) - .3 ), by = -.2)
		, cex.axis  = .8, font.axis = 2, col.axis  = fitColours[5])

	box() 

	abline( h = as.list( seq( from = 1.0, to = ( max( fitEst ) + .3 ), by = .2))
	      , lwd = .5,col = fitColours[5])
	
	abline( h = as.list( seq( from = 0.8, to = ( max( min( fitEst ) - .3, 0) ), by = -.2))
	      , lwd = .5,col = fitColours[5])

	lines( fitEst ~ item
	     , type = "p"
	     , pch  = pch
	     , ylim = c( .5, 1.5)
	     , col  = fitColours[5]
	     , bg   = fitColours[4]
	     , xlab = 'Items'
	     , cex  = cex)
    
	text( item, fitEst, as.list( itemLabels), pos = 3, cex = .75)

	par(mar=c(5.1,4.1,4.1,2.1))
}
