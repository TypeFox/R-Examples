"tree.demo" <-
function(x,y){
	old.opt <- options(locatorBell = FALSE)
        on.exit( options(old.opt) )
	cuts <- range(x)
	
	repeat {
		
		cut2 <- numeric(0)
		repeat {
			plot(x,y,xlab=deparse(substitute(x)),
                             ylab=deparse(substitute(y)))
			abline( v=cuts, col='blue' )
			abline( v=cut2, col='red' )
			cuts3 <- sort( c(cuts,cut2) )
			cats <- cut( x, cuts3, include.lowest=T)
			means <- tapply(y, cats, mean )
			index <- tapply(y, cats )
			segments(cuts3[-length(cuts3)], means, cuts3[-1], means, col='green' )
			resid <- y-means[index]
			ss <- round(resid %*% resid)
			title( paste( "Residual sum of squares =", ss ) )
			tempx <- locator(1)$x
			if (length(tempx) < 1) break
			cut2 <- tempx
		}
		if(length(cut2) < 1) break
		cuts <- sort( c(cuts,cut2) )
	}
	
	
}

