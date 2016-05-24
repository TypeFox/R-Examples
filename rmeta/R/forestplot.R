

 
meta.colors<-function(all.elements,box="black",lines="gray",summary="black",zero="lightgray",
                      mirror="lightblue",text="black", axes="black",background=NA){

    if (missing(all.elements)){
        return(list(box=box, lines=lines, summary=summary,
             zero=zero, mirror=mirror, text=text,
             axes=axes, background=background))
    }
    
    if (is.null(all.elements)) 
        all.elements<-par("fg")
    
    return(list(box=all.elements, lines=all.elements,
             summary=all.elements, zero=all.elements,
             mirror=all.elements, text=all.elements,
             axes=all.elements, background=NA))

}



metaplot <- function( mn, se, nn=NULL, labels=NULL, conf.level = .95,
		      xlab = "Odds ratio", ylab = "Study Reference",
		       xlim = NULL, summn = NULL,
		      sumse = NULL, sumnn = NULL, 
		      summlabel = "Summary", logeffect = FALSE,
		      lwd = 2, boxsize = 1, 
		      zero = as.numeric(logeffect),
                      colors=meta.colors(), xaxt="s", logticks=TRUE,
		      ... ) {
    nth<-function(x,i){
        x[ (i-1) %% length(x) +1]
    }
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    ok <- is.finite( mn + se )
    if ( is.null( xlim ) ) 
        xlim <- c( min( mn[ok] - ci.value * se[ok], na.rm = TRUE ),
      		     max( mn[ok] + ci.value * se[ok], na.rm = TRUE ) )
    ##par( pty="s" )
    n <- length( mn )
    if ( logeffect ) {
        xlog <- "x"
        nxlim <- exp( xlim )
    }
    else {
        xlog <- ""
        nxlim <- xlim
    }

    leftedge<-nxlim[1]
    
    if ( !is.null( labels ) ) {
        if ( logeffect )  
            nxlim[1] <- nxlim[1] / sqrt( nxlim[2] / nxlim[1] )
        else
          nxlim[1] <- nxlim[1] - 0.5 * ( nxlim[2] - nxlim[1] )

        labels<-as.character(labels)
        
    }
    par( xaxt = "n",yaxt = "n", bg=colors$background )
    plot( nxlim,c( 1,-n-2-3 * !is.null( summn ) ),
          type = "n", bty = "n", xaxt = "n", yaxt = "n",
          log = xlog, xlab=xlab, ylab=ylab,..., col.lab=colors$axes )

    par( xaxt = "s" )
    if (xaxt=="s"){
        if (logeffect) {
            if (logticks){
                ats<-round( 10 ^ pretty( log( exp( xlim ),10), 8,min.n=6  ), 2 )
                ats<-ats[ats> exp(xlim[1]) & ats< 10^(par("usr")[2])]
                axis( 1, at = ats, col= colors$axes, col.axis= colors$axes)
            } else {
                ats<-pretty(exp(xlim),8, min.n=6)
                ats<-ats[ats> exp(xlim[1]) & ats <10^(par("usr")[2])]
                axis( 1, at=ats, col= colors$axes, col.axis= colors$axes)
            }
        }  else {
            ats<-pretty(xlim, 6)
            ##ats<-ats[ats> xlim[1] & ats <xlim[2]]
            axis( 1, at=ats, col= colors$axes, col.axis= colors$axes)
        }
    }
    
    if ( !is.null( zero )&& zero>leftedge )
        abline( v = zero, lty = 2, lwd = 2 ,col=colors$zero)

    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    lower <- mn - ci.value * se
    upper <- mn + ci.value * se
    if ( logeffect ){
        lower <- exp( lower )
        upper <- exp( upper )
    }
    for ( i in 1:n ){
        if ( is.na( lower[i]+upper[i] ) ) 
            next
        lines( c( lower[i], upper[i] ), c( -i, -i ), lwd = lwd, col=nth(colors$lines,i),... )
    }

    if ( !is.null( labels ) )
        text( rep( nxlim[1], n ), -( 1:n ), labels,..., col=rep(colors$text,length.out=n),adj=0 )

    if ( is.null( nn ) ) 
        nn <- se ^ -2
    yscale <- 0.3 * boxsize / max( sqrt( nn ), na.rm = TRUE )

    if ( logeffect ) { 
        scale <- ( nxlim[2] / nxlim[1] ) ^ ( yscale / ( 4 + n ) )
        xl <- exp( mn ) * ( scale ^ -sqrt( nn ) )
        xr <- exp( mn ) * ( scale ^ sqrt( nn ) )
    }
    else {
        scale <- yscale * ( nxlim[2] - nxlim[1] ) / ( 4 + n )
        xl <- mn - scale * sqrt( nn )
        xr <- mn + scale * sqrt( nn )
    }
    yb <- ( 1:n ) - yscale * sqrt( nn )
    yt <- ( 1:n ) + yscale * sqrt( nn )
    for ( i in 1:n ) {
        if ( !is.finite( mn[i] ) ) 
            next  
        rect( xl[i], -yb[i], xr[i], -yt[i], col = nth(colors$box,i),border=nth(colors$box,i))
    }
    if ( !is.null( summn ) ) {
        if ( logeffect ) {
            x0 <- exp( summn )
            xl <- exp( summn - ci.value * sumse )
            xr <- exp( summn + ci.value * sumse )
        }
        else{
            x0 <- summn
            xl <- summn - ci.value * sumse
            xr <- summn + ci.value * sumse
        }
        y0 <- n + 3
        yb <- n + 3 - sqrt( sumnn ) * yscale
        yt <- n + 3 + sqrt( sumnn ) * yscale
        polygon( c( xl, x0, xr, x0 ), -c( y0, yt, y0, yb ),
    	         col = colors$summary, border = colors$summary )
        text( nxlim[1], -y0, labels = summlabel, adj = 0,col=colors$text )
    }
}



funnelplot.default<-function(x,se,size=1/se,summ=NULL,xlab="Effect",ylab="Size",
		colors=meta.colors(),
		conf.level=0.95,plot.conf=FALSE,zero=NULL,mirror=FALSE,...)
{
   finite<-function(x) x[is.finite(x)]
    
   if (mirror && is.null(summ))
	stop("Can't do a mirror plot without a summary value")

   if (plot.conf){
	ci<--qnorm((1-conf.level)/2)
	xlim<-range(finite(c(zero,x-ci*se,x+ci*se)))
	if (mirror)
	  xlim<-range(finite(c(xlim,2*summ-x-ci*se,2*summ-x+ci*se)))
   }else{
      xlim<-range(finite(c(zero,x)))
      if (mirror)
	xlim<-range(finite(c(xlim,2*summ-x,2*summ+x)))
   }
   plot(x,size,ylim=c(0,max(size)*1.1),xlim=xlim,xlab=xlab,
	ylab=ylab,col=if(is.null(colors$points)) par("fg") else colors$points)

   if (plot.conf)
       segments(x-ci*se,size,x+ci*se,size,col=if(is.null(colors$conf)) par("fg") else colors$conf,lwd=2)
   if (!is.null(summ))
	   abline(v=summ,col=if(is.null(colors$summary)) par("fg") else colors$summary,lty=2,lwd=2)

   if(!is.null(zero))
	abline(v=zero,col=if(is.null(colors$zero)) par("fg") else colors$zero,lwd=2)

   if(mirror){
	points(2*summ-x,size,col=if(is.null(colors$mirror)) par("fg") else colors$mirror)
       if (plot.conf)
	segments(2*summ-x-ci*se,size,2*summ-x+ci*se,size,col=if(is.null(colors$mirror)) par("fg") else colors$mirror,lwd=2)
  }
}




