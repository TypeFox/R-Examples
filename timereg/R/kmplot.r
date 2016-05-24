
kmplot<- function(x,loc=NULL,col=NULL,lty=NULL,...)
{ ## {{{ 
	### default location if loc not given 
	if (is.null(loc)) { 
		if (min(x$surv)>0.7) loc <- "bl" else loc <- "bl"
	}
	if (loc=="bl") loc <- "bottomleft"
	else if (loc=="br") loc <- "bottomright"
	else if (loc=="tr") loc <- "topright"
	else if (loc=="tl") loc <- "topleft"
	else loc <- "bottomleft"
	nn <- names(x$strata)
	ll <- length(nn)
	if (is.null(col)) cols <- seq(1,ll) else cols <- col
	if (is.null(lty)) ltys <- seq(1,ll)  else ltys <- lty
	plot(x,col=cols,lty=ltys,...)
	legend(loc,legend=names(x$strata),col=cols,lty=ltys)
} ## }}}

