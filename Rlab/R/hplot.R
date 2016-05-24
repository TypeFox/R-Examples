"hplot" <-

function(x, breaks="Sturges", freq=FALSE, nclass=NULL, col=8, ...)



{



	if ( is.null(nclass) )



		hist(x, breaks=breaks, freq, col=col, ...)



	else

	

		hist(x, breaks=seq(min(x),max(x),length=nclass+1), freq, col=col, ...)



}

