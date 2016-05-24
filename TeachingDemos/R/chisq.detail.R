"chisq.detail" <-
function(tab){
	d <- dim(tab)
	ct <- colSums(tab)
	rt <- rowSums(tab)
	tt <- sum(rt)
	ev <- ( rt %o% ct )/ tt
	ch2 <- (tab - ev)^2 / ev
	out1 <- matrix( "", ncol = d[2] + 1, nrow = d[1]*3 + 1)
	if( is.null( dimnames(tab) ) ){
		dimnames(out1) <- list( c( rep("",d[1]*3), "Total"), c( rep("", d[2]), "Total") )
	} else {
		dimnames(out1) <- list( c( rbind(dimnames(tab)[[1]],"",""), "Total" ), 
									c( dimnames(tab)[[2]], "Total" ) ) 
	}	
	out1[ 3*(1:d[1])-2, 1:d[2] ] <- paste(tab,"   ", sep="")
	out1[ 3*(1:d[1])-1, 1:d[2] ] <- format(round(ev,2), nsmall=2)
	out1[ 3*(1:d[1])-2, d[2]+1] <- rt
	out1[ 3*d[1]+1, 1:d[2] ] <- paste(ct,"   ",sep="")
	out1[ 3*d[1]+1, d[2]+1 ] <- tt
	
	cat("\n\nobserved\nexpected\n\n")
	print(out1, quote=FALSE, right=TRUE)
	
	out2 <- matrix("", nrow=d[1], ncol= 2*d[2]+1)
	if( is.null( dimnames(tab) ) ){
		dimnames(out2) <- list( rep("",d[1]), rep("", d[2]*2+1) )
	} else {
		dimnames(out2) <- list( dimnames(tab)[[1]], 
									c( rbind(dimnames(tab)[[2]],""), "" ) ) 
	}
	out2[ 1:d[1], 2*(1:d[2])-1 ] <- format(round(ch2,2), nsmall=2)
	out2[ 1:d[1], 2*(1:d[2])   ] <- "+"
	out2[   d[1], 2*   d[2]    ] <- "="
	out2[   d[1], 2*   d[2] +1 ] <- round( sum( ch2 ),2 )
	
	cat("\n\nCell Contributions\n")
	print(out2, quote=FALSE, right=TRUE) 
	cat("\ndf =", (d[1]-1)*(d[2]-1), " P-value =", 
	    round( 1 - pchisq(sum(ch2), (d[1]-1)*(d[2]-1)), 3),
	    "\n\n" )
	
	invisible(list( obs = tab, expected = ev, chi.table = ch2, chi2 = sum(ch2) ) )
}

