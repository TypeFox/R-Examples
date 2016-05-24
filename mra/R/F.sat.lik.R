"F.sat.lik" <-
function( ch ){

#	Currently, only works if no removals
if( any(ch == 2) ) stop("Current routine cannot compute saturated likelihood in presence of censored animals.\ni.e., no 2's allowed.")
	
#	Calculate 'n' the number of animals first released on each occasion: 
temp <- col(ch) * (ch >= 1)
temp <- apply( temp, 1, function(x) { min( x[x>0] ) } )
n <- sapply( 1:ncol(ch), function( x, temp ) { sum( temp == x, na.rm = TRUE ) }, temp ) 

#	Convert to characters and use table() to count uniques
ch.char <- NULL
for( i in 1:ncol(ch)){
	ch.char <- paste( ch.char, ch[,i], sep="" )
}
njj <- table( ch.char )

#	Need vector indicating the occasion of first release for each of the unique
#	capture histories in 'x.unique'.  
x.unique <- names( njj )
ind <- sapply( x.unique, function(x){ 
		tmp <- as.numeric(substring( x, 1:nchar(x), 1:nchar(x) )) # convert to integer vector
		tmp <- (1:length(tmp)) * tmp
		tmp <- min( tmp[ tmp>0 ] )
		tmp} )


sat.lik <- -2 * sum( njj * log( njj / n[ind] ) )

return(sat.lik)

}

