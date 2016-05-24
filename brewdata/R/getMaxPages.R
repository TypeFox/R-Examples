getMaxPages <-
function ( url ) {
	# Finds the number of pages of data to ensure getGradCafeData exits.
	webpage = getURL( url )
	webpage = readLines( tc <- textConnection(webpage) ); close(tc)
	
	#defines the line of data needed to find the number of pages
	pgs= "\t\t\t\tShowing <strong>[0-9]* results</strong> over [0-9]* pages"
	max_index = grep( pgs, webpage )
	
	#return the number of data pages
	as.numeric( strsplit( webpage[max_index], " " )[[1]][5] )
}