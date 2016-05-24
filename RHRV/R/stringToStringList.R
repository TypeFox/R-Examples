stringToStringList <- function(string) {
#-------------------------------
# convert one string in a list with the single strings 
# each string is 16 characters long, and trailing whitespaces are removed
# this function is used in LoadBeatEDFPlus
#-------------------------------
 	
	limit = nchar(string)
	num=limit/16
	final=c()
	for (i in 1:num) {
		start=(i-1)*16+1
		# cat("String no.",i,"\n")
		substring = substr(string,start,start+15)
		# cat("  ..",substring,"..\n",sep="")
		substring2 = gsub(" *$","", substring, perl=T)
		# cat("  ..",substring2,"..\n",sep="")
		final=c(final,substring2)
	}

	return(final)
}
