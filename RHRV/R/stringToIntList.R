stringToIntList <- function(string) {
#-------------------------------
# convert one string in a list with the numbers in the string
# example: "111 222 333" -> c(111,222,333)
#-------------------------------
 
	limit = nchar(string)
	aux=""
	final=c()
	for(i in 1:limit)
	{
		if(substr(string,i,i) == " ")
		{
			if(nchar(aux) > 0)
			{
				final=c(final,as.integer(aux))
				aux=""
			}
		}
		else
		{
			aux = paste(aux,substr(string,i,i),sep="")
		}
	}
	return(final)
}
