intToBinary<-function(number) {
#-------------------------------
# pass integer to binary
#-------------------------------

	if(number == 0)
		return (0)
	else
	{
		if(number == 1)
			return (1)
		else
		{
			aux = number%%2	
			number = number %/% 2
			return(c(intToBinary(number),aux))
		}
	}
}
