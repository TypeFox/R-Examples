binaryToInt<-function(number) {
#-------------------------------
# pass binary to integer
#-------------------------------

	count = length(number)
	result = 0
	for(i in 1:count)
	{
		result = result + (number[i]*2**(count-i))
	}
	return(result)
}
