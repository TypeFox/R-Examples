#' Add Leading Zeroes to a Numeric Vector
#' 
#' Appends leading zeroes to a vector of numbers based on a string length or a
#' maximum number.
#' 
#' 
#' @param number A numeric vector.
#' @param number_length The length of the output string.
#' @param max_number A number to base the length of the output string on.
#' @return A character vector.
#' @author Jonathan A. Greenberg
#' @keywords format
#' @examples
#' 
#' x=c(1:10)
#' add_leading_zeroes(x,number_length=4)
#' add_leading_zeroes(x,max_number=10000)
#' 
#' @export

add_leading_zeroes=function(number,number_length,max_number)
{
	if(missing(max_number))
	{
		max_number=max(number)	
	}
	if(missing(number_length))
	{
		number_length=floor(log10(max_number))+1
	}
	
	
	fmt=paste("%0",number_length,"d",sep="")
	return(sprintf(fmt,number))
}
