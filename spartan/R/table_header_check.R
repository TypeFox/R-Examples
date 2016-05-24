table_header_check <- function(VECTOR)
{
	# Simple function to replace any occurences of a hyphen or space in a parameter or measure with a dot
	for(ENTRY in 1:length(VECTOR))
	{
		VECTOR[ENTRY]<-gsub("-",".",VECTOR[ENTRY])
		VECTOR[ENTRY]<-gsub(" ",".",VECTOR[ENTRY])
	}
	return(VECTOR)
}
