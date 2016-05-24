#' ':' operator which won't allow you to go "backwards"
#' 
#' 
`%>%` <-
	function (low, high) 
	{
		if (low <= high) {
			return(low:high)
		}
		else {
			return(seq_len(0))
		}
	}


#' Concat operator
#' 
#' 
`%&%` <-
	function (x, y) 
	{
		paste(x, y, sep = "")
	}
