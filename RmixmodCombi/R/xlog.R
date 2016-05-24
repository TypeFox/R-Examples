xlog <-
function(x) 
	{
		xlog1d <- function (xi) if (xi == 0) 0 else (xi*log(xi))
		
		if (is.null(dim(x)))
			{
				return(sapply(x,xlog1d))
			}
		else
			{
				return(matrix(sapply(x,xlog1d),dim(x)))
			}
	}
