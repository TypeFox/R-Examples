#' @title Defining labels of a path coefficient.  
#' @details
#' Internal function. \code{element} is called by \code{F.data.pls}.
#' @param x matrix containing information of casual relationship of latent variables.
#' @param \dots Further arguments passed on to \code{\link{element}}. 
#' @return vector path coefficients labels.
#' @keywords internal
#' @export

element	<-	function(x,...)
{
	y = mat.or.vec(1,sum(ifelse(x==0,0,1)))
	k = 1
	for (i in 1:dim(x)[1]) 
	{
   		for (j in 1:dim(x)[1])
   		{
      		if(x[i,j] != 0)
      		{
				v	= paste(colnames(x)[j],rownames(x)[i],sep = "->")
           		y[k]	= v
           		k	= k+1
      		}
   		}
	}
	as.vector(y)
}

