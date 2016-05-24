#' Compute 1 to k-step transition proportions
#' 
#' Computes 1 to k-step forward transition proportions in each state with a single transition matrix or a 3-d array of transition matrices.
#' 
#' @param x either a transition matrix, list of transition matrices or 3-d array ( a set of transition matrices)
#' @param k if x is a transition matrix, this is number of steps 1 to k
#' @param labels labels for states except for last which is always dead and is added at end
#' @author Jeff Laake
#' @export 
#' @importFrom expm %^%
#' @keywords utility
omega=function(x,k=NULL,labels=NULL)
{
	if(is.list(x))
	{
		y=x
		x=array(NA,dim=c(length(y),nrow(y[[1]]),ncol(y[[1]])))
		for(i in 1:length(y))
			x[i,,]=y[[i]]
	}	
	array_dim=dim(x)
	if(length(array_dim)==2)
	{
		if(is.null(k))stop("k - number of transitions must be specified")
		mat=t(sapply(1:k,function(k,x) (x%^%k)[1,] ,x=x))
	} else
	{
		xmat=x[1,,]
		k=array_dim[1]
		mat=matrix(NA,nrow=k,ncol=nrow(xmat))
		mat[1,]=xmat[1,]
		for(j in 2:k)
		{
		  xmat=xmat%*%x[j,,]
		  mat[j,]=xmat[1,]
		}
	}
	if(!is.null(labels))colnames(mat)=c(labels,"Dead")
	rownames(mat)=1:k
	return(mat)
}
