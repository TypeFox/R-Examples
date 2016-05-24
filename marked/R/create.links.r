#' Creates a 0/1 vector for real parameters with sin link
#' 
#' For each row in a given design matrix it assigns a value 1, if the columns
#' used in the design matrix are only used as an identity matrix (only one 1
#' and remaining columns all 0.
#' 
#' 
#' @param dm design matrix
#' @return A vector of length=nrow(dm) with value of 0 except for rows that can
#' accommodate a sin link (identity design matrix).  This function is not
#' currently used because it has not been thoroughly tested.
#' @author Jeff Laake
create.links=function(dm)
{
	links=rep(0,nrow(dm))
	possible.sin=apply(dm,1,function(x){
				w=which(x==1)
				ifelse(length(w)==1,w,NA)
			})								 
	valid.columns=(1:ncol(dm))[sapply(1:ncol(dm),function(x) all(dm[dm[,x]==1,-x]==0))]
	if(length(valid.columns)>0)
		for (i in valid.columns)
			links[possible.sin==i]=1
	return(links)
}
