#' Set fixed real parameter values in ddl
#' 
#' Merges fixed real parameter values in parameter specification with those specified with fix in the 
#' design data list and returns the fixed values in the design data list.
#' 
#' @param ddl design data list
#' @param parameters parameter specififcation list
#' @return design data list with fixed values
#' @author Jeff Laake
#' @keywords utility
set.fixed=function(ddl,parameters)
{
   for(parx in names(parameters))
   {
	   if(!is.null(parameters[[parx]]$fixed))
		   ddl[[parx]]$fix[ddl[[parx]]$id%in%parameters[[parx]]$fixed[,1]&as.numeric(ddl[[parx]]$time)%in%parameters[[parx]]$fixed[,2]]=parameters[[parx]]$fixed[,3]
   }
   return(ddl)
}
#' Create parameters with fixed matrix
#' 
#' Creates fixed matrix in parameters from ddl$fix values
#' 
#' @param ddl design data list
#' @param parameters parameter specififcation list
#' @return parameters with fixed matrix set
#' @author Jeff Laake
#' @keywords utility
create.fixed.matrix=function(ddl,parameters)
{
	for(parx in names(parameters))
	{
		if(!is.null(ddl[[parx]]$fix))
		{
			select=!is.na(ddl[[parx]]$fix)
			parameters[[parx]]$fixed=cbind(as.numeric(ddl[[parx]]$id[select]),as.numeric(ddl[[parx]]$time[select])+1,ddl[[parx]]$fix[select])
		} else
		    parameters[[parx]]$fixed=matrix(c(-1,-1,0),nrow=1,ncol=3)
	}	
	return(parameters)
}
	