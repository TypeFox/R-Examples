#' Merge two parameterDef objects. 'others' of two parameterDef objects must be the same. 'scalars' of two parameterDef objects must have same name and the value must not overlap.
#' @name merge.parameterDef
#' @aliases merge.parameterDef
#' @title Merge two parameterDef objects
#' @method merge parameterDef
#' @param x A parameterDef to merge with 
#' @param y A parameterDef to merge with
#' @param ... unused
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export 
#' @seealso \code{\link{createParDef}}, \code{\link{createParDef}}
merge.parameterDef<-function(x,y,...){
	if (class(y)!='parameterDef')
		stop('y must be an parameterDef object')
		
	if ( digest(x$banker) != digest(y$banker) )
		stop(' \'banker\' of two parameterDef must be the same.')
		
	if (!identical(names(x$selection),names(y$selection)))
		stop('selection of x and y must have same set of names')
		
	i=j=name_i=NULL
	new_selection<-
		foreach(i = x$selection, j =y$selection,name_i=names(x$selection),.combine=c ) %do% {
			out<-list(sort(unique(c(i,j))))
			names(out)<-name_i
			out
		}
		
	createParDef(selection=new_selection, banker=x$banker)
}
