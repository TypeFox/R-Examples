#' @title Defining labels of the binary partions.  
#' @details
#' Internal function. \code{bin.levels} is called by \code{splitopt.pls}.
#' @param variable segmentation variables (factor). 
#' @param spl matrix containing the posible combination of levels of factor.
#' @param \dots Further arguments passed on to \code{\link{bin.levels}}. 
#' @return list of matrices containing for each dependend latent variable the dependent LV and the explicative latent variables.
#' @keywords internal
#' @export

bin.levels	<-	function(variable,spl,...)
{
	new.lab = list()
	for (i in 1:2)
	{
		mod = levels(variable)[spl==i]
		for (i in 1:length(mod))
		{
			v.temp = paste(mod,collapse = '/')
		}
		new.lab[[length(new.lab)+1]] = v.temp
		}
	unlist(new.lab)
}
