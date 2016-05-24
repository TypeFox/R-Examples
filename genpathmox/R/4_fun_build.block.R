#' @title Defining the linear relations between the latent variables. 
#' @details
#' Internal function. \code{build.block} is called by \code{splitopt.pls}.
#' @param inner matrix containing information of the latent causal relationships.
#' @param latent matrix containing the estimated latent variables.
#' @param \dots Further arguments passed on to \code{\link{build.block}}.
#' @return list of matrices containing for each dependend latent variable the dependent LV and the explicative latent variables.
#' @keywords internal
#' @export

build.block	<-	function(inner,latent,...)
{
	constructs			= list() 
	constructs.label	= list()
	
	for (i in 1:ncol(inner))
	{
		y = NULL
		y = as.matrix(latent[,i])
		colnames(y) = rownames(inner)[i]
		x = NULL
		for (j in 1:ncol(inner))
		{
			xj = NULL
			if (inner[i,j] == 0) next
			xj = as.matrix(latent[,j])
			colnames(xj) = rownames(inner)[j]
			x = cbind(x,xj)						
		}
		constructs[[length(constructs)+1]] = cbind(y,rep(1,nrow(y)),x)
		constructs.label[[length(constructs.label)+1]] = c("int",colnames(x))		
	}

	new.constructs		= list()
	new.constructs.label	= list()	
	k=1
    while (k <= length(constructs))
    {
    	if (ncol(constructs[[k]]) > 2)
    	{
    		 new.constructs[[length(new.constructs)+1]] = constructs[[k]]
    		 new.constructs.label[[length(new.constructs.label)+1]] = constructs.label[[k]]
    	}
    	k=k+1	
    }

	label.block = NULL
 
	for (i in 1:length(new.constructs)) label.block[i] = colnames(new.constructs[[i]])[1]
	for (i in 1:length(new.constructs.label)) names(new.constructs.label)[i] = label.block[i]

	resp = NULL
	pred = list()
	
	for (i in 1:length(new.constructs))
	{
		yi		= NULL
		yi		= as.matrix(new.constructs[[i]][,1])	
		resp	= rbind(resp,yi)
		pred[[length(pred)+1]]	= as.matrix(new.constructs[[i]][,-1])
	}
		
	x.block = blockdiag(pred)	
	colnames(x.block) = c(unlist(new.constructs.label))

	list(x.block=x.block,
		resp=resp,
		constructs.label=new.constructs.label)
}

