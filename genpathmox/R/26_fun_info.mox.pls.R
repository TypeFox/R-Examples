#' @title Giving general information about the pathmox algorithm
#' @details
#' Internal function. \code{info.mox.pls} is called by \code{pls.pathmox}.
#' @param signif stop condition 1: significance of the p-value
#' @param size stop condition 2: minimum number of individuals  in a node
#' @param deep stop condition 3: maximum tree depth level
#' @param y: set of segmentation variables
#' @param \dots Further arguments passed on to \code{\link{info.mox.pls}}. 
#' @keywords internal

#' @export

info.mox.pls	<-	function(signif,size,deep,y,...)
{
	cat("\n")
	cat("PLS-PM SEGMENTATION TREE","\n")
	cat("\n")
	cat("---------------------------------------------")
	cat("\n")
	cat("Info Parameters Algorithm","\n")
	info.value = rbind(signif,size,deep)
	dimnames(info.value) = NULL
	info.name = c("Threshold signif","Node size limit(%)","Tree depth level")
	info.tree = data.frame(info.name,info.value)
	names(info.tree) = c("parameters Algorithm", "value")
	print(info.tree)
	cat("\n")
	cat("---------------------------------------------")
	cat("\n")
	cat("Info Segmentation Variables","\n")
	type.y = rep(0, ncol(y))
	treat.y = rep("binary", ncol(y))
	for (i in 1:length(type.y))
	{
        type.y[i] = ifelse(is.ordered(y[, i]), "ord","nom")
       
        if (nlevels(y[, i]) > 2) 
            if (is.ordered(y[, i])) 
                treat.y[i] = "ordinal"
            else treat.y[i] = "nominal"
    }
    df.y = data.frame(Nlevels = unlist(lapply(y, nlevels)),Ordered = unlist(lapply(y, is.ordered)), Treatment = treat.y)
   	print(df.y)
}
