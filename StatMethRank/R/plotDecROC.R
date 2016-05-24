#' ROC Curves Plot of Decision Tree
#'
#' This function will draw ROC curves based on the tree information list.
#'
#' @param mytree a list returned by the decTreexxx function contains the information of the decision tree.
#' @export
#' @author Li Qinglong <liqinglong0830@@163.com>
#' @seealso plotDecTree
plotDecROC <- function(mytree) 
{
	plot(0:1, 0:1, type = 'l', xlab = "1-Specificity", ylab = "Sensitivity", main = "ROC Curve")
	for (iter in 1:length(mytree[[1]]))
	{
	    lines(mytree[[1]][[iter]], lty = iter)
	} 
	getAUC <- function(matSpeSen)
	{
		nRow = nrow(matSpeSen)
	    return(
	      round(sum(diff(matSpeSen[, 1]) * (matSpeSen[1:(nRow-1), 2] + matSpeSen[2:nRow, 2])) / 2, 3)
	      )
	}
	AUC = lapply(mytree[[1]], getAUC)
	legend.label = paste(names(mytree[[1]]), "(", AUC, ")", sep = "")
	legend("bottomright", legend = legend.label, lty = 1:length(mytree[[1]]), cex = 0.75)
}
