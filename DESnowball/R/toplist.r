#' select the top list of genes
#'
#' Report the top list based on p values.
#' @param fs an object output from function \code{select.features}
#' @return a data.frame with two columns \code{RD} and \code{pvalue} (see references for details)
#' @export
#' @references
#' Xu, Y., Guo, X., Sun, J. and Zhao. Z. Snowball: resampling combined with distance-based regression to discover transcriptional consequences of driver mutation, manuscript.
toplist <- function(fs)
{
   ret <- fs$selectedList[,c('rd', 'pval')]
   names(ret) <- c("RD", "pvalue")
   invisible(ret)
}
