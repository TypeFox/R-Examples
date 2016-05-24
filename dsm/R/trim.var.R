#' Trimmed variance
#'
#' Trim the variance estimates from the bootstrap. This is defined as the 
#' percentage defined as amount necessary to bring median and trimmed mean 
#' within 8% of each other these are defined as 'outliers'. 
#'
#' Code originally by Louise Burt.
#' 
#' @param untrimmed.bootstraps (usually the \code{$study.area.total} element
#'        of a returned \code{dsm} bootstrap object.
#' @param boxplot.coef the value of \code{coef} used to calculate the outliers
#'        see \code{\link{boxplot}}.
#'
#' @return trimmed variance
#'
#' @export
#' @importFrom stats var
#' @importFrom grDevices boxplot.stats
#' @author Louise Burt
trim.var <- function(untrimmed.bootstraps, boxplot.coef=1.5){
  outliers <- boxplot.stats(untrimmed.bootstraps, coef=boxplot.coef)$out
  bootstrap.abund <-untrimmed.bootstraps[!(untrimmed.bootstraps %in% outliers)]

  ret <- var(bootstrap.abund)

  attr(ret,"trim.prop") <- length(outliers) / length(untrimmed.bootstraps)
  attr(ret,"untrimn") <- length(untrimmed.bootstraps)
  attr(ret,"outliers") <- length(outliers) 
  attr(ret,"trim.ind") <- !(untrimmed.bootstraps %in% outliers)
  attr(ret,"boxplot.coef") <- boxplot.coef

  return(ret)
}
