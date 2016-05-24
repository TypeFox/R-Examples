#' Gene abundance proportion distribution of RNA-seq data.
#'
#' An example of gene abundance proportion distribution function of RNA-seq data,
#' generated from a real dataset. See supplmenetary material of Law \emph{et al}.
#'
#'@name qAbundanceDist
#'@rdname qAbundanceDist
#'@aliases qAbundanceDist
#'
#'@usage data(qAbundanceDist)
#'
#'@references Law CW, Chen Y, Shi W & Smyth GK, voom: Precision
#'weights unlock linear model analysis tools for RNA-seq read counts, \emph{Genome
#'Biology}, 15(2), R29, 2014.
#'
#' @format A function: \code{qAbundanceDist}.
#'
#' @examples
#' \dontrun{
#' # Get distribution function of abundance proportions
#' # This distribution was generated from a real dataset
#' #load(url("http://bioinf.wehi.edu.au/voom/qAbundanceDist.RData"))
#' data("qAbundanceDist")
#' curve(qAbundanceDist, from=0, to =0.99)
#'
#' # Generate baseline proportions for desired number of genes
#' ngenes <- 10000
#' baselineprop <- qAbundanceDist( (1:ngenes)/(ngenes+1) )
#' baselineprop <- baselineprop/sum(baselineprop)
#' }
#'
#'
#' @source \url{http://bioinf.wehi.edu.au/voom/}
#'
#' @keywords datasets
#' @docType data
NULL

