#' Aggregate data by given taxonomic rank
#'
#' @export
#' @param data A data.frame. Column headers must have capitalized ranks (e.g.,
#'    Genus, Tribe, etc.) (data.frame)
#' @param datacol The data column (character)
#' @param rank Taxonomic rank to aggregate by (character)
#' @param fxn Arithmetic function or vector or functions (character)
#' @examples
#' library("vegan")
#' data(dune.taxon, package='vegan')
#' dat <- dune.taxon
#' set.seed(1234)
#' dat$abundance <- round(rlnorm(n=nrow(dat),meanlog=5,sdlog=2),0)
#' rankagg(data=dat, datacol="abundance", rank="Genus")
#' rankagg(data=dat, "abundance", rank="Family")
#' rankagg(data=dat, "abundance", rank="Genus", fxn="mean")
#' rankagg(data=dat, "abundance", rank="Class")
#' rankagg(data=dat, "abundance", rank="Class", fxn="sd")
rankagg <- function(data=NULL, datacol=NULL, rank=NULL, fxn="sum") {
  if (is.null(data) | is.null(rank))
    stop("You must specify your data.frame and taxonomic rank")
  rank <- match.arg(rank, choices = c('Specices', 'Genus', 'Tribe', 'Family',
                                      'Order', 'Subclass', 'Class', 'Phylum',
                                      'Kingdom'))
  fn <- match.fun(fxn)
  out <- aggregate(data[ , datacol], list(data[ , rank]), fxn)
  names(out) = c(rank, fxn)
  return(out)
}
