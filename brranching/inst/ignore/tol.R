#' Query Treebase for phylogenetic trees
#'
#' @importFrom rotl tnrs_match_names tol_induced_subtree
#' @export
#' @param taxa A list or character vector of taxon names
#'
#' @param context_name name of the taxonomic context to be searched (length-one character
#' vector)
#' @param do_approximate_matching A boolean indicating whether or not to perform
#' approximate string (a.k.a. "fuzzy") matching. Will greatly improve speed if this is
#' turned OFF (false). By default, however, it is on (true).
#' @param ids An array of ids to use for identifying names. These will be assigned to each
#' name in the names array. If ids is provided, then ids and names must be identical in length.
#' @param include_deprecated A boolean indicating whether or not to include deprecated taxa in
#' the search.
#' @param include_dubious Whether to include so-called 'dubious' taxa, those which are not
#' accepted by OTT.
#' @param parser
#' @param ... Further args passed on to \code{\link[rotl]{tol_induced_subtree}}
#' @return One or more trees of class \code{phylo}
#' @references \url{https://github.com/OpenTreeOfLife/opentree/wiki/Open-Tree-of-Life-APIs}
#' @examples \dontrun{
#' taxa <- c("Poa annua", "Phlox diffusa", "Helianthus annuus")
#' tree <- tol(taxa)
#' plot(tree, no.margin=TRUE)
#' }

tol <- function(taxa, context_name=NULL, do_approximate_matching=TRUE, ids=NULL,
  include_deprecated=FALSE, include_dubious=FALSE, parser="rncl", ...)
{
  res <- tnrs_match_names(taxa, context_name=context_name, do_approximate_matching=do_approximate_matching,
                          ids=ids, include_deprecated=include_deprecated, include_dubious=include_dubious)
  tol_induced_subtree(ott_ids=res$ott_id, parser=parser)
}
