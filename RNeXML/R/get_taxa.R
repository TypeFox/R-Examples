
#' get_taxa
#'
#' Retrieve names of all species/otus otus (operational taxonomic units) included in the nexml 
#' @aliases get_taxa get_otu
#' @param nexml a nexml object
#' @return the list of taxa
#' @export get_taxa get_otu
#' @examples
#' comp_analysis <- system.file("examples", "comp_analysis.xml", package="RNeXML")
#' nex <- nexml_read(comp_analysis)
#' get_taxa(nex)
#' @seealso  \code{\link{get_item}}
get_taxa <- function(nexml){
  get_level(nexml, "otus/otu")
}
get_otu <- get_taxa

#' get_taxa_list
#' 
#' Retrieve names of all species/otus otus (operational taxonomic units) included in the nexml 
#' @aliases get_taxa_list get_otus_list
#' @param nexml a nexml object
#' @return the list of taxa
#' @export get_taxa_list get_otus_list
#' @seealso  \code{\link{get_item}}
get_taxa_list <-
          function(nexml){
            out <- lapply(nexml@otus, function(otus){
              out <- sapply(otus@otu, function(otu) otu@label)
              names(out) <- name_by_id(otus@otu)
              out
              })
            names(out) <- name_by_id(nexml@otus)
            out
          }
get_otus_list <- get_taxa_list


