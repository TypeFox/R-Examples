
#' taxize nexml
#'
#' Check taxanomic names against the specified service and 
#' add appropriate semantic metadata to the nexml OTU unit 
#' containing the corresponding identifier. 
#' @param nexml a nexml object
#' @param type the name of the identifier to use
#' @param ... additional arguments (not implemented yet)
#' @import taxize
#' @export 
#' @examples \dontrun{
#' data(bird.orders)
#' birds <- add_trees(bird.orders)
#' birds <- taxize_nexml(birds, "NCBI")
#' }
taxize_nexml <- function(nexml, type = c("NCBI"), ...){
          type <- match.arg(type)
          if(type == "NCBI"){
            for(j in 1:length(nexml@otus)){
              for(i in 1:length(nexml@otus[[j]]@otu)){
                id <- get_uid(nexml@otus[[j]]@otu[[i]]@label)
                if(is.na(id))
                  warning(paste("ID for otu", nexml@otus[[j]]@otu[[i]]@label, "not found. Consider checking the spelling or alternate classification"))
                else 
                  nexml@otus[[j]]@otu[[i]]@meta <- new("ListOfmeta", list(
                                 meta(href = paste0("http://ncbi.nlm.nih.gov/taxonomy/", id),
                                      rel = "tc:toTaxon")))

              }
            }
          }
  nexml
}
