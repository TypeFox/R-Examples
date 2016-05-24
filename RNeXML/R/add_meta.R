#' Add metadata to a nexml file
#' 
#' @param meta a meta S4 object, e.g. ouput of the function \code{\link{meta}}, or a list of these meta objects
#' @param nexml (S4) object
#' @param level the level at which the metadata annotation should be added.
#' @param namespaces named character string for any additional namespaces that should be defined.  
#' @param i for otus, trees, characters: if there are multiple such blocks, which one should be annotated?  Default is first/only block.  
#' @param at_id the id of the element to be annotated.  Optional, advanced use only. 
#' @seealso \code{\link{meta}} \code{\link{add_trees}} \code{\link{add_characters}} \code{\link{add_basic_meta}}
#' @return the updated nexml object
#' @examples
#' ## Create a new nexml object with a single metadata element: 
#' modified <- meta(property = "prism:modificationDate", content = "2013-10-04")
#' nex <- add_meta(modified) # Note: 'prism' is defined in nexml_namespaces by default.  
#' 
#' ## Write multiple metadata elements, including a new namespace:  
#' website <- meta(href = "http://carlboettiger.info", 
#'                 rel = "foaf:homepage")              # meta can be link-style metadata
#' nex <- add_meta(list(modified,  website), 
#'                 namespaces = c(foaf = "http://xmlns.com/foaf/0.1/"))
#'
#' ## Append more metadata, and specify a level: 
#' history <- meta(property = "skos:historyNote",
#'                  content = "Mapped from the bird.orders data in the ape package using RNeXML")
#' nex <- add_meta(history, 
#'                 nexml = nex,
#'                 level = "trees",
#'                 namespaces = c(skos = "http://www.w3.org/2004/02/skos/core#"))
#' 
#' @export add_meta
#' @include classes.R
#' 
add_meta <- function(meta, 
                     nexml=new("nexml"), 
                     level=c("nexml", "otus", "trees", "characters"), 
                     namespaces = NULL,
                     i = 1, 
                     at_id = NULL){
  level <- match.arg(level)
  if(is(meta, "meta"))
    meta <- list(meta)
  if(!all(sapply(meta, is, "meta")))
    stop("All elements in list must be of class 'meta'")
 
  if(!is.null(at_id)){
    stop("function does not yet handle at_id assignments")
    # case not written yet
  } else if(level =="nexml"){ 
    nexml@meta <- new("ListOfmeta", c(unlist(nexml@meta), unlist(meta)))
  } else if(level =="otus"){ 
    nexml@otus[[i]]@meta <- new("ListOfmeta", c(nexml@otus[[i]]@meta, meta))
  }  else if(level =="nexml"){ 
    nexml@trees[[i]]@meta <- new("ListOfmeta", c(nexml@trees[[i]]@meta, meta))
  }  else if(level =="nexml"){ 
    nexml@characters[[i]]@meta <- new("ListOfmeta", c(nexml@characters[[i]]@meta, meta))
  } 

  ## append additional namespaces
  nexml <- add_namespaces(namespaces, nexml)

  nexml
}


