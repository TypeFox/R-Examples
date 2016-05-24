
#' Concatenate nexml files 
#' 
#' Concatenate nexml files 
#' @param x,... nexml objects to be concatenated, e.g. from 
#'  \code{\link{write.nexml}} or \code{\link{read.nexml}}. 
#'  Must have unique ids on all elements
#' @param recursive  logical.  If 'recursive = TRUE', the function recursively
#'        descends through lists (and pairlists) combining all their
#'        elements into a vector. (Not implemented).  
#' @return a concatenated nexml file
#' @examples 
#' \dontrun{
#' f1 <- system.file("examples", "trees.xml", package="RNeXML")
#' f2 <- system.file("examples", "comp_analysis.xml", package="RNeXML")
#' nex1 <- read.nexml(f1)
#' nex2 <- read.nexml(f2)
#' nex <- c(nex1, nex2)
#' }
setMethod("c", 
          signature("nexml"), 
          function(x, ..., recursive = FALSE){
              elements = list(x, ...)
              nexml <- new("nexml")
  ## Check that ids are unique
  if(!do.call(unique_ids,elements))
    stop("ids are not unique across nexml files. 
          Consider regenerating ids")
  else {

  nexml@otus <- new("ListOfotus", 
                    unlist(lapply(elements, 
                                  function(n) n@otus), 
                           recursive=FALSE))
  nexml@characters <- new("ListOfcharacters", 
                    unlist(lapply(elements, 
                                  function(n) n@characters), 
                           recursive=FALSE))
  nexml@trees <- new("ListOftrees", 
                    unlist(lapply(elements, 
                                  function(n) n@trees), 
                           recursive=FALSE))
  }
  nexml
})

get_ids <- function(nexml){
  doc <- xmlDoc(as(nexml, "XMLInternalNode"))
  out <- unname(xpathSApply(doc, "//@id"))
  free(doc)
  out
}

unique_ids <- function(...){
  set <- list(...)
  counts <- table(unlist(lapply(set, get_ids)))
  !any(counts > 1)
}

