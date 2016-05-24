#' Write nexml files
#' 
#' @param x a nexml object, or any phylogeny object (e.g. phylo, phylo4) 
#' that can be coerced into one. Can also be omitted, in which case a new 
#' nexml object will be constructed with the additional parameters specified.
#' @param file the name of the file to write out
#' @param trees phylogenetic trees to add to the nexml file (if not already given in x)
#' see \code{\link{add_trees}} for details.  
#' @param characters additional characters
#' @param meta A meta element or list of meta elements, see \code{\link{add_meta}}
#' @param ... additional arguments to add__basic_meta, such as the title.  See \code{\link{add_basic_meta}}.   
#' @return Writes out a nexml file
#' @import ape
#' @import XML 
#' @import methods
#' @aliases nexml_write write.nexml
#' @export nexml_write write.nexml
#' @seealso \code{\link{add_trees}} \code{\link{add_characters}} \code{\link{add_meta}} \code{\link{nexml_read}}

#' @examples
#'  ## Write an ape tree to nexml, analgous to write.nexus:
#'  library(ape); data(bird.orders)
#'  write.nexml(bird.orders, file="example.xml")
#' 
#' \dontrun{ # takes > 5s
#'  ## Assemble a nexml section by section and then write to file:
#'  library(geiger)
#'  data(geospiza)
#'  nexml <- add_trees(geospiza$phy) # creates new nexml
#'  nexml <- add_characters(geospiza$dat, nexml = nexml) # pass the nexml obj to append character data
#'  nexml <- add_basic_meta(title="my title", creator = "Carl Boettiger", nexml = nexml)
#'  nexml <- add_meta(meta("prism:modificationDate", format(Sys.Date())), nexml = nexml)
#'
#'  write.nexml(nexml, file="example.xml")
#'
#'  ## As above, but in one call (except for add_meta() call).  
#'  write.nexml(trees = geospiza$phy, 
#'              characters = geospiza$dat, 
#'              title = "My title", 
#'              creator = "Carl Boettiger",
#'              file = "example.xml")
#'  
#'  ## Mix and match: identical to the section by section: 
#'  nexml <- add_meta(meta("prism:modificationDate", format(Sys.Date())))
#'  write.nexml(x = nexml,
#'              trees = geospiza$phy, 
#'              characters = geospiza$dat, 
#'              title = "My title", 
#'              creator = "Carl Boettiger",
#'              file = "example.xml")
#' 
#' }
nexml_write <- function(x = new("nexml"),
                        file = NULL,
                        trees = NULL,
                        characters = NULL,
                        meta = NULL, 
                        ...){
  
  nexml <- as(x, "nexml")
  if(!is.null(trees))
    nexml <- add_trees(trees, nexml = nexml)
  if(!is.null(characters))
    nexml <- add_characters(characters, nexml = nexml)
  if(!is.null(meta))
    nexml <- add_meta(meta, nexml = nexml)
  nexml <- do.call(add_basic_meta, c(list(...), list(nexml=nexml)))
  
  out <- as(nexml, "XMLInternalNode")
  saveXML(out, file = file)
}

write.nexml <- nexml_write


############## Promotion methods ########
## FIXME -- Coercion is not the way to go about any of this


## want generator methods that can handle id creation better
# consider:
# setMethod("promote", 
#           signature("tree", "character"),
#           function(object, target_type)

setAs("tree", "nexml", function(from){
  trees = as(from, "trees")
  otus = as(from, "otus")
  otus@id = "tax1" #UUIDgenerate()
  trees@id = "Trees" #UUIDgenerate()
  trees@otus = otus@id
  new("nexml", 
      trees = new("ListOftrees", list(trees)),
      otus = otus)
})


setAs("ListOfnode", "otus", function(from)
  new("otus", otu = from))

setAs("tree", "trees", function(from)
  new("trees", tree = new("ListOftree", list(from))))


