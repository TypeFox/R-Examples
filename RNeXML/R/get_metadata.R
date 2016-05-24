
## FIXME might want to define this for sub-nodes.  e.g. so we can get all metadata on "nodes" in tree2...

#' get_metadata
#' 
#' get_metadata 
#' @param nexml a nexml object
#' @param level the name of the level of element desired, see details
#' @return the requested metadata as a data.frame. Additional columns
#' indicate tha parent element of the return value.
#' @details 'level' should be either the name of a child element of a NeXML document 
#' (e.g. "otu", "characters"), or a path to the desired element, e.g. 'trees/tree'
#' will return the metadata for all phylogenies in all trees blocks.
#' @import XML
#' @examples \dontrun{
#' comp_analysis <- system.file("examples", "primates.xml", package="RNeXML")
#' nex <- nexml_read(comp_analysis)
#' get_metadata(nex)
#' get_metadata(nex, "otus/otu")
#' }
#' @export
get_metadata <- function(nexml, level = "nexml"){
  
#  level = c("nexml", "otus", "trees", "characters", 
#            "otus/otu", "trees/tree", "characters/format", "characters/matrix",
#            "characters/format/states")
#  level <- match.arg(level)

  ## Handle deprecated formats
  if(level =="otu")
    level <- "otus/otu"
  if(level =="tree")
    level <- "trees/tree"

  
  
  if(level == "nexml")
    level <- "meta"
  else
    level <- paste(level, "meta", sep="/") 
 
  get_level(nexml, level)
  

}
