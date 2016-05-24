#' generate a table of all available metadata for TreeBASE entries 
#'
#' @param phylo.md cached phyloWS (tree) metadata, (optional)
#' @param oai.md cached OAI-PMH (study) metadata (optional)
#' @return a data frame of all available metadata, (as a data.table object)
#' columns are: "Study.id", "Tree.id", "kind", "type", "quality", "ntaxa"    
#' "date", "publisher", "author", "title".  
#' @examples \dontrun{
#' meta <- metadata()
#' library("data.table")
#' meta[publisher %in% c("Nature", "Science") & ntaxa > 50 & kind == "Species Tree",]
#' }
# @importFrom reshape2 melt dcast
# @suggest data.table
#' @export
metadata <- function(phylo.md = NULL, oai.md=NULL){

  if(is.null(phylo.md))
    phylo.md <- cache_treebase(only_metadata=TRUE, save=FALSE)
  if(is.null(oai.md))
    oai.md <- download_metadata() 

  #  phylo_data <- melt(phylo.md) # This is very slow
  #  phylo <- dcast(phylo_data, ... ~ L2)

  phylo <- 
    data.frame(Study.id = phylo_metadata("Study.id", phylo.md), 
    Tree.id = phylo_metadata("Tree.id", phylo.md), 
    kind = phylo_metadata("kind", phylo.md), 
    type = phylo_metadata("type", phylo.md), 
    quality = phylo_metadata("quality", phylo.md),
    ntaxa = as.numeric(phylo_metadata("ntaxa", phylo.md)))

  oai <- 
    data.frame(date = as.integer(oai_metadata("date", oai.md)), 
    publisher = oai_metadata("publisher", oai.md),
    author = oai_metadata("author", oai.md),
    title = oai_metadata("title", oai.md),
    Study.id = oai_metadata("Study.id", oai.md))
  both <- merge(phylo, oai)
  class(both$ntaxa) <- "numeric"
  both
}

# Helper functions 

#' Search the PhyloWS metadata 
#' 
#' @param x one of "Study.ids", "Tree.ids", "kind", "type", "quality", "ntaxa"
#' @param metadata returned from \code{search_treebase} function. 
#' if not specified will download latest copy of PhyloWS metadata from treebase.  Pass
#' in search results value during repeated calls to speed function runtime substantially
#' @param ... additional arguments to \code{search_treebase}
#' @return a list of the values matching the query
#' @keywords internal
# @examples
# \dontrun{
#      # calls will work without a metadata object, but requires longer to download data
#      kind <- phylo_metadata("kind")
#      type <- phylo_metadata("type") 
#      table(kind, type)
#      }
#      # but are much faster if the data object is provided, see cache_treebase():
#      data(treebase)
#      kind <- phylo_metadata("kind", metadata=treebase)
#      type <- phylo_metadata("type", metadata=treebase) 
#      table(kind, type)
# 
phylo_metadata <- function(x =  c("Study.id", "Tree.id", "kind", "type", "quality", "ntaxa"), metadata=NULL, ...){
  x = match.arg(x)
# Aliases
  x <- gsub("Study.id", "S.id", x)
  x <- gsub("Tree.id", "Tr.id", x)
  x <- gsub("ntaxa", "ntax", x)
  if(is.null(metadata))
    metadata <- cache_treebase(only_metadata=TRUE, save=FALSE)
  sapply(metadata, `[[`, x)
}


#' Search the OAI-PMH metadata by date, publisher, or identifier
#' 
#' @param x one of "date", "publisher", "identifier" for the study
#' @param metadata returned from \code{download_metadata} function. 
#' if not specified will download latest copy from treebase.  Pass
#' in the value during repeated calls to speed function runtime substantially
#' @param ... additional arguments to \code{download_metadata}
#' @return a list of values matching the query
#' @keywords internal
# @examples
# \dontrun{
#     # automatically search each time (note: this calls internal function now)
#     dates <- oai_metadata("date") 
#     pub <- oai_metadata("publisher")
#     table(dates, pub)
# }
#    # Using cached data from an earlier download
#     data(metadata) 
#     dates <- oai_metadata("date", metadata=metadata) 
#     pub <- oai_metadata("publisher", metadata=metadata)
#     table(dates, pub)
oai_metadata <- function(x = c("date", "publisher", "author", "title", "Study.id", "attributes"), metadata=NULL, ...){
  x = match.arg(x)
# Aliases
  x <- gsub("attributes", ".attr", x)
  x <- gsub("author", "creator", x) # gets first author only
  x <- gsub("Study.id", "identifier", x)
  if(is.null(metadata))
    metadata <- download_metadata(...)
  out <- sapply(metadata, `[[`, x)
  if(x == "identifier")
    out <- gsub(".*TB2:S(\\d*)", "\\1", out)
  as.character(out) # avoid list object returns
}



