#' A function to cache the phylogenies in treebase locally
#' 
#' @param file filename for the cache, otherwise created with datestamp
#' @param pause1 number of seconds to hesitate between requests
#' @param pause2 number of seconds to hesitate between individual files
#' @param attempts number of attempts to access a particular resource
#' @param max_trees maximum number of trees to return (default is Inf)
#' @param only_metadata  option to only return metadata about matching trees
#' @param save logical indicating whether to save a file with the resuls.
#' @return saves a cached file of treebase
#' @details it's a good idea to let this run overnight
#' @examples \dontrun{
#'  treebase <- cache_treebase()
#' }
#' @export
cache_treebase <- function(file=paste("treebase-", Sys.Date(), ".rda",sep=""),
                           pause1 = 3, pause2 = 3, attempts = 10, 
                           max_trees = Inf, only_metadata = FALSE, save=TRUE){
  treebase <- search_treebase(c("Consensus or Single"), 
                              by=c("type.tree", "type.tree"), 
                              max_trees = max_trees, pause1 = pause1, 
                              pause2 = pause2, attempts = attempts, 
                              only_metadata = only_metadata)
  if(save)
    save("treebase", file=file)
  treebase
}

