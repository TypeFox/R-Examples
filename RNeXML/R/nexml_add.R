

#' add elements to a new or existing nexml object
#' 
#' add elements to a new or existing nexml object
#' @param x the object to be added
#' @param nexml an existing nexml object onto which the object should be appended
#' @param type the type of object being provided.  
#' @param ... additional optional arguments to the add functions
#'
#' @return a nexml object with the additional data
#' @seealso \code{\link{add_trees}} \code{\link{add_characters}} \code{\link{add_meta}} \code{\link{add_namespaces}}    
#' @export 
#' @examples
#' library("geiger")
#' data(geospiza)
#' geiger_nex <- nexml_add(geospiza$phy, type="trees")
#' geiger_nex <- nexml_add(geospiza$dat, nexml = geiger_nex, type="characters")
nexml_add <- function(x, 
                      nexml = new("nexml"), 
                      type = c("trees", 
                               "characters", 
                               "meta", 
                               "namespaces"),
                      ...){

  switch(type,
         "trees" = add_trees(x, nexml, ...),
         "characters" = add_characters(x, nexml, ...),
         "metadata" = add_meta(x, nexml, ...), 
         "namespaces" = add_namespaces(x, nexml, ...))

}
