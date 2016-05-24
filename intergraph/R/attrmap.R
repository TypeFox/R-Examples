#' Network attribute copying/renaming table
#' 
#' Setting and retrieving rules through which network/edge/vertex attributes
#' are copied or renamed when converting network objects.
#' 
#' Different classes for network data use different attribute names to store
#' the same information. Some of the classes define attributes not used by
#' other classes. This function is used to retrieve or set the rules in which
#' these attributes are copied, renamed or dropped when converting network
#' objects.
#' 
#' The rules are stored in a data frame with the following columns (all of mode
#' character):
#' \describe{
#' \item{type}{type, or level of the attribute, one of "vertex", "edge" or
#' "network"}
#' \item{fromcls}{name of the class which is being coerced}
#' \item{fromattr}{name of the "source" attribute}
#' \item{tocls}{name of the class to which coercion is being done}
#' \item{toattr}{name of the attribute in the result of coercion, or NA}
#' }
#' When converting network \code{x}, of class \code{xclass} to network \code{y}
#' of class \code{yclass} the coercion method looks up the rows in this table
#' for which \code{fromcls=xclass} and \code{tocls=yclass}. If network \code{x}
#' has any of the attributes listed in the \code{fromattr} in the resulting
#' subset of the table then, it is renamed to the corresponding value of
#' \code{toattr}.  If \code{toattr} is \code{NA} the attribute is dropped from
#' the result.
#' 
#' @param newdf data.frame, new set of copy/rename rules, see Details
#'
#' @return If \code{newdf} is NULL, a data frame with the current copy/rename
#' rules. Otherwise \code{newdf} are set as new rules and the old ones are
#' returned invisibly, much like \code{\link{par}}.
#'
#' @seealso This is used by \code{\link{asIgraph}} and \code{\link{asNetwork}}
#'
#' @export
#'
#' @examples
#'
#'# Current values
#'attrmap()
#'
attrmap <- function(newdf=NULL)
{
  cur <- utils::getFromNamespace(".attrmap", "intergraph")
  if( is.null(newdf) )
  {
    # return current
    return(cur)
  } else
  {
    # assign new
    utils::assignInNamespace(".attrmap", newdf, "intergraph")
    # return old invisibly
    invisible(cur)
  }
}

# subset of attrmap depending on type and from/to classes
attrmapmat <- function(from, to, atype, db=attrmap())
{
  i <- with(db, type==atype & fromcls==from & tocls==to)
  as.matrix( db[i, c("fromattr", "toattr")] )
}

# attribute renaming table
.attrmap <- matrix(ncol=5, byrow=TRUE, c(
        "network" , "network" , "directed"     , "igraph"  , NA             , 
        "network" , "network" , "bipartite"    , "igraph"  , NA             , 
        "network" , "network" , "loops"        , "igraph"  , NA             , 
        "network" , "network" , "mnext"        , "igraph"  , NA             , 
        "network" , "network" , "multiple"     , "igraph"  , NA             , 
        "network" , "network" , "n"            , "igraph"  , NA             , 
        "network" , "network" , "hyper"        , "igraph"  , NA             , 
        "vertex"  , "igraph"  , "name"         , "network" , "vertex.names"
        ) )
.attrmap <- as.data.frame(.attrmap, stringsAsFactors=FALSE)
names(.attrmap) <- c("type", "fromcls", "fromattr", "tocls", "toattr")
