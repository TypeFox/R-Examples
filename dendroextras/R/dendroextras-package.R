#' Extra functions to cut, label and colour dendrogram clusters
#' 
#' @section Cutting clusters:
#'   
#'   \code{dendroextras} provides the \code{\link{slice}} function as an 
#'   alternative to the base \code{\link{cut}} function. In contrast to 
#'   \code{cut}, \code{slice} returns group membership \emph{in dendrogram 
#'   order} i.e. the first element in the group vector that is returned will be 
#'   the leftmost member of the leftmost cluster (cluster #1).
#'   
#' @section Colouring clusters:
#'   
#'   \code{dendroextras} provides \code{\link{colour_clusters}} to colour all of
#'   the edges forming clusters cut by height or number of groups. You can also 
#'   set and retrieve the leaf colours (i.e. the terminal nodes) using 
#'   \code{\link{set_leaf_colours}} and \code{\link{leaf_colours}}.
#'   
#' @section Labels:
#'   
#'   \code{dendroextras} provides \code{\link{labels}} and
#'   \code{\link{labels<-}} methods to get and set the labels of cluster
#'   members.
#'   
#' @name dendroextras-package
#' @aliases dendroextras
#' @seealso \code{\link{dendrogram}, \link{hclust}} in \code{\link{stats}} 
#'   package.
#' @docType package
#' @keywords package, dendrogram
NULL
