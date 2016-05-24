#' Cut a tree-like object into groups numbered in tree order
#' 
#' In comparison with cutree, the groups are numbered from left to right as per
#' the tree when plotted in its standard horizontal form. Note also that the
#' return value will have the leaves sorted in dendrogram order.
#' @param x tree like object
#' @param k an integer scalar with the desired number of groups
#' @param h numeric scalar with height where the tree should be cut
#' @param ... Additional parameters passed to methods
#' @return a named vector with group memberships
#' @author jefferis
#' @export
#' @seealso \code{\link{cutree},\link{cut.dendrogram},\link{rect.hclust}}
#' @examples
#' hc <- hclust(dist(USArrests), "ave")
#' # return groups, leaves ordered by dendrogram
#' slice(hc,k=5)
#' # return groups, leaves ordered as originally passed to hclust
#' slice(hc,k=5)[order(hc$order)]
slice<-function(x,k=NULL,h=NULL,...){
  UseMethod("slice")
}

#' @method slice hclust
#' @export
slice.hclust<-function(x,k=NULL,h=NULL,...){
  # start with results of built-in cutree (where groups are numbered according
  # to original sort order of tree elements)
  chc=cutree(x,k=k,h=h)
  # find the numeric ids of the first listed members of each group
  unique_ids=which(!duplicated(chc))
  # create a numeric vector where the position in the vector indicates new group
  # number i.e. if position 1 contains group 4 then old group 4 -> new group 1
  # NB x$order contains original indices of tree elements ordered by dendrogram
  xtable=order(match(unique_ids, x$order))
  # now map old groups to new
  newgroups=structure(match(chc,xtable),.Names=names(chc))
  # finally return them in dendrogram order
  newgroups[x$order]
}

#' @method slice dendrogram
#' @export
slice.dendrogram<-function(x,k=NULL,h=NULL,...){
  # TODO something a bit more efficient than this!
  slice(as.hclust(x),k=k,h=h,...)
}
