#######################################################################
##
## Function: trim.data() [ formerly subset.anchors() ]
## Author  : Jonathan Wand (wand(at)stanford.edu)
## Created :  2002-09-02
##
## DESCRIPTION:
##   Trim a dataset to have same cases present in an anchors.data object 
##
## INPUT:
##   data : data.frame
##   anchors: an anchors.$method object
##   
## OUTPUT:
##   new data.frame, subsetted
##
## MODIFIED:
##   2008-05-01
##   - trims for all anchors.data objects (except chopit, which is ambiguous
##     without delete=listwise)
##   - removes 'insert' option
## 
#######################################################################

trim.data <- function( data , anchors) {

  if ( !( is.data.frame(data) | is.matrix(data) ) ) 
    stop("trim.data needs first argument to be a data.frame or matrix")

  if (class(anchors)!="anchors.data" && !is.null(anchors$data)){
    anchors <- anchors$data
  }

  
#  if (!( anchors$method %in% c("order","rank","cpolr","minentropy","entropy") || anchors$delete == "listwise" ))
  if (( anchors$method == "chopit" && anchors$delete != "listwise" ))
  stop(paste("trim.data() not defined for anchors.data objects created with \n",
             "combination method='chopit' and options$delete='minimal'"))

  idx <- rownames(data) %in% rownames(anchors$z0)

  if (sum(idx)==0) {
    warning("no rownames match between original data and anchors.data object")
  }
  
  return(data[idx,])
}

