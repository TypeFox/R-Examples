#' Set the labels of an object
#' @param x Object on which to set labels
#' @param value New labels
#' @param ... Additional parameters passed to specific methods
#' @usage labels(x,...) <- value
#' @export
#' @rdname labels-assign
`labels<-`<-function(x, ..., value)
  UseMethod("labels<-")

#' Set the labels of a dendrogram
#' 
#' @method labels<- dendrogram
#' @usage labels(x,...) <- value
#' @return object of class dendrogram
#' @author jefferis
#' @seealso \code{\link{dendrogram},\link{labels}}
#' @rdname labels-assign
#' @export
#' @examples
#' hc <- hclust(dist(USArrests), "ave")
#' dend <- as.dendrogram(hc)
#' labels(dend)<-abbreviate(labels(dend),minlength=2)
`labels<-.dendrogram`<-function(x, ..., value){
  i=1
  replaceLabel<-function(n){
    if(is.leaf(n)){
      attr(n,'label')=value[i]
      i<<-i+1
    }
    n
  }
  dendrapply(x,replaceLabel)
}

#' Find labels of hclust object (in dendrogram order)
#' 
#' NB will return labels in dendrogram order, not in the
#' order of the original labels retained in object$labels
#' ususally corresponding to the row or column names of 
#' the \code{\link{dist}} object provided to \code{\link{hclust}}.
#' @method labels hclust
#' @param object hclust object from which to extract labels
#' @param ... Additional arguments (ignored)
#' @return character vector of labels in dendrogram order
#' @author jefferis
#' @export
#' @seealso \code{\link{labels},\link{hclust}}
#' @examples
#' hc <- hclust(dist(USArrests), "ave")
#' dend <- as.dendrogram(hc)
#' stopifnot(all.equal(labels(hc),labels(dend)))
labels.hclust<-function(object,...){
  object$labels[object$order]
}
