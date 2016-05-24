#' Colour sub-clusters of a tree (dendrogram/hclust) object
#' 
#' The distinctive feature of this function is to colour both the 
#' terminal leaves of a cluster and the edges leading to those leaves.
#' The edgePar attribute of nodes will be augmented by a new list item col.
#' The groups will be defined by a call to \code{\link{slice}} using the k or h
#' parameters.
#' @details If \code{groupLabels=TRUE} then numeric group labels will added to
#'   each cluster. If a vector is supplied then these entries will be used as 
#'   the group labels. If a function is supplied then it will be passed a 
#'   numeric vector of groups (e.g. 1:5) and must return the formatted group
#'   labels.
#' @param d A \code{dendrogram} or \code{hclust} tree object
#' @param k number of groups (passed to \code{slice})
#' @param h height at which to cut tree (passed to \code{slice})
#' @param col Function or vector of colours
#' @param groupLabels If TRUE add numeric group label - see Details for options
#' @return a tree object of class dendrogram.
#' @aliases color_clusters
#' @author jefferis
#' @export
#' @seealso \code{\link{slice},\link{cutree},\link{dendrogram}}
#' @examples
#' d5=colour_clusters(hclust(dist(USArrests), "ave"),5)
#' plot(d5)
#' d5g=colour_clusters(hclust(dist(USArrests), "ave"),5,groupLabels=TRUE)
#' plot(d5g)
#' d5gr=colour_clusters(hclust(dist(USArrests), "ave"),5,groupLabels=as.roman)
#' plot(d5gr)
colour_clusters<-function(d,k=NULL,h=NULL,col=rainbow,groupLabels=NULL){
  # TODO make this more modular
  if(!inherits(d,'dendrogram') && !inherits(d,'hclust'))
    stop("Expects a dendrogram or hclust object")
  g=slice(d,k=k,h=h)
  if(inherits(d,'hclust')) d=as.dendrogram(d)
  
  # we need unique labels later on, so store the original labels if there are
  # duplicates and make the active ones unique
  if(fix_labels<-any(duplicated(orig_labels<-labels(d)))){  
    unique_labels=make.unique(orig_labels)
    labels(d)<-unique_labels
    names(g)=unique_labels
  }
  
  k=max(g)
  if(is.function(col)) col=col(k)
  else if(length(col)!=k) stop("Must give same number of colours as clusters")
  
  if(!is.null(groupLabels)){
    if(length(groupLabels)==1){
      if(is.function(groupLabels))
        groupLabels=groupLabels(seq.int(length.out=k))
      else if(is.logical(groupLabels)){
        if(groupLabels)
          groupLabels=seq.int(length.out=k)
        else groupLabels=NULL
      }
    }
    if(!is.null(groupLabels) && length(groupLabels)!=k)
      stop("Must give same number of group labels as clusters")
  }
  
  addcol<-function(n,col) {
    attr(n,'edgePar')=c(attr(n,'edgePar'),list(col=col))
    n
  }
  
  descendTree<-function(sd){
    groupsinsubtree=unique(g[labels(sd)])
    if(length(groupsinsubtree)>1){
      # keep descending 
      for(i in seq(sd))
        sd[[i]]<-descendTree(sd[[i]])
    } else {
      # else assign colours
      # sd=dendrapply(sd,addcol,col[groupsinsubtree],groupsinsubtree)
      sd=dendrapply(sd,addcol,col[groupsinsubtree])
      if(!is.null(groupLabels)){
        attr(sd,'edgetext')=groupLabels[groupsinsubtree]
        attr(sd,'edgePar')=c(attr(sd,'edgePar'),list(p.border=col[groupsinsubtree]))
      }
    }
    sd
  }
  d=descendTree(d)
  if(fix_labels) labels(d)=orig_labels
  d
}

#' @export
#' @rdname colour_clusters
color_clusters<-colour_clusters

#' Return the leaf colours of a dendrogram
#' 
#' @details The returned colours will be in dendrogram order.
#' @param d the dendrogram
#' @param col_to_return Character scalar - kind of colour attribute to return
#' @return named character vector of colours, NA_character_ where missing
#' @author jefferis
#' @export
#' @aliases leaf_colors
#' @seealso \code{\link{slice},\link{colour_clusters}}
#' @examples
#' d5=colour_clusters(hclust(dist(USArrests), "ave"),5)
#' leaf_colours(d5)
leaf_colours<-function(d,col_to_return=c("edge",'node','label')){
  if(!inherits(d,'dendrogram')) stop("I need a dendrogram!")
  col_to_return=match.arg(col_to_return)
  leaf_col<-function(n,col_to_return) {
    if(is.leaf(n)) {
      col=switch(col_to_return,
          edge=attr(n,'edgePar')$col,
          node=attr(n,'nodePar')$col,
          label=attr(n,'nodePar')$lab.col)
      if(is.null(col)) col=NA_character_
      structure(col,.Names=attr(n,'label'))
    } else NULL
  }
  
  unlist(dendrapply(d,leaf_col,col_to_return))
}

leaf_colors<-leaf_colours

#' Set the leaf colours of a dendrogram
#' 
#' @param d the dendrogram
#' @param col Single colour or named character vector of colours. When NA no
#'   colour will be set.
#' @param col_to_set Character scalar - kind of colour attribute to set
#' @author jefferis
#' @export
#' @aliases set_leaf_colors
#' @seealso \code{\link{slice},\link{colour_clusters}}
#' @examples
#' d5=colour_clusters(hclust(dist(USArrests), "ave"),5)
#' dred=set_leaf_colours(d5,'red','edge')
#' stopifnot(isTRUE(all(leaf_colours(dred)=='red')))
#' d52=set_leaf_colours(d5,leaf_colours(d5),'edge')
#' stopifnot(all.equal(d5,d52))
set_leaf_colours<-function(d,col,col_to_set=c("edge",'node','label')){
  if(!inherits(d,'dendrogram')) stop("I need a dendrogram!")
  d_labels=labels(d)
  if(length(col)==1) {
    # special case, we want to set all the colour the same - doesn't need to be
    # named (indeed makes no sense to be named)
    col=structure(rep(col,length(d_labels)),.Names=d_labels)
  } else {
    # check the colour vector is sensibly named
    if(is.null(names(col)) || !all(names(col) %in% d_labels))
      stop("col must be a vector of colours _named_ by existing node labels")
  }
  col_to_set=match.arg(col_to_set)
  leaf_col<-function(n,col,col_to_set) {
    if(is.leaf(n)) {
      this_col=unname(col[attr(n,'label')])
      if(is.na(this_col)) {
        # do nothing
      } else if(col_to_set=='edge'){
        attr(n,'edgePar')$col=this_col
      } else if(col_to_set=='node') {
        attr(n,'nodePar')$col=this_col
      } else if(col_to_set=='label') {
        attr(n,'nodePar')$lab.col=this_col
      }
    }
    n
  }
  dendrapply(d,leaf_col,col,col_to_set)
}

#' @export
#' @rdname set_leaf_colours
set_leaf_colors<-set_leaf_colours
