#'Extract a list from the dendrogram object
#'
#'This function extract list of data.frames for drawing dendrograms
#'
#' @param d dendrogram class object
#' @param type either "triangular" or "rectangular". It determines the same of branches.
#' @param segments_df data.frame storing the segment information
#' @param labels_df data.frame storking the label positions
#' @param ... ignored
#' @export extract_list
#' @aliases extract_list
#' @return the extracted list
#' @keywords internal
#' 

#get list of drawing parameter, like gg.plotNode
extract_list <- function(d, type, segments_df=NULL, labels_df=NULL,...){
  inner <- !is.leaf(d) #check if it is subtree
  yTop <- attr(d, "height") #height of subtree
  xMid <- attr(d, "xmid")
  
  if (is.leaf(d)) { #singleton cluster
    Y <- yTop 
    X <- attr(d,"xpos")
    nodeText <- as.character(attr(d, "label"))
    labels_df <- rbind(labels_df, data.frame(x=X, y=0, text=nodeText))
    
  }else if (inner) { #subtree
    for (k in seq_along(d)) {
      child <- d[[k]]
      yBot <- attr(child, "height")
      if (is.null(yBot)) {
        yBot <- 0
      }
      xBot <- 0
      if(is.leaf(child)){
        xBot <- attr(child, "xpos")
      }else{
        xBot <- attr(child, "xmid")
        if(is.null(attr(child,"xmid"))) print("xmid is nulll")
      }  
      if (type == "triangle") {
        #adding lines
        segments_df <- rbind(segments_df, get_segment_df(xMid, yTop, xBot, yBot))
      } else {
        segments_df <- rbind(segments_df, get_segment_df(xMid, yTop, xBot, yTop))
        segments_df <- rbind(segments_df, get_segment_df(xBot, yTop, xBot, yBot))
      }
      #call recursively
      temp_list <- extract_list(d = child, type, segments_df, labels_df)
      segments_df <- temp_list$segments
      labels_df <- temp_list$labels
    }
  }
  l =list(segments=segments_df, labels=labels_df) 
  return(l)
}