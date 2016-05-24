#'Function to draw a gapped dendrogram
#'
#'This function draws a gapped dendrogram using the ggplot2 package. The input for the
#'function is the gapdata class object, generated from gap_data() function.
#'
#' @param data gapdata class object
#' @param leaf_labels a logical to show labels or not
#' @param rotate_label a logical to rotate labels or not
#' @param orientation a character to set the orientation of dendrogram. Choices are "top", "right", "bottom", "left".
#' @param ... ignored
#' @export gap_dendrogram
#' @aliases gap_dendrogram
#' @return a ggplot object
#' 


gap_dendrogram <- function (data, leaf_labels = TRUE, rotate_label=FALSE, orientation =c("top", "right", "bottom", "left"), ...) {
  dataClass <- class(data)
  orientation <- match.arg(orientation)
  #if data is not gapmap
  if (!is.gapdata(data)){
    stop("data is not a gapdata object")
  } 
  #start ggplot
  p <- ggplot2::ggplot()
  #dendrogram lines
  if (is.null(data$segments)) {
    stop("no segments information in the gapmap object.")
  }
  
  #position labels
  x_min = min(data$labels$x) - 0.5
  x_max = max(data$labels$x) + 0.5
  
  if(leaf_labels){
    if(orientation == "top"){
      p <- p + geom_segment(data = data$segments, aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) + 
        scale_x_continuous(breaks = data$labels$x, labels = data$labels$label)
    }else if(orientation =="right"){
      p <- p + geom_segment(data = data$segments, aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) +
        scale_x_continuous(breaks = data$labels$x, labels = data$labels$label)+ coord_flip()
    }else if(orientation =="bottom"){
      p <- p + geom_segment(data = data$segments, aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) + 
        scale_x_continuous(breaks = data$labels$x, labels = data$labels$label) + scale_y_reverse()
    }else if(orientation == "left"){
      p <- p + geom_segment(data = data$segments, aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) +
        scale_x_continuous(breaks = data$labels$x, labels = data$labels$label)+ coord_flip() +scale_y_reverse()
    }else{
      stop("invalid orientation parameter.")
    }
  }else{
    #no labels
    if(orientation == "top"){
      p <- p + geom_segment(data = data$segments, aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) + 
        scale_x_continuous(expand=c(0,0), limits = c(x_min, x_max))+ scale_y_continuous(expand=c(0,0)) +
        theme(plot.margin = grid::unit(c(0.2, 0, -0.4, -0.4), "lines"))
    }else if(orientation =="right"){
      p <- p + geom_segment(data = data$segments, aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) +
        scale_x_continuous(expand=c(0,0), limits = c(x_min, x_max))+ scale_y_continuous(expand=c(0,0))+ coord_flip()+
        theme(plot.margin = grid::unit(c(0.2, 0, -0.4, -0.4), "lines"))
    }else if(orientation =="bottom"){
      p <- p + geom_segment(data = data$segments, aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) + 
        scale_x_continuous(expand=c(0,0), limits = c(x_min, x_max))+ scale_y_continuous(expand=c(0,0), trans="reverse") +
        theme(plot.margin = grid::unit(c(0.2, 0, -0.4, -0.4), "lines"))
    }else if(orientation == "left"){
      p <- p + geom_segment(data = data$segments, aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) +
        scale_x_continuous(expand=c(0,0), limits = c(x_min, x_max))+ scale_y_continuous(expand=c(0,0), trans="reverse")+ coord_flip()  +
        theme(plot.margin = grid::unit(c(0, 0, -0.4, -0.2), "lines"))
    }else{
      stop("invalid orientation parameter.")
    }
  }
  
  #styling
  if(leaf_labels){
    p <- p + theme( 
      plot.margin = grid::unit(c(0,0,0,0), "lines"),               
      panel.margin = grid::unit(c(0,0,0,0), "lines"),
      panel.background = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(),
      axis.text.x =  element_text(color = "black"), 
      axis.text.y =  element_text(color = "black"), 
      axis.line.y = element_line(colour = "black"),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank()) + labs(x=NULL, y=NULL)
    if(rotate_label){
      p <- p + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
    }
  }else{
    p <- p + theme( 
#       plot.margin = grid::unit(c(0.5,0, grid::unit(-0.3, "line"), grid::unit(-0.5, "line")), "lines"),               
      panel.margin = grid::unit(c(0,0,0,0), "lines"),
      panel.background = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(),
      axis.line=element_blank(),
      axis.text.x=element_blank(),axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),axis.title.y=element_blank()) + labs(x=NULL, y=NULL)
  }
  p #return
}