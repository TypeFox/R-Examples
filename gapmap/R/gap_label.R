#'Function to draw a gapped labels
#'
#'This function draws a gapped labels using the ggplot2 package. The input for the
#'function is the gapdata class object, generated from gap_data() function.
#'
#' @param data gapdata class object
#' @param orientation orientation of the labels, "left", "top", "right", or "bottom"
#' @param label_size a numeric to set the label text size
#' @export gap_label
#' @aliases gap_label
#' @return a ggplot object
#' 

gap_label <- function(data, orientation, label_size = 5){
  #if data is not gapmap
  if (!is.gapdata(data)){
    stop("data is not a gapdata object")
  } 
  if(is.null(data)){
    stop("you need to provide a gapdata object.")
  }
  if(is.null(orientation)){
    stop("you need to indicate the orientation: left, top, right or bottom")
  }
  
  #position labels
  x_min = min(data$labels$x) - 0.5
  x_max = max(data$labels$x) + 0.5
  0
  if(orientation  == "left"){
    p <- ggplot(data=data$label, aes_string(x="y", y="x", label="label") ) +
      scale_y_continuous(expand=c(0,0), limits = c(x_min, x_max))+ scale_x_continuous(expand=c(0,0), limits=c(-1, 0)) +
      geom_text(hjust=1, size = label_size) + theme(plot.margin = grid::unit(c(0, 0, -0.4, -0.2), "lines"))
  }else if(orientation =="top"){
    p <- ggplot(data=data$label, aes_string(x="x", y="y", label="label") ) +
      scale_x_continuous(expand=c(0,0), limits = c(x_min, x_max))+ scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
      geom_text(angle = 90, hjust=0, size = label_size) + theme(plot.margin = grid::unit(c(0.2, 0, -0.4, -0.4), "lines"))
  }else if(orientation=="right"){
    p <- ggplot(data=data$label, aes_string(x="y", y="x", label="label") ) +
      scale_y_continuous(expand=c(0,0), limits = c(x_min, x_max))+ scale_x_continuous(expand=c(0,0), limits=c(0,1)) +
      geom_text(hjust=0, size = label_size) + theme(plot.margin = grid::unit(c(0.2, 0, -0.4, -0.2), "lines"))
  }else if(orientation =="bottom"){
    p <- ggplot(data=data$label, aes_string(x="x", y="y", label="label") ) +
      scale_x_continuous(expand=c(0,0), limits = c(x_min, x_max))+ scale_y_continuous(expand=c(0,0), limits=c(-1,0)) +
      geom_text(angle = 90, hjust=1, size = label_size) + theme(plot.margin = grid::unit(c(0.2, 0, -0.4, -0.4), "lines"))
  }
 
  p <- p + theme( 
      #plot.margin = grid::unit(c(0,0,0,0), "lines"),               
      panel.margin = grid::unit(c(0,0,0,0), "lines"),
      panel.background = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(),
      axis.text.x =  element_blank(),
      axis.text.y =  element_blank(),
      axis.line=element_blank(),
      axis.line.x=element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.x=element_blank(),
      axis.ticks.y=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank()) + labs(x=NULL, y=NULL)
  
  p #return
}
