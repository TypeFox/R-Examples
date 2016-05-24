#'Function to draw a gapped heatmap
#'
#'This function draws a gapped heatmap using the ggplot2 package. The input for the
#'function are the gapdata class objects, generated from gap_data() function, and the data matrix.
#'
#' @param m data matrix
#' @param row_gap a gapdata class object for rows
#' @param col_gap a gapdata class object for columns
#' @param row_labels a logical to show labels for rows
#' @param col_labels a logical to show lables for columns
#' @param rotate a logical to rotate row labels
#' @param col colors used for heatmap
#' @param ... ignored
#' @export gap_heatmap
#' @aliases gap_heatmap
#' @return a ggplot object
#' 


gap_heatmap <- function(m, row_gap=NULL, col_gap=NULL, row_labels=TRUE, col_labels=TRUE, rotate=FALSE,
                        col=c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", 
                              "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")){
  sort_row = FALSE
  r_order = rownames(m)
  if(!is.null(row_gap)){
    #check row size
    if(nrow(m) == nrow(row_gap$labels)){
      sort_row = TRUE
      #reorder matrix row
      r_order <- match(row_gap$labels$label, rownames(m))
    }else{
      stop(paste0("The numbers of rows between the matrix and row_gapdata are different. matrix:",nrow(m)," row:", nrow(row_gap$labels)))
      return(NULL)
    }
  }
  sort_col = FALSE
  c_order = colnames(m) 
  if(!is.null(col_gap)){
    #check col size
    if(ncol(m) == nrow(col_gap$labels)){
      sort_col = TRUE
      #reorder matrix col
      c_order <- match(col_gap$labels$label, colnames(m))
    }else{
      stop(paste0("The numbers of columns between the matrix and col_gapdata are different. matrix:",ncol(m)," col:", nrow(col_gap$labels)))
      return(NULL)
    }
  }
  #reorder matrix
  M <- m[r_order,c_order]
  
  #change the colum names based on the actual position
  if(sort_col){
    colnames(M) <- as.character(col_gap$labels$x)
  }else{
    colnames(M) <- colnames(m)
  }
  
  if(sort_row){
    rownames(M) <- as.character(row_gap$labels$x)
  }else{
    rownames(M) <- rownames(m)
  }
  M <- melt(M) #flatten matrix
  #change the column names
  colnames(M) <- c("y", "x", "value")
  
  if(!row_labels & !col_labels){
    #no decoration
    p <- ggplot(data = M, aes_string(x="x",y="y")) + geom_tile(aes_string(fill="value"))+ scale_fill_gradientn(colours = col) + 
      theme(
        plot.margin = grid::unit(c(0,0, grid::unit(-0.4, "line"), grid::unit(-0.4, "line")), "lines"), 
        #plot.margin = grid::unit(c(0,0,0,0), "lines"), 
        #panel.margin = grid::unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        line = element_blank(),
        legend.position="none")+ 
      labs(x=NULL, y=NULL, title=NULL) +
      scale_x_continuous(expand=c(0,0))+ 
      scale_y_continuous(expand=c(0,0)) 
  }else{
    p <- ggplot(data = M, aes_string(x="x",y="y")) + geom_tile(aes_string(fill="value"))+ scale_fill_gradientn(colours = col) + 
      labs(title=NULL) +
      theme(
        plot.margin = grid::unit(c(0,0,0,0), "lines"), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        line = element_blank(),
        legend.position="none")
      
    if(row_labels & !is.null(row_gap)){
      p <- p + scale_y_continuous(breaks = row_gap$labels$x, labels = row_gap$labels$label)
      if(rotate){
        p <- p + theme(axis.text.y =element_text(color = "black", angle=90))
      }else{
        p <- p + theme(axis.text.y =  element_text(color = "black"))
      }
    }else if(!row_labels | is.null(row_gap)){
      p <- p + labs(y=NULL) + theme(axis.text.y=element_blank())
    }
    
    if(col_labels & !is.null(col_gap)){
      p <- p +scale_x_continuous(breaks = col_gap$labels$x, labels = col_gap$labels$label)  + theme(axis.text.x =  element_text(color = "black"))
    }else{
      p <- p + labs(x=NULL) + theme(axis.text.x = element_blank())
    }
  }
  p #return
}