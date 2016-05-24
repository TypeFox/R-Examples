#'Function to draw a gapped cluster heatmap
#'
#'This function draws a gapped cluster heatmap using the ggplot2 package. The input for the
#'function is the a matrix, two dendrograms, and parameters for gaps.
#'
#' @param m matrix
#' @param d_row a dendrogram class object for rows
#' @param d_col a dendrogram class object for columns
#' @param mode gap mode, either "threshold" or "quantitative"
#' @param mapping in case of quantitative mode, either "linear" or "exponential" mapping
#' @param ratio the percentage of width allocated for the sum of gaps.
#' @param scale the sclae log base for the exponential mapping
#' @param threshold the height at which the dendrogram is cut to infer clusters
#' @param row_threshold the height at which the row dendrogram is cut
#' @param col_threshold the height at which the column dendrogram is cut
#' @param rotate_label a logical to rotate column labels or not
#' @param verbose logical for whether in verbose mode or not
#' @param left a character indicating "label" or "dendrogram" for composition
#' @param top a character indicating "label" or "dendrogram" for composition
#' @param right a character indicating "label" or "dendrogram" for composition
#' @param bottom a character indicating "label" or "dendrogram" for composition
#' @param col colors used for heatmap
#' @param h_ratio a vector to set the horizontal ratio of the grid. It should add up to 1. top, center, bottom.
#' @param v_ratio a vector to set the vertical ratio of the grid. It should add up to 1. left, center, right.
#' @param label_size a numeric to set the label text size
#' @param show_legend a logical to set whether to show a legend or not
#' @param ... ignored
#' @export gapmap
#' @aliases gapmap
#' @return a ggplot object
#' @examples
#' set.seed(1234)
#' #generate sample data
#' x <- rnorm(10, mean=rep(1:5, each=2), sd=0.4)
#' y <- rnorm(10, mean=rep(c(1,2), each=5), sd=0.4)
#' dataFrame <- data.frame(x=x, y=y, row.names=c(1:10))
#' #calculate distance matrix. default is Euclidean distance
#' distxy <- dist(dataFrame)
#' #perform hierarchical clustering. default is complete linkage.
#' hc <- hclust(distxy)
#' dend <- as.dendrogram(hc)
#' #make a cluster heatmap plot
#' gapmap(m = as.matrix(distxy), d_row= rev(dend), d_col=dend)
#' 


gapmap <-function(m, d_row, d_col, mode=c("quantitative", "threshold"), mapping=c("exponential", "linear"), 
                  ratio= 0.2, scale = 0.5, threshold = 0, row_threshold = NULL, col_threshold = NULL, rotate_label = TRUE, verbose=FALSE, 
                  left="dendrogram", top="dendrogram", right="label", bottom="label",
                  col=c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"), 
                  h_ratio=c(0.2,0.7,0.1), v_ratio=c(0.2,0.7,0.1),label_size=5, show_legend = FALSE,...){
  #check input
  if(is.null(m) | !inherits(m, "matrix")){
    stop("You need to provide a matrix object.")
  }
  noRowGaps <- FALSE
  if(!is.null(d_row) & !inherits(d_row, "dendrogram")){
    stop("You need to provide a dendrogram object for d_row")
  }else if(is.null(d_row)){
    noRowGaps <- TRUE #no gap for rows
  }
  noColGaps <- FALSE
  if(!is.null(d_col) & !inherits(d_col, "dendrogram")){
    stop("You need to provide a dendrogram object for d_col")
  }else if(is.null(d_col)){
    noColGaps <- TRUE # no gap for columns
  }
  
  #check row and column size
  if(nrow(m) != attr(d_row, "members")){
    stop("The number of rows between the matrix and the d_row is not consistent.")
  }
  if(ncol(m) != attr(d_col, "members")){
    stop("The number of columns between the matrix and the d_col is not consistent.")
  }
  
  mode <- match.arg(mode)
  mapping <- match.arg(mapping)
  
  if(is.null(row_threshold)){
    row_threshold = threshold
  }
  if(is.null(col_threshold)){
    col_threshold = threshold
  }
  
  left_item = NULL
  right_item = NULL
  top_item = NULL
  bottom_item = NULL
  
  #parse Row
  if(noRowGaps){
    row_data <- NULL
  }else{
    row_data <- gap_data(d=d_row, mode=mode, mapping = mapping, ratio = ratio, threshold = row_threshold, verbose = verbose, scale=scale)
  }
  
  #parse Columns
  if(noColGaps){
    col_data <- NULL
  }else{
    col_data <- gap_data(d=d_col,mode=mode, mapping = mapping, ratio = ratio, threshold = col_threshold, verbose = verbose, scale=scale)
  }
  #get all elements
  if(left == "dendrogram"){
    left_item = gap_dendrogram(data = row_data, leaf_labels = FALSE, orientation="left")
  }else if(left =="label"){
    left_item = gap_label(row_data, "left", label_size)    
  }
  
  if(top == "dendrogram"){
    top_item = gap_dendrogram(data = col_data, leaf_labels = FALSE, orientation="top")
  }else if(left =="label"){
    top_item = gap_label(col_data, "top", label_size)    
  }
  
  if(right == "dendrogram"){
    right_item = gap_dendrogram(data = row_data, leaf_labels = FALSE, orientation="right")
  }else if(right =="label"){
    right_item = gap_label(row_data, "right", label_size)    
  }
  
  if(bottom == "dendrogram"){
    bottom_item = gap_dendrogram(data = col_data, leaf_labels = FALSE, orientation="bottom")
  }else if(bottom =="label"){
    bottom_item = gap_label(col_data, "bottom", label_size)    
  }
  
  center_item = gap_heatmap(m, row_gap = row_data, col_gap=col_data, row_labels = FALSE, col_labels = FALSE, col=col)
  
  #extract legend 
  hm <-center_item +theme(legend.position=c(0.5, 0.5), legend.direction = "horizontal")
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(hm)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]]   

  
  #composition of grid layout
  top.height = v_ratio[1]
  center.height = v_ratio[2]
  bottom.height = v_ratio[3]
  left.width = h_ratio[1]
  center.width = h_ratio[2]
  right.width = h_ratio[3]
  
  row.n = 3
  col.n = 3
  
#   if(is.null(left_item)){
#     col.n <- col.n-1
#   }
#   if(is.null(right_item)){
#     col.n <- col.n-1
#   }
#   if(is.null(top_item)){
#     row.n = row.n -1
#   }
#   if(is.null(bottom_item)){
#     row.n = row.n -1
#   }
  
  #start a new blank page
  grid::grid.newpage() 
  #set a Grid layout
  grid_layout <- grid::grid.layout(nrow = row.n, ncol = col.n, widths = grid::unit(c(left.width, center.width, right.width), "null"), heights = grid::unit(c(top.height, center.height, bottom.height), "null"))
  #push viewport to use the layout
  grid::pushViewport(grid::viewport(layout=grid_layout))
  if(left.width>0 & !is.null(left_item))
    print(left_item, vp=grid::viewport(layout.pos.col=1, layout.pos.row=2))
  if(top.height>0& !is.null(top_item))
    print(top_item, vp=grid::viewport(layout.pos.col=2, layout.pos.row=1))
  if(right.width>0& !is.null(right_item))
    print(right_item, vp=grid::viewport(layout.pos.col=3, layout.pos.row=2))
  if(bottom.height>0 & !is.null(bottom_item))
    print(bottom_item, vp=grid::viewport(layout.pos.col=2, layout.pos.row=3))
  ## print centre without legend
  print(center_item, vp=grid::viewport(layout.pos.col=2, layout.pos.row=2))

  if(show_legend){
    grid::pushViewport(grid::viewport(layout.pos.col=3, layout.pos.row=1))
    grid::grid.draw(legend)
  }

  grid::upViewport(0)
} 