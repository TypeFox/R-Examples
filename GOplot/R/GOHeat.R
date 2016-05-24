#' 
#' @name GOHeat
#' @title Displays heatmap of the relationship between genes and terms.
#' @description The GOHeat function generates a heatmap of the relationship 
#'   between genes and terms. Biological processes are displayed in rows and
#'   genes in columns. In addition genes are clustered to highlight groups of
#'   genes with similar annotated functions. The input can be generated with the
#'   \code{\link{chord_dat}} function.
#' @param data The matrix represents the binary relation (1= is related to, 0= 
#'   is not related to) between a set of genes (rows) and processes (columns)
#' @param nlfc Defines the number of logFC columns (default = 0)
#' @param fill.col Defines the color scale break points
#' @details The heatmap has in general two modes which depend on the \code{nlfc}
#'   argument. If \code{nlfc = 0}, so no logFC values are available, the 
#'   coloring encodes for the overall number of processes the respective gene is
#'   assigned to. In case of \code{nlfc = 1} the color corresponds to the logFC 
#'   of the gene.
#' @import ggplot2
#' @examples
#' \dontrun{
#' # Load the included dataset
#' data(EC)
#' 
#' # Generate the circ object
#' circ <- circle_dat(EC$david, EC$genelist)
#' 
#' # Generate the chord object
#' chord <- chord_dat(circ, EC$genes, EC$process)
#' 
#' # Create the plot with user-defined colors
#' GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))
#' }
#' @export
#' 

GOHeat <- function(data, nlfc, fill.col){
  x <- y <- z <- NULL
  if(missing(nlfc)) nlfc <- 0 else nlfc <- nlfc
  if(missing(fill.col)) fill.col <- c('firebrick', 'white', 'dodgerblue') else fill.col <- fill.col
  
  distance <- dist(data)
  cluster <- hclust(distance)
  M <- dim(data)[2]
  nterm <- M - nlfc
  if(nlfc == 0){
    s <- rowSums(data[,1:nterm])
    tmp <- NULL
    for(r in 1:nrow(data)){
      tmp <- c(tmp, as.numeric(gsub(1, s[r], data[r, 1:nterm])))
    }
  }else{
    tmp <- NULL
    for(r in 1:nrow(data)){
      tmp <- c(tmp, as.numeric(gsub(1, data[r, (nterm + 1)], data[r, 1:nterm])))
    }
  }
  df <- data.frame(x = rep(cluster$order, each = nterm), y = rep(colnames(data[,1:nterm]), length(rownames(data))), z = tmp, 
                   lab = rep(rownames(data), each = nterm))
  df_o <- df[order(df$x),]
  
  g <- ggplot() + 
        geom_tile(data = df_o, aes(x = x, y = y, fill = z))+
        scale_x_discrete(breaks = 1:length(unique(df_o$x)), labels = unique(df_o$lab)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.title.x=element_blank(), axis.title.y=element_blank(),
              axis.text.y = element_text(size = 14), panel.background=element_blank(), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank())
  if(nlfc == 0){
    g + scale_fill_gradient2('Count', space = 'Lab', low=fill.col[2], mid=fill.col[3], high=fill.col[1])
  }else{
    g + scale_fill_gradient2('logFC', space = 'Lab', low=fill.col[3], mid=fill.col[2], high=fill.col[1])
  }
}



