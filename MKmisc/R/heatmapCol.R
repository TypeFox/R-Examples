#########################################################
## function for computation of colors used for heatmaps
## data: matrix or data.frame; data which shall be displayed in a heatmap
##         ranging from negative to positive numbers
## col: vector of colors used for heatmap
## lim: constant colors are used for data below -lim resp. above lim
heatmapCol <- function(data, col, lim, na.rm = TRUE){
  nrcol <- length(col)
  data.range <- range(data, na.rm = na.rm)

  if(diff(data.range) == 0)
    stop("data has range 0")
  if(lim <= 0)
    stop("lim has to be positive")

  if(lim > min(abs(data.range))){
    warning("specified bound 'lim' is out of data range\n
                hence 'min(abs(range(data)))' is used")
    lim <- min(abs(data.range))
  }

  nrcol <- length(col)
  reps1 <- ceiling(nrcol*(-lim-data.range[1])/(2*lim))
  reps2 <- ceiling(nrcol*(data.range[2]-lim)/(2*lim))
  col1 <- c(rep(col[1], reps1), col, rep(col[nrcol], reps2))

  return(col1)
}
