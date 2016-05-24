##' a square plot, ala old NY Times
##'
##' Used as an alternative to a segmented barplot when the actual count is of interest.
##' @param x data
##' @param col color
##' @param border show border
##' @param nrows number of rows
##' @param ncols number of columns
##' @param ... passed on
##' @return NULL
##'
##' @export
squareplot <- function(x,
                       col = gray(seq(.5,1,length=length(x))),
                       border=NULL,
                       nrows=ceiling(sqrt(sum(x))),
                       ncols=ceiling(sum(x)/nrows),
                       ...
                       ) {
  ## create a squareplot ala the New York Times. Used as an
  ## alternative to a segmented barplot when the actual 
  ## count is of interest.

  ## helper function
  draw.square <- function(x,y,w=1,...) {
    ## draw a square with lower left corner at (x,y)
    polygon(x+c(0,0,w,w,0),y+c(0,w,w,0,0),...)
  }

  ## size of big square
  square.size = max(nrows,ncols)

  ## setup window with plot.new() and plot.window()
  ## arguments to ... are passed along here
  plot.new()
  plot.window(xlim=c(0,square.size),ylim=c(-square.size,0))
  title(...)
  
  ## vector with colors
  cols = rep(col,x)

  for(i in 1:sum(x)) {
    x.pos = floor((i-1)/nrows)          # adjust by 1
    y.pos = (i-1) %% nrows

    draw.square(x.pos,-y.pos -1,col=cols[i])
  }
}
