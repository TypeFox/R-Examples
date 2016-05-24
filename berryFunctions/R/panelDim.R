#' Arrange panels in a multipanel plot (par mfrow)
#' 
#' Returns the optimum where deviation from ncol=nrow and number of panels left empty have a minimum sum.
#' 
#' @details There probably are other ways to find the optimal way to arrange panels, so if you find anything, please give me a hint.
#' 
#' @return vector with 2 values, can be passed to par(mfrow), or a data.frame if \code{all=TRUE}.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Dec 2014, Jan 2015
#' @seealso \code{\link{groupHist}}, which is using this function
#' @keywords dplot array optimize
#' @export
#' @examples
#' 
#' # basic usage
#' op <- par(mfrow=panelDim(6))
#' for(i in 1:6) plot(i:10, main=i)
#' par(op)
#' 
#' # Adcanced options
#' panelDim(7)
#' g <- panelDim(7, all=TRUE)
#' panelDim(7, plot=TRUE)
#' panelDim(7, plot=TRUE, all=TRUE) # compares 4 best options
#' 
#' panelDim(26, all=TRUE)
#' panelDim(26, plot=TRUE, all=TRUE) # compares 4 best options
#' panelDim(26, plot=TRUE, all=TRUE, weight=c(3,0) ) # fewer empty panels
#' 
#' # effect of maxempty:
#' panelDim(13, plot=TRUE)             #  4 x 4
#' panelDim(13, maxempty=2, plot=TRUE) #  5 x 3
#' panelDim(13, maxempty=1, plot=TRUE) #  7 x 2
#' panelDim(13, maxempty=0, plot=TRUE) # 13 x 1
#' 
#' panelDim(45, plot=TRUE) # no empty panels
#' # focus on aspect ratio of each panel (make it as square as possible):
#' panelDim(45, weight=c(1,3), plot=TRUE) # better aspect for each panel
#' 
#' # Orientation of plot:
#' panelDim(45, plot=TRUE) # good for portrait orientation of plot
#' panelDim(45, landscape=TRUE, plot=TRUE) # better if plot width > height
#' 
#' \dontrun{
#' ## Rcmd check --as-cran doesn't like to open external devices,
#' ## so this example is excluded from running in the checks.
#' plot of several n with defaults
#' dev.new(record=TRUE)
#' for(i in 1:50)  panelDim(i, plot=TRUE)
#' }
#' 
#' @param n Number of panels to be arranged
#' @param weight Weights to avoid \emph{empty panels} and \emph{discrepancy between ncol and nrow}, respectively. DEFAULT: c(1,1)
#' @param maxempty Maximum number of panels that are allowed to be left empty.
#'        If \code{maxempty=0}, no panel is left blank, so 11 plots would be beneath each other 
#'        instead of in a 4x3 grid with one panel left blank. DEFAULT: round(n/4)
#' @param landscape Use landscape orientation instead of portrait? DEFAULT: FALSE
#' @param all Show all reasonable possibilities in a data.frame? DEFAULT: FALSE
#' @param plot Show the panel layout result? (the 4 best options are compared if \code{all=TRUE}). DEFAULT: FALSE
#' @param mfcol use mfcol instead of mfrow. DEFAULT: FALSE
#' 
panelDim <- function(
n,
weight=c(1,1),
maxempty=round(n/4),
landscape=FALSE,
all=FALSE,
plot=FALSE,
mfcol=FALSE
)
{
# Input control:
n <- as.integer(n[1])
if(n<1) stop("n must be an integer larger than 1")
weight <- rep(weight, length=2)
maxempty <- maxempty[1]
#
# matrix of possibilities:
n2 <- n+maxempty 
m <- outer(1:n2, 1:n2)
# Get the possibilities that would yield n (or up to n+maxempty) panels
getpossibilities <- function(nblank)
  {
  x <- y <- NA;   i <- 0
  while(is.na(x))
    {
    d <- if(i==0) diag(m) else diag(m[, -(1:i)])
    if(any(d==n+nblank))
      {
      y <- which(d==n+nblank)
      x <- y+i
      }   else  x <- y <- NA
    i <- i+1
    }
  c(x, y) # function output
  }
# Now actually get them:
g <- as.data.frame(t(sapply(0:maxempty, getpossibilities)))
colnames(g) <- c("x","y")
# Number of empty panels:
g$empty <- g$x * g$y - n
# difference between ncol and nrow, = distance from diagonal:
g$diff <- g$x - g$y
# ranking of possibilities via weigted penalty points: 
penalty <- weight[1]*g$empty + weight[2]*g$diff
g <- g[order(penalty),]
rownames(g) <- NULL
# use the combination with a minimum of penalty points:
xy <- as.numeric( g[1, c("x","y")] )
if(landscape) xy <- rev(xy)
#
# Show resulting panel layout:
if(plot & !all)
  {
  op <- par(mfrow=xy, mar=rep(0,4))
  on.exit( par(op) )
  if(mfcol) par(mfcol=xy)
  for(i in 1:n)
    {
    plot(1, type="n", ann=F, axes=F)
    text(1,1, i, cex=2)
    box("figure")
    }
  }
if(plot & all) 
  {
  op <- par(mfrow=c(2,2), mar=rep(0.3,4))
  on.exit( par(op) )
  # 4 best plots (b):
  for(b in 1:pmin(4, nrow(g)))
    {
    nr <- g[b,1]  ; nc <- g[b,2] 
    if(landscape) { nr <- g[b,2]  ; nc <- g[b,1]  }
    #cat("b: ", b, ", nr: ", nr, ", nc: ", nc, "\n")
    plot(1, ylim=c(nr,0), xlim=c(0,nc), type="n", yaxs="i", xaxs="i", ann=F, axes=F)
    #if(mfcol) 
    for(j in 1:n) rect(xleft=(j-1)%%nc, ybottom=ceiling(j/nc), 
                      xright=(j-1)%%nc+1,  ytop=ceiling(j/nc)-1, col=8)
    #box("figure")
    }
  }
# Function output
if(all) return(g) else return(xy)
} # End of function
