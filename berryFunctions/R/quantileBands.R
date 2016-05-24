#' Quantile bands
#' 
#' Quantile bands with optional smoothing, e.g. for visualizing simulations
#' 
#' @return Quantiles of each column, invisible. Smoothed if \code{smooth} is given!
#' @note This is the first version and is not tested very well yet.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Sept 2014
#' @seealso \code{\link{quantile}}, \code{\link{quantileMean}}, \code{\link{ciBand}}, \code{\link{polygon}}, \url{cran.r-project.org/package=fanplot}
#' @keywords dplot aplot univar
#' @export
#' @examples
#' 
#' neff <- t(replicate(n=30, sapply(1:400, function(nn) max(rnorm(nn)))   ))
#' qB <- quantileBands(neff, x=1:400)
#' qB[,1:9]
#' quantileBands(neff, smooth=19, meanargs=list(col=2), txi=NA)
#' 
#' library(RColorBrewer)
#' 
#' quantileBands(neff, smooth=35, ylab="max of rnorm(n)",
#'   xlab="sample size (n)", probs=0:10/10, col=brewer.pal(5,"BuGn"),
#'   medargs=list(lwd=2), meanargs=list(col=2, lty=1), txi=c(40,50,60),
#'   main="Maximum is an unsaturated statistic:\n it rises with sample size")
#' 
#' neff2 <- t(replicate(n=50, sapply(1:400, function(nn) mean(rnorm(nn)))   ))
#' quantileBands(neff2, x=1:400, smooth=35, ylab="mean of rnorm(n)",
#'   xlab="sample size (n)", probs=0:10/10, col=brewer.pal(5,"BuGn"),
#'   txi=c(40,50,60), textargs=list(col="yellow"), medargs=list(lwd=2),
#'   meanargs=list(col=2, lty=1), main="Mean converges to true population mean")
#'    
#' @param mat Matrix or data.frame with columns of data
#' @param x X-axis positions for each column. DEFAULT: 1:ncol(mat)
#' @param col Vector of colors for each quantile group, recycled reversively if necessary. DEFAULT: rgb(0,0,1, alpha=c(0.5, 0.7))
#' @param add Add to existing plot? Allows to add to highly customized plot. DEFAULT: FALSE
#' @param main,xlab,ylab plot labels. DEFAULT: "Quantile Bands", ""
#' @param probs Probabilities passed to \code{\link{quantile}}. DEFAULT: 0:4/4
#' @param na.rm Remove NAs before computing \code{\link{quantile}s}, \code{\link{median}} and \code{\link{mean}}? DEFAULT: FALSE
#' @param type Which of the 9 \code{\link{quantile}} algorithms should be used. DEFAULT: 7
#' @param smooth If(!is.na), \code{width} passed to \code{\link{movAv}} smoothing quantiles. DEFAULT: NA
#' @param medargs List of arguments passed to lines drawing \code{\link{median}}. Not drawn if NULL. DEFAULT: NULL
#' @param meanargs List of arguments passed to lines drawing \code{\link{mean}}. Not drawn if NULL. DEFAULT: NULL
#' @param txi Text x position index (along columns of mat), recyled if necessary. NA to suppress. INTERNAL DEFAULT: middle of the plot for all.
#' @param textargs List of arguments passed to \code{\link{text}}, like col, adj, ... DEFAULT: NULL
#' @param \dots Further arguments passed to \code{\link{polygon}}, like border, lty, ... 
#' 
quantileBands <- function(
mat,
x=1:ncol(mat),
col=rgb(0,0,1, alpha=c(0.5, 0.7)),
add=FALSE,
main="Quantile Bands", ylab="", xlab="",
probs=0:4/4,
na.rm=FALSE,
type=7,
smooth=NA,
medargs=NULL,
meanargs=NULL,
txi,
textargs=NULL,
...)
{
# input check:------------------------------------------------------------------
if(length(x) != ncol(mat)) stop("length(x) (",length(x),") must equal ncol(mat) (",ncol(mat),").")
if(!is.vector(x)) stop("x must be a vector.")
if(!is.numeric(x)) stop("x must be numeric.")
# number of classes (and thus, polygon calls)
nclass <- length(probs)-1
# recycle color vectors (with reverse):
Ncol <- length(col)
if(Ncol > nclass) col <- col[1:nclass]
if(Ncol==1)
  col2 <- rep(col, nclass)
else
if(Ncol == nclass)
  col2 <- col
else
if(Ncol <= ceiling(nclass/2))
  {
  col2 <- rep(NA, nclass)
  col2[1:Ncol] <- col
  col2[nclass:(nclass-Ncol+1)] <- col # this is reverse
  col2[is.na(col2)] <- tail(col,1)
  }
else stop("col has wrong length (",Ncol,"). mus either be ==1, ==length(probs)-1 (==",
         nclass,") or <=length(probs)/2 (==", ceiling(nclass/2),").")
#
# Quantiles for each column
qc <- apply(mat, MARGIN=2, quantile, probs, na.rm, type)
# smooth:
if(!is.na(smooth))
  {
  qc <- t(apply(qc, MARGIN=1, movAv, smooth))
  nna <- !is.na(qc[1,])
  qc <- qc[,nna]
  x <- x[nna]
  }
else nna <- 1:ncol(qc)
# plot initiation:
if(!add) plot(x, qc[1,], type="n", ylim=range(qc, na.rm=TRUE),
              main=main, ylab=ylab, xlab=xlab)
# actual polygon drawing:
for(i in 1:nclass)
   polygon(x=c(x, rev(x)), y=c(qc[i,], rev(qc[i+1,])), col=col2[i], ...)
# add median line:
if(!is.null(medargs))
   {
   medians <- apply(mat, MARGIN=2, median, na.rm=na.rm)[nna]
   do.call(lines, owa(list(x=x, y=medians), medargs))
   }
# add mean line:
if(!is.null(meanargs))
   {
   means <- apply(mat, MARGIN=2, mean, na.rm=na.rm)[nna]
   do.call(lines, owa(list(x=x, y=means), meanargs))
   }
# text
if(missing(txi)) txi <- round(ncol(qc)/2)
if(any(!is.na(txi)))
  {
  txi <- rep(txi, length=nrow(qc))
  utxi <- na.omit(unique(txi))
  if(any(utxi < 1)) stop("txi (",utxi,") must be a vector of positive integers.")
  if(any(utxi > ncol(qc))) stop("txi (",utxi,") must be smaller than ncol(mat)-smooth (",ncol(qc),").")
  do.call(text, args=owa(list(x=x[txi], y=diag(qc[,txi]), labels=rownames(qc)), textargs))
  }
# output
colnames(qc) <- x
return(invisible(qc))
} # Function end
