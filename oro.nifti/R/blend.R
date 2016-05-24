setGeneric("blend", function(x, y, ...) standardGeneric("blend"))
#' Merge Two NIfTI or ANALYZE Volumes
#' 
#' @description Two volumes of medical imaging data are merged together in the 
#' superior-inferior (or $z$) direction.  One assumes that there is at least one
#' slice that overlaps between the two volumes.   
#'
#' @param x,y are objects of class \code{nifti} or \code{anlz}.
#' @param seqX,seqY are vectors that provide the $z$-coordinate values for the 
#' two imaging volumes.
#' @param method is the type of weighing to use when combining information where 
#' there is an overlap (default = \code{"linear"}).  
#' @return A single volume that blends the voxel-wise information from \code{x} 
#' and \code{y}.
#' @section Methods: 
#' \describe{ 
#' \item{x = "nifti", y = "nifti"}{Merge \code{x} and \code{y}.} 
#' \item{x = "anlz", y = "anlz"}{Merge \code{x} on \code{y}.} 
#' \item{x = "nifti", y = "anlz"}{Merge \code{x} on \code{y}.} 
#' \item{x = "anlz", y = "nifti"}{Merge \code{x} and \code{y}.} 
#' }
#' @name blend
#' @rdname blend-methods
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{image-methods}}, \code{\link{overlay-methods}}
#' @keywords methods
#' @docType methods
#' @export
blendVolumes <- function(x, y, seqX, seqY, method = "linear") {
  ## Test for overlap
  if (min(seqX) > max(seqY) || max(seqX) < min(seqY)) {
    stop("There is no overlap between the two volumes.")
  }
  ## Test dimensions
  nr <- nrow(x) # must be the same as y
  nc <- ncol(x) # must be the same as y
  if (nr != nrow(y) || nc != ncol(y)) {
    stop("Number of rows/columns for the two volumes must be equal.")
  }
  ## Test order of input volumes
  if (min(seqY) > min(seqX) && max(seqY) > max(seqX)) {
    ## Swap definitions of x and y
    cat("# SWAP!", fill = TRUE)
    tmp <- x
    x <- y
    y <- tmp
    tmp <- seqX
    seqX <- seqY
    seqY <- tmp
    rm(tmp)
  }
  ns <- list(x = nsli(x), y = nsli(y))
  dist <- list (XY = seqX - max(seqY), YX = seqY - min(seqX))
  overlap <- list(start = which(abs(dist$XY) == min(abs(dist$XY))),
                  end = which(abs(dist$YX) == min(abs(dist$YX))))
  blendRange <- range(1:ns$x, overlap$start + 0:(ns$y - 1))
  ##
  index <- list(top = 1:(overlap$start - 1), # top?
                middle = overlap$start:(overlap$start + overlap$end - 1),
                bottom = (overlap$start + overlap$end):max(blendRange)) # bottom?
  index <- list(top = max(blendRange) - (overlap$start - 2):0,
                middle = max(blendRange) - (overlap$start + overlap$end - 2):(overlap$start - 1),
                bottom = max(blendRange) - max(blendRange - 1):(overlap$start + overlap$end - 1))
  ##
  blendVolume <- array(0, c(nr, nc, max(blendRange)))
  blendVolume[, , index$top] <- x@.Data[, , ns$x - (overlap$start - 2):0] # x@.Data[, , 1:(overlap$start - 1)]
  blendVolume[, , index$bottom] <- y@.Data[, , ns$y - (ns$y - 1):overlap$end] # y@.Data[, , (overlap$end + 1):ns$y]
  if (method == "linear") {
    mix <- array(rep(0:(overlap$end-1) / overlap$end, each = nr * nc),
                 dim = c(nrow(y), ncol(y), overlap$end))
  } else {
    stop("Only method = \"linear\" has been implemented.")
  }
  blendVolume[, , index$middle] <-
    mix * x@.Data[, , ns$x - (ns$x - 1):(overlap$start - 1)] + (1-mix) * y@.Data[, , ns$y - (overlap$end - 1):0]
  blendSequence <- c(seqX[1:(overlap$start-1)], 
                     rowMeans(cbind(seqX[overlap$start:ns$x], seqY[1:overlap$end])), 
                     seqY[(overlap$end+1):ns$y])
   if (is.nifti(x)) {
     as(blendVolume, "nifti") <- x
   } else {
     if (is.nifti(y)) {
       as(blendVolume, "nifti") <- y
     } else {
       as(blendVolume, "anlz") <- x
     }
   }
  list(img = blendVolume, seq = blendSequence)
}
#' @export
#' @rdname blend-methods
#' @aliases blend,nifti,nifti-methods
setMethod("blend", signature(x="nifti", y="nifti"), blendVolumes)
#' @export
#' @rdname blend-methods
#' @aliases blend,anlz,anlz-methods
setMethod("blend", signature(x="anlz", y="anlz"), blendVolumes)
#' @export
#' @rdname blend-methods
#' @aliases blend,anlz,nifti-methods
setMethod("blend", signature(x="anlz", y="nifti"), blendVolumes)
#' @export
#' @rdname blend-methods
#' @aliases blend,nifti,anlz-methods
setMethod("blend", signature(x="nifti", y="anlz"), blendVolumes)
