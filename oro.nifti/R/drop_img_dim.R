#' @title Drop Image Dimension
#' @return Object of class nifti
#' @name dropImageDimension
#' @param img nifti object
#' @param onlylast is a logical variable (default = \code{TRUE}).  Drop the 
#' dimension only if it is the last dimension.  For example, if \code{dim} is 
#' 10x10x1x10 then no dimension is dropped, but if \code{dim} is 10x10x10x1 then 
#' it will be changed to 10x10x10.
#' @param warn produces a text output if the number of dimensions is under 
#' three.
#' @description Drops a dimension of an image that has one-dimension and 
#' sets respective values to 0 in \code{pixdim} or 1 in \code{dim}.
#' @importFrom abind adrop
#' @examples
#' 
#' nim <- nifti(array(rnorm(10^3), dim = rep(10, 3)))
#' nim2 <- nifti(array(rnorm(10^3), dim = c(10, 10, 1, 10)))
#' dropImageDimension(nim2)
#' dropImageDimension(nim2, onlylast = FALSE)
#' nim3 <- nifti(array(rnorm(10^3), dim = c(10, 10, 10, 1)))
#' dropImageDimension(nim3)
#' dropImageDimension(nim3, onlylast = FALSE) # the same as above
#' nim4 <- nifti(array(rnorm(10^3), dim = c(10, 10, 10, 1, 10)))
#' dim(nim4[,,,1,])
#' dim(nim4[,,,1,,drop=TRUE])
#' dropImageDimension(nim4)
#' 
#' nim5 <- nifti(array(rnorm(10^4), dim = c(1, 10, 10, 10, 1, 10)))
#' dropImageDimension(nim5)
#' dropImageDimension(nim5, onlylast = FALSE)
#' 
#' nim6 <- nifti(array(rnorm(10^3), dim = c(1, 10, 10, 10, 1, 1)))
#' dropImageDimension(nim6)
#' \dontrun{
#' ## 27 scans of Colin Holmes (MNI) brain co-registered and averaged
#' ## NIfTI two-file format
#' URL <- "http://imaging.mrc-cbu.cam.ac.uk/downloads/Colin/colin_1mm.tgz"
#' urlfile <- file.path(tempdir(), "colin_1mm.tgz")
#' download.file(URL, dest=urlfile, quiet=TRUE)
#' untar(urlfile, exdir=tempdir())
#' colin <- readNIfTI(file.path(tempdir(), "colin_1mm"))
#' dim(colin)
#' dim_(colin)
#' pixdim(colin)
#' # this will error
#' writeNIfTI(colin, filename = tempfile())
#' colin <- dropImageDimension(colin)
#' writeNIfTI(colin, filename = tempfile())
#' }
#' @rdname drop_img_dim
#' @export
dropImageDimension <- function(img, onlylast=TRUE, warn=TRUE) {
  dim_  <- dim_(img)
  imgdim <- dim(img)
  ####
  ndim <- length(imgdim) + 1
  dim_[seq(2, ndim)] <- imgdim
  if (ndim +1 <= length(dim_)) {
    dim_[seq(ndim+1, length(dim_))] <- 1
  }
  pdim <- pixdim(img)
  no.data <- dim_ <= 1
  no.data <- no.data | pdim == 0
  no.data[1] <- FALSE
  ## keeping cases like 10,10,1,10
  ## if onlylast is FALSE - drop anything that's a 1
  if (onlylast) {  
    maxdim <- max(which(! no.data))
    no.data[seq(maxdim)] <- FALSE
  } else {
    no.data[1] <- FALSE
  }
  ### subtract 1 for first observation
  ndim <- sum(! no.data) - 1
  dim_[1] <- ndim
  ### need the if statement in case 1x1x1 array (as is default)
  ### Must also if the dimensions are less than 3 then not an array
  pdim <- pdim[! no.data]
  pdim <- c(pdim, rep(1, 8 - length(pdim)))
  dim_ <- dim_[! no.data]
  dim_ <- c(dim_, rep(1, 8 - length(dim_)))
  #     dim_[no.data] = 1
  if (length(imgdim) > ndim) {
    pixdim(img) <- pdim
    dim_(img) <- dim_    
    if (onlylast) {
      ############# code for last only
      ## cs - so first must be a 1, then 2, for all TRUE, b/c reversed
      cs <- cumsum(rev(no.data[1 + seq(length(imgdim))]))
      ### if cs[1] = 1, and cs[2] = 2, then last cols
      dropcols <- cs == seq(length(imgdim))
      ### reverse it back to the correct order
      dropcols <- rev(dropcols)        
      dropcols <- which(dropcols)
      D <- adrop(img@.Data, drop = dropcols)        
    } else {
      D <- drop(img@.Data)
    }
  } else {
    return(img)
  }
  
  # dim_ must have values >= 1
  checkdim = dim_(img)
  checkdim[checkdim < 1] = 1
  dim_(img) <- checkdim
  
  if (ndim >= 3) {
    img@.Data <- D
    return(img)
    # img@.Data = drop(img@.Data)
  } else {
    if (warn) {
      warning("Dropping under 3 dimensions - returning non-nifti object array.")
    }
    return(D)
  }
}
#' @rdname drop_img_dim
#' @export
drop_img_dim <- function(img, onlylast=TRUE, warn=TRUE) {
  dropImageDimension(img=img, onlylast=onlylast, warn=warn)
}
