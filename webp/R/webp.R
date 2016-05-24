#' Webp image format
#'
#' Read and write webp images into a bitmap array. The bitmap array uses the same
#' conventions as the \code{png} and \code{jpeg} package.
#'
#' @export
#' @useDynLib webp R_webp_decode
#' @rdname read_webp
#' @param source raw vector or path to webp file
#' @param numeric convert the image to 0-1 real numbers to be compatible with
#' images from the jpeg or png package.
#' @examples # Convert to webp
#' library(png)
#' img <- readPNG(system.file("img", "Rlogo.png", package="png"))
#' write_webp(img, "rlogo.webp")
#' browseURL("rlogo.webp")
#' rm(img)
#'
#' # Convert from webp
#' library(jpeg)
#' img <- read_webp("rlogo.webp")
#' writeJPEG(img, "rlogo.jpeg")
#' browseURL("rlogo.jpeg")
read_webp <- function(source, numeric = TRUE) {
  if(is.character(source))
    source <- readBin(source[1], raw(), file.info(source)$size)
  stopifnot(is.raw(source))
  out <- .Call(R_webp_decode, source)
  if(isTRUE(numeric)){
    out <- structure(as.numeric(out)/255, dim = dim(out))
    out <- aperm(out)
  } else {
    class(out) <- c("rawimg", class(out))
  }
  out
}

#' @export
#' @rdname read_webp
#' @useDynLib webp R_webp_encode
#' @param image array of 3 dimensions (width * height * channel) with real numbers
#' between 0 and 1.
#' @param target path to a file or \code{NULL} to return the image as a raw vector
#' @param quality value between 0 and 100
write_webp <- function(image, target = NULL, quality = 80) {
  if(is.numeric(image)){
    image <- structure(as.raw(image * 255), dim = dim(image))
    image <- aperm(image)
  }
  channels = dim(image)[1]
  stopifnot(channels == 3 || channels == 4)
  buf <- .Call(R_webp_encode, image, quality)
  if(is.character(target))
    writeBin(buf, target)
  else
    structure(buf, class = "webp")
}

#' @useDynLib webp R_webp_get_info
webp_dims <- function(buf) {
  stopifnot(is.raw(buf))
  .Call(R_webp_get_info, buf)
}

#' @export
print.rawimg <- function(x, ...){
  dims <- dim(x)
  cat(sprintf("raw image (%d x %d) with %d channels\n", dims[2], dims[3], dims[1]))
  invisible()
}

#' @export
print.webp <- function(x, ...){
  dims <- webp_dims(x)
  cat(sprintf("webp buffer (%d x %d)\n", dims[1], dims[2]))
  invisible()
}
