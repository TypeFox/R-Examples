#' readbitmap: read standard bitmap image formats
#' 
#' The readbitmap package enables users to read the three main general purpose 
#' bitmap image formats (jpeg, png, bmp) without having to specify the image 
#' type directly. This is provided by a single function 
#' \code{\link{read.bitmap}}, which uses a second function 
#' \code{\link{image_type}}, which is also exported for users, to identify the 
#' image type of a file using appropriate \emph{magic} values encoded in the 
#' first few bytes of the image header. Images can therefore have any file
#' extension.
#' 
#' @name readbitmap-package
#' @aliases readbitmap
#' @docType package
#' @keywords package
#' @seealso \code{\link{image_type}}, \code{\link{read.bitmap}}
NULL
