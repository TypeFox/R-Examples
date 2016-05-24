#' An Image of a Corrupted Face Image.
#' 
#' A 192 by 168 data matrix storing gray-scale pixel intensities of a corrupted face image. 
#' The original image was taken from Yale Face Database B [Georghiades et al. 2001] 
#' and then manually corrupted.
#' Run the following line to see the image: 
#' image(t(face[ncol(face):1, ]), axes = FALSE, col = grey(seq(0, 1, length = 256)), asp = 1)
#' 
#' @docType data
#' @keywords datasets
#' @format A data matrix containing 192 rows and 168 columns.
#' @name face
NULL


