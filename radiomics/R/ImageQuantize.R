#' Image Discretization.
#'
#' \code{discretizeImage} Scales the grey values of an image into a specified number of values.
#'
#' This function is called in \code{glcm}, \code{glrlm}, \code{glszm}, and \code{mglszm}.
#' 
#' If n_grey is greater than the number of unique grey levels in the matrix then no action
#' is taken.
#' 
#' @param data A numeric 2D matrix.
#' @param n_grey an integer value, the number of grey levels the image should
#'   be quantized into.
#' @param verbose Logical, a message is given when the user 
#'  supplies more grey values than exist in the image. Setting this value to FALSE
#'  will suppress this message.       
#' @return A matrix of the same dimensions as the input matrix. The entries of the matrix
#'  will be set to begin at 1, and go up to the specified value. There is no guarantee
#'  that each gray level between 1 and n_grey will have pixels of that value (for example, 
#'  although n_grey = 32 may be specified, certain images may contain fewer than 32 grey levels).
#'     
#'
#' @examples
#' image(psf)
#' image(discretizeImage(psf, n_grey=5, verbose=F))
#' 
#' image(tumor)
#' image(discretizeImage(tumor, n_grey=8, verbose=F))
#' @export
 
discretizeImage <- function(data, n_grey=32, verbose=TRUE){
  #Error checking
  #Check validity of input
  if (!is.matrix(data)) {
    stop(paste0("Object of class ", class(data), ".  is.matrix(object) must evaluate TRUE."))
  }
  if (any(data < 0, na.rm=TRUE)) {
    stop("Object contains negative values. All values must be greater than 0.")
  }
  
  #Not a perfect solution. Makes n_grey breaks, but doesn't necessarily populate all of them
  # eg. n_gey could be 100, but only 75 of the levels are used by pixels
  l_unique <- length(unique(c(data)))
  
  #Make sure discretization is valid
  if(l_unique == 0 ) stop("Function not valid for empty input")
  if(sum(is.na(data)) == length(data)){
    if(verbose) warning("Matrix composed entirely of NA's\nReturning Matrix")
    return(data)
  } 
  
  #If we don't need to do anything, we don't do anything
  if(n_grey == l_unique){
    return(data)
  } else if(n_grey > l_unique){
    if(verbose) message(sprintf("n_grey (%d) cannot be larger than the number of gray levels in the image (%d)", n_grey, l_unique))
    if(verbose) message(sprintf("n_grey set to %d", l_unique))
    return(data)
  } 
  
  discretized <- cut(data, breaks=seq(min(data, na.rm=TRUE), max(data, na.rm=TRUE), length.out=(n_grey + 1)),
                     labels = seq(1, n_grey, 1),
                     include.lowest=TRUE, right=FALSE) 
  return(matrix(as.numeric(discretized), nrow=nrow(data)))
}

#' Image Discretization.
#' 
#'#' \code{discretizeImage2} Scales the grey values of an image into a specified number of values.
#'
#' Not currently used. Different methods of discretizing the image will be explored in the future.
#'
#'@param image A numeric image matrix. 
#'@param n_grey The grey levels the output image will have
discretizeImage2 <- function(image, n_grey=32){
  #Different from discretizeImage2 in that it preserves the grey levels of the image (roughly)
  l_unique <- length(unique(c(image)))
  if(identical(n_grey, l_unique)){
    return(image)
  } else if(n_grey > l_unique){
    message(sprintf("n_grey (%d) cannot be larger than the number of gray levels in the image (%d)", n_grey, l_unique))
    message(sprintf("n_grey set to %d", l_unique))
    return(image)
  } 
  
  
  break_vals <- seq(min(image, na.rm=TRUE), max(image, na.rm=TRUE), length.out=(n_grey + 1))
  label_vals <- round(break_vals[-1], 4)
  #need to remove the upper bound label, as we are using the floor of the values
  #label_vals <- label_vals[-length(label_vals)] 
  
  discretized <- as.numeric(as.character(cut(image, breaks=break_vals,
                                             labels = label_vals,
                                             include.lowest=TRUE, right=FALSE))) 
  return(matrix(discretized, nrow=nrow(image)))
}
