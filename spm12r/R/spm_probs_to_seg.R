#' @title SPM Probability to Hard Segmenation
#'
#' @description Converts probability images from SPM 
#' segmentation to a hard, choose-one segmentation
#' @param img list of images for probabilities for each class
#' @param ties.method a character string specifying how ties are handled.  
#' See \code{\link{max.col}}.  Note, order of ties is different than 
#' \code{\link{max.col}}.
#' @export
#' @return Object of class nifti
#' @examples \dontrun{
#' spm_seg = spm12_segment(image)
#' seg = spm_probs_to_seg(spm_seg)
#'}
spm_probs_to_seg <- function(img,
    ties.method = c("first", "last", "random") 
    ){
    stopifnot(inherits(img, "list"))
    xmax = sapply(img, c)
    ties.method = match.arg(ties.method, c("first", "last", "random") )

    maxs = max.col(xmax, 
        ties.method = ties.method)
    res = niftiarr(img[[1]], maxs)   
    res
}