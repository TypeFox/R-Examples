#' @title OASIS Prediction
#' @description This function creates the OASIS probability map from a single MRI study with FLAIR, T1, T2, and PD volumes. 
#' @param flair flair volume of class \code{\link{nifti}}
#' @param t1 t1 volume of class \code{\link{nifti}}
#' @param t2 t2 volume of class \code{\link{nifti}}
#' @param pd pd volume of class \code{\link{nifti}}
#' @param brain_mask brain mask of class \code{\link{nifti}}, 
#' if \code{NULL} a brain mask will be created using \code{\link{fslbet}}.  
#' Note that provided brain masks should be in the same space as the T1 volume 
#' if \code{preproc = TRUE}, as all volumes will be registered to this space 
#' @param preproc is a logical value that determines whether to 
#' call \code{\link{oasis_preproc}} and performs the necessary preprocessing steps 
#' for OASIS
#' @param normalize is a logical value that determines whether to 
#' perform z-score normalization using \code{\link{zscore_img}} 
#' of the image over the brain mask, 
#' should be \code{TRUE} unless you train model
#' using an alternative normalization 
#' @param model an object of class \code{\link{glm}} used to make the OASIS predictions 
#' @param return_preproc is a logical value that indicates whether the 
#' preprocessed images should be returned, if \code{NULL} 
#' then the model from the OASIS paper will be used 
#' @param binary logical indicating whether a binary map 
#' should be returned by thresholding the probability map
#' @param threshold numeric indicating the threshold value 
#' for the probability map, with default of 0.16 for the OASIS paper
#' @param cores numeric indicating the number of cores to be 
#' used (no more than 4 is useful for this software implementation)
#' @import fslr
#' @import parallel
#' @return If \code{return_preproc = FALSE} the function returns a 
#' volume of class \code{\link{nifti}} containing the OASIS probability for each voxel. 
#' Otherwise, the function returns a list of volumes: 
#' the OASIS probability map, the FLAIR volume, the T1 volume, the T2 volume,
#' the PD volume, the brain mask for the subject, and the voxel selection mask. 
#' @examples \dontrun{
#' library(fslr)
#' flair <- readnii('path/to/flair', reorient = FALSE) 
#' t2 <- readnii('path/to/t2', reorient = FALSE) 
#' t1 <- readnii('path/to/t1', reorient = FALSE) 
#' pd <- readnii('path/to/pd', reorient = FALSE) 
#' oasis_map <- oasis_predict(flair = flair, t1 = t1, t2 = t2, pd = pd) }
#' @export 
oasis_predict <- function(flair, ##flair volume of class nifti
                          t1, ##t1 volume of class nifti
                          t2, ##t2 volume of class nifti
                          pd, ##pd volume of class nifti
                          brain_mask = NULL, ##brain mask of class nifti
                          preproc = FALSE, ##option to preprocess the data
                          normalize = TRUE, ##option to normalize 
                          model = NULL, ##an OASIS model of class glm
                          return_preproc = FALSE, ##option to return the preprocessed data
                          binary = FALSE,
                          threshold = 0.16, 
                          cores = 1
) {
  flair = check_nifti(flair)
  t1 = check_nifti(t1)
  t2 = check_nifti(t2)
  pd = check_nifti(pd)
  
  ##correct image dimmension
  flair <- correct_image_dim(flair)
  t1 <- correct_image_dim(t1)
  t2 <- correct_image_dim(t2)
  pd <- correct_image_dim(pd)
  
  
  ##image preproceesing 
  if (preproc == TRUE) {
    ## the image preproceesing 
    preprocess <- oasis_preproc(flair = flair, t1 = t1, t2 = t2, pd = pd, 
                                cores = cores,
                                brain_mask = brain_mask)
    oasis_study <- preprocess[c("flair","t1", "t2", "pd")]
    brain_mask <- preprocess$brain_mask
  } else{
    ## no preprocessing  
    oasis_study <- list(flair = flair, t1 = t1, t2 = t2, pd = pd)
  }
  if (is.null(brain_mask) & !preproc){
    ## create a brain mask if not supplied
    brain_mask <- fslbet(infile = oasis_study$t1, retimg = TRUE)
  } 
  brain_mask = check_nifti(brain_mask)
  brain_mask <- brain_mask > 0
  brain_mask <- datatyper(brain_mask, trybyte = TRUE)  
  
  
  ##adjust brain mask for OASIS 
  brain_mask <- correct_image_dim(brain_mask)
  brain_mask <- fslerode(brain_mask, kopts = "-kernel box 5x5x5", retimg = TRUE)
  cutpoint <- quantile(oasis_study$flair[brain_mask == 1], probs = .15)
  brain_mask[oasis_study$flair <= cutpoint] <- 0 
  
  ## the image normalization 
  if (normalize == TRUE) {
    oasis_study <- lapply(oasis_study, zscore_img, 
                          mask = brain_mask, 
                          margin = NULL)
  }
  
  orig_study = oasis_study
  ## smooth the images using fslsmooth from the fslr package 
  smooth <- mclapply(orig_study, fslsmooth,
                     sigma = 10, 
                     mask = brain_mask, 
                     retimg = TRUE, 
                     smooth_mask = TRUE,
                     mc.cores = cores)
  names(smooth) = paste0(names(smooth), "_10")
  oasis_study = c(oasis_study, smooth)
  
  
  smooth <- mclapply(orig_study, fslsmooth,
                     sigma = 20, 
                     mask = brain_mask, 
                     retimg = TRUE, 
                     smooth_mask = TRUE,
                     mc.cores = cores)
  names(smooth) = paste0(names(smooth), "_20")
  oasis_study = c(oasis_study, smooth)  
  
  rm(list = c("orig_study", "smooth"))
  
  ##create and apply the voxel selection mask 
  top_voxels <- voxel_selection(flair = oasis_study$flair,
                                brain_mask = brain_mask, 
                                cutoff = .85)
  
  ## create a dataframe to make oasis predictions on    
  oasis_study <- lapply(oasis_study, function(x) x[top_voxels == 1])
  oasis_dataframe <- do.call(cbind.data.frame, oasis_study)
  
  names <- c("FLAIR", "T1", "T2", "PD")
  colnames(oasis_dataframe) <-  c(names, 
                                  paste0(names, "_10"),  
                                  paste0(names, "_20"))
  
  ## make the model predictions 
  if (is.null(model)) {
    predictions <- predict( oasis::oasis_model, 
                            newdata = oasis_dataframe, 
                            type = 'response')    
  } else { 
    predictions <- predict(model, 
                           newdata = oasis_dataframe, 
                           type = 'response')    
  }
  
  ##put the predictions onto the brain 
  predictions_nifti <- niftiarr(brain_mask, 0) 
  predictions_nifti[top_voxels == 1] <- predictions
  
  
  ##smooth the probability map 
  prob_map <- fslsmooth(predictions_nifti, sigma = 1.25, 
                        mask = brain_mask, retimg = TRUE,
                        smooth_mask = TRUE)
  
  if (binary == TRUE) {
    binary_map <- prob_map
    binary_map[prob_map > threshold] <- 1
    binary_map[prob_map <= threshold] <- 0
  }
  if (return_preproc == TRUE & binary == FALSE) {
    return(list(oasis_map = prob_map, flair = flair, 
                t1 = t1, t2 = t2,
                pd = pd, brain_mask = brain_mask, 
                voxel_selection = top_voxels,
                unsmoothed_map = predictions_nifti))
  } 
  if (return_preproc == TRUE & binary == TRUE) {
    return(list(oasis_map = prob_map, flair = flair, 
                t1 = t1, t2 = t2,
                pd = pd, brain_mask = brain_mask, 
                voxel_selection = top_voxels, 
                binary_map = binary_map,
                unsmoothed_map = predictions_nifti))
  }   
  if (return_preproc == FALSE & binary == TRUE) {
    return(list(oasis_map = prob_map, 
                binary_map = binary_map,
                unsmoothed_map = predictions_nifti))
  } else{
    return(prob_map)
  }
  
}


