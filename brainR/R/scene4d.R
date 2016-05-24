#' Wrapper to write a 4D scene
#'
#' This function takes in filenames, levels, and creates an output html file,
#' with 4D elements.  The html is based on XTK (https://github.com/xtk/X#readme)
#' 
#' @param files (character) vector of filenames (first being a brain file if useTemp=FALSE)
#' @param fnames (character) filenames for the 3D surfaces in the scene - needs to 
#' be the same length as files
#' @param outfile (character) html filename 
#' @param levels (numeric/list) levels to make contours/surfaces for each file. 
#' Either a numeric vector may be passed, one level for each file.  Or a list of numeric vectors of multiple
#' levels for each file. Will be coerced to a list.
#' @param alpha (numeric/list) alpha opacities for each contours/surface for each file. Will be coerced to list
#' similarly as levels
#' @param color (character/list) colors for each contours/surface for each file. Will be coerced to list
#' similarly as levels
#' @param useTemp (logical) whether to use template from brainR as the brain figure
#' @param MNITemp (character) if (useTemp = TRUE) either "1mm" or "2mm" denoting the resolution 
#' of the template used
#' @param objtype (character) object type to write the files to.  Either "stl", "obj", or "ply" to write.
#' @param ... other options to be passed to \link{write4D}
#' @export
#' @examples
#' ### Faster - 8mm resampled but very coarse
#' imgs <- paste("Visit_", 1:5, "_8mm.nii.gz", sep="") 
#'  files <- sapply(imgs, system.file, package='brainR')
#' scene4d(files, levels=rep(0.99, length(files)), color= rep("blue", length(files)), useTemp=TRUE, 
#' MNITemp = "8mm", alpha = rep(1, length(files)), rescale=TRUE  )
#' \dontrun{
#' imgs <- paste("Visit_", 1:5, ".nii.gz", sep="") 
#'  files <- sapply(imgs, system.file, package='brainR')
#' scene4d(files, levels=rep(0.99, length(files)), color= rep("blue", length(files)), useTemp=TRUE, 
#' MNITemp = "8mm", alpha = rep(1, length(files)), rescale=TRUE  )
#' }
#' @return NULL


scene4d <- function(files, fnames=NULL, outfile = "index_4D_stl.html", levels=NULL, alpha=NULL, color="white", 
                    useTemp=FALSE, MNITemp= c("1mm", "2mm"), objtype = "stl", ...){
  
  ### checking object type
  objtype <- tolower(objtype)
  objtype <- gsub(".", "", objtype, fixed=TRUE)
  stopifnot(objtype %in% c("stl", "obj", "ply"))
  
    
  ## make list of these
  levels <- as.list(levels)
  alpha <- as.list(alpha)
  color <- as.list(color)
  
  ## make output image names from image names
  if (useTemp) {
    mnifile <- system.file(paste0("MNI152_T1_", MNITemp, "_brain.nii.gz"), package="brainR")
    files <- c(mnifile, files)
    levels <- c(4500, levels)
    alpha <- c(0.8, alpha)
    color <- c("white", color)
  }
  
  ## makes fnames from the filenames - checking for nii.gz, .nii, .img, .hdr
  if (is.null(fnames)){
    fnames <- basename(files)
    nii.gz <- grepl("\\.nii\\.gz$", fnames)
    nii <- grepl("\\.nii$", fnames)
    img <- grepl("\\.img$", fnames)
    hdr <- grepl("\\.hdr$", fnames)
    types <- cbind(nii.gz, nii, img, hdr)
    checktypes <- apply(types, 1, any)
    if (!all(checktypes)) stop("Not type in nii.gz, .nii, .img, .hdr")
    types <- c("\\.nii\\.gz$", "\\.nii$", "\\.img$", "\\.hdr$")[apply(types, 1, which)]
    for (ifile in seq_along(fnames)){
      fnames[ifile] <- gsub(paste0("(.*)", types[ifile]), "\\1", fnames[ifile])
    }
  }


  
  ## first file is always assumed as the template
  template <- readNIfTI(files[1], reorient=FALSE) 
  dtemp <- dim(template)
 ### 4500 - value that empirically value that presented a brain with gyri
 ### lower values result in a smoother surface
  brain <- contour3d(template, x=1:dtemp[1], y=1:dtemp[2], 
    z=1:dtemp[3], level = levels[[1]], alpha = alpha[[1]], 
                     color=color[[1]], draw = FALSE)

 ### Example data courtesy of Daniel Reich 
 ### Each visit is a binary mask of lesions in the brain
 scene <- list(brain)
 ## loop through images and thresh
 nimgs <- length(files)
 cols <- rainbow(nimgs)
  ## start by 2 because brain always first file
 for (iimg in 2:nimgs) {
   mask <- readNIfTI(files[iimg], reorient=FALSE)
   if (length(dim(mask)) > 3) mask <- mask[,,,1] 
    ### use 0.99 for level of mask - binary
      activation <- makeScene(mask, cutoffs = levels[[iimg]], alpha = alpha[[iimg]], cols=color[[iimg]])  
    ## add these triangles to the list
    scene <- c(scene, list(activation))
 }

  ## adding stl/obj  
  objtype <- paste0(".", objtype)
  
  fnames <- paste0(fnames, objtype)
  
  write4D(scene=scene, fnames=fnames, outfile=outfile, ...)
}