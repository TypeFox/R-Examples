#########################################################################
#
# Copyright Weierstrass Institute for Applied Analysis and 
#           Stochastics (WIAS) & Humboldt Universitaet zu Berlin, 
#           Institut fuer Mathematik, Germany 2006
# *********************************************************
#
# Name:          eval.r
#                ---------------
# Author:        Joern Schulz
# Stand:         25.08.2006
#
#########################################################################


norm <- function( orgImage, testImage, mode="L2")
{
# ========================================================================
# Calculates the $L_1$--norm or $L_2$--norm between the Original image
# 'orgImage' and the test image. 
# Note, it is assumed that the data have the same size. 'orgImage' and 
# 'testImage' have to be of type 'vector' or 'matrix'.
# Possibililties for 'mode': 'L1', 'L2', 'PT-L1', 'PT-L2'
# 

  if (is.vector(orgImage)){
     if (!(is.vector(testImage)))
        stop("The data have to be of the same type, but 'orgImage' is of type 'vector' and 'testimage' not.") 
     n <- length(orgImage)
     if(n != length(testImage)){
        stop("The vectors should have the same size.")
     }
  } else if(is.matrix(orgImage)){
    if (!(is.matrix(testImage)))
        stop("The data have to be of the same type, but 'orgImage' is of type 'matrixr' and 'testimage' not.") 
    n1 <- nrow(orgImage)
    n2 <- ncol(orgImage)
    n <- n1 * n2
    if ( n1 != nrow(testImage) || n2 != ncol(testImage))
       stop("The matrices should have the same size.")
  } else
    stop("Invalid type of data.")

  if ( mode == "L1" ){
     LP <- (sum(abs(testImage-orgImage)))/(n)
  } else if ( mode == "L2" ){
     LP <- sqrt(sum((testImage-orgImage)^2)/n)
  } else if ( mode == "PT-L1" ){
     mdata <- mean(orgImage)
     DC <- mean(testImage)- mdata
     LP <- (sum(abs(testImage-orgImage+DC)))/(n*mdata)
  } else if ( mode == "PT-L2" ){
     mdata <- mean(orgImage)
     DC <- mean(testImage)- mdata
     LP <- sqrt((sum((testImage-orgImage-DC)^2))/n)
     DV <- sqrt((sum((orgImage-mdata)^2))/n)
     LP <- LP/DV
  } else
    stop("The 'mode' = ", mode, " is not supported.")

  return(LP)
}


scaleImage <- function( image, mode="mean" )
{
# ========================================================================
# Scales the minimum value of an image to 0 and the mean to 1.
# 

  imMin <- min(image)
  image <- image - imMin
  
  if (mode=="mean"){
    if (imMin < 0){
        #image <- image + abs(imMin)
        image <- image/mean(image)
    } else if (imMin > 0){ 
        #image <- image - imMin
        image <- image/mean(image)
    } else
        image <- image/mean(image)
  } else if (mode=="max"){
    if (imMin < 0){
        #image <- image + abs(imMin)
        image <- image/max(image)
    } else if (imMin > 0){ 
        #image <- image - imMin
        image <- image/max(image)
    } else
        image <- image/max(image)
  } else {
      cat("'mode=",mode,"' is not supported. Choose mode='mean' or 'max'. \n", sep="")
      image <- image + imMin
  }

  return(image)
}


infoImage <- function(image)
{
# ========================================================================
# Print different info statments about the image.
#
cat("Number of pixels of image:", prod(dim(image)), "\n")
cat("Min:", min(image), "Max:", max(image), "Mean:", mean(image), "\n")
cat("Quantile: \n")
print( quantile(image)) #"0.75-Quantile:", quantile(image,3/4), 
#cat("IQR:", IQR(image), "\n")
# sd "Standard Abweichung, var, "variance""

}