# demonstrate usage of the complex image class
# 
# Author: Allison Irvine
###############################################################################
#load image
data(earth)
#initialize object
obj<- new("CmplxIm", img=img)
#set the moment type to generalized Pseudo-Zernike
momentType(obj)<- "gpzm"
#set the order
setOrder(obj)<- 100
#set the parameter  
setParams(obj)<- 1
#calculate moments of the image
Moments(obj)<- NULL
#reconstruct the image from moments
Reconstruct(obj)<- c(80,80)

#display the original image, moments, and reconstructed image
displayImg(list(obj@I,abs(obj@moments),obj@reconstruction))

