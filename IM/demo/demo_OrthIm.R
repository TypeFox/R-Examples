# demonstrate usage of the orthogonal image class
# 
# Author: Allison Irvine
###############################################################################
#load image
data(mandril)
#initialize object
obj <- new("OrthIm",img = img,filename = "");
#set the moment type to bivariate - continuous Chebyshev and Legendre
momentType(obj) = c("chebycont","legend");
#set the order
setOrder(obj) = c(150,150);
# no parameters are required for continuous Chebyshev of Legendre polynomials.
#calculate moments of the image
Moments(obj) = NULL;
#reconstruct the image from moments up to order 125 in both the x and y directions
Reconstruct(obj) = c(125,125);
#display the original image, moments, and reconstructed image
displayImg(list(obj@I,obj@moments,obj@reconstruction))

#display the polynomials of order 100 and 150 used to calculate the moments
plotPoly(obj,order=c(100,150));

