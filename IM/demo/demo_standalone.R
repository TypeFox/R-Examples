# demonstrate moment calculation using the standalone function, which returns an object
# 
# Author: Allison Irvine
###############################################################################
#load image
data(pirate)
#calculate discrete Chebyshev moments up to order 100 in the x and y directions
obj = momentObj(I=img,type="cheby",order=c(100,100));
#display moments
plotMoment(obj)
#display the polynomials used to calculate the moments, of orders 11 through 15
plotPoly(obj,order=11:15)
#reconstruct the image using moments up to order 75 in the x and y directions
Reconstruct(obj) = c(75,75);
#display reconstructed image
displayImg(obj@reconstruction)
