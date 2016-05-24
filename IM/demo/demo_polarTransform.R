# demonstrate the polar transform applied to an image
# 
# Author: Allison Irvine
###############################################################################

#load image
data(circles)
#convert to grayscale
img = rowSums(img,dim=2)
#perform polar unwrapping on image
transf = polarTransform(img, resolution=20, scale=100, center= calcCentroid(img))
#reverse transform
itransf= revPolar(dim(img), transf);
#display the original image, transformed image, and reverse transformed image
displayImg(list(img,transf[[1]],itransf))



