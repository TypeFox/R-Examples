# compare reconstruction error of different real orthogonal moment types
# 
# Author: Allison Irvine
###############################################################################

data(earth_200)
orders=c(25,50,100)

#####orthogonal moments#####
#discrete
#"cheby"		Chebyshev
#"krawt"		Krawtchouk					param = 0.5
#"hahn"			Hahn						param = 0,0.5
#CONTINUOUS:
#"gegen"		Gegenbauer					param = 1  (not 0)
#"chebycont"	Continuous Chebyshev	
#"legend"		Legendre

types = c("cheby","krawt","hahn","gegen","chebycont","legend")
params = list(0,0.5,c(dim(img)[1]/2,dim(img)[1]/5),2,0,0)
obj = new("OrthIm",img=img)
displayImg(obj@I)

#results will be stored in a matrix, rows are moment type, cols are order
#sum of absolute error
error = mat.or.vec(length(types),length(orders))
#rms error
error1 = mat.or.vec(length(types),length(orders))

#total number of pixels
nPix = dim(obj@I)[1]*dim(obj@I)[2]

#for each moment type
for (i in 1:length(types)) {
	momentType(obj) = types[i];
	#for each order in the list
	for (j in 1:length(orders)) {
		setOrder(obj) = c(orders[j],orders[j])
		setParams(obj) = list(params[[i]],params[[i]]);
		Moments(obj) = 0;
		Reconstruct(obj) = c(orders[j],orders[j]);
		#display reconstructed image
		#displayImg(obj@reconstruction)
		#title(sprintf("%s order %d",types[i],orders[j]))
		
		#normalize both image pixel values between 0 and 1
		I1 = (obj@I-min(obj@I)) / (max(obj@I)-min(obj@I))
		I2 = (obj@reconstruction-min(obj@reconstruction)) / (max(obj@reconstruction)-min(obj@reconstruction))
		
		#calculate percent reconstruction error relative to original image
		error[i,j] = sum(abs(I1-I2))/sum(I1);
		#calculate rms error of reconstructed image
		error1[i,j] = sqrt(sum((I1-I2)^2)/nPix);
	}
}


