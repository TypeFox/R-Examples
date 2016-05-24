# classify bacteria light scatter images based on generalized pseudo-zernike moment invariants 
# same as the cmplxClassification demo
# Not using MultiIm class, instead using a loop and individual CmplxIm objects for each image
# Author: Allison Irvine
###############################################################################


#load images
data(bacteria);

#### feature extraction 
#analyze each image as a separate CmplxIm object
invariants = list();
for (i in 1:length(img)) {
	#create a CmplxIm object for each image
	obj=new("CmplxIm",img=img[[i]]);
	#set moment type
	momentType(obj) <- "gpzm";
	#set order
	setOrder(obj) <- 50;
	#set scaling parameter for gpzm
	setParams(obj) <- 2;
	#calculate moment invariants
	Moments(obj) <-NULL
	Invariant(obj) <- NULL
	#save invariants as vector
	#everything above upper triangular in matrix is zeros, only take lower triangular and diagonal 
	invariants[[i]] = c(obj@invariant[lower.tri(obj@invariant, diag = TRUE)]);
}
invariants=t(do.call(cbind,invariants));


#### classification/validation
#separate training and testing data (stratified random sampling) 3 out of 5 from every category
trainIndex = logical(length(labels));
uniqueY = unique(labels);
for (i in 1:length(uniqueY)) {
	trainIndex[sample(which(labels==uniqueY[i]),3)]=TRUE;
}
trX = invariants[trainIndex,];
trY = labels[trainIndex];
tstX = invariants[!trainIndex,];
tstY = labels[!trainIndex];

#apply svm classifier
library(e1071)
model = svm(x=trX,y=as.factor(trY))
pred = predict(model,tstX)

#calculate error
err = sum(tstY!=pred)/length(pred)



