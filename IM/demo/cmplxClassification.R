# classify bacteria light scatter images based on generalized pseudo-zernike moment invariants 
# 
# Author: Allison Irvine
###############################################################################


#load images
data(bacteria);

#create MultiIm object to store all images
obj=new("MultiIm",images=img)

#set parameters
momentType(obj) <- "gpzm"
setOrder(obj) <- 50;
setParams(obj) <- 2;

#calculate moment invariants
Invariant(obj) <- NULL

#vectorize moment invariants
#everything above upper triangular in matrix is zeros, only take lower triangular and diagonal 
invariants=list();
for (i in 1:length(labels)) {
	invariants[[i]] = c(obj@invariant[[i]][lower.tri(obj@invariant[[i]], diag = TRUE)]);
}
invariants=t(do.call(cbind,invariants));


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



