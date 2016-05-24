# classify Tamil characters using discrete Chebyshev moments
# 
# Author: Allison Irvine, Tan Dang
###############################################################################

#load images
data(characters)

#create MultiIm object to store all images
obj=new("MultiIm",images=img)

#set type and order
momentType(obj) <- "cheby"
setOrder(obj) <- c(5,5);

#calculate moments
Moments(obj) <- NULL

#vectorize moments
moments=list();
for (i in 1:length(labels)) {
	moments[[i]] = c(obj@moments[[i]]);
}
moments=t(do.call(cbind,moments));

#separate training and testing data (stratified random sampling) About 2/3 from every category will be in the training set.
trainIndex = logical(length(labels));
uniqueY = unique(labels);
for (i in 1:length(uniqueY)) {
	trainIndex[sample(which(labels==uniqueY[i]),26)]=TRUE;
}
trX = moments[trainIndex,];
trY = labels[trainIndex];
tstX = moments[!trainIndex,];
tstY = labels[!trainIndex];

#apply svm classifier
library(e1071)
model = svm(x=trX,y=as.factor(trY))
pred = predict(model,tstX)

#calculate error
err = sum(tstY!=pred)/length(pred)




