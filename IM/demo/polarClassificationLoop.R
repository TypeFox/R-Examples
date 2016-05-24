# classify bacteria light scatter images based on discrete chebyshev moment invariants of the polar unwrapped images 
# Not using MultiIm class, instead using a loop and individual OrthIm objects for each image
# Author: Allison Irvine
###############################################################################


#load images
data(bacteria);

#### feature extraction 
#analyze each image as a separate OrthIm object
invariants = list();
for (i in 1:length(img)) {
	#create a OrthIm object for each image
	obj=new("OrthIm",img=img[[i]]);
	#perform polar unwrapping on images
	transform(obj) <- 8;
	#set moment type
	momentType(obj) <- "cheby";
	#set order
	setOrder(obj) <- c(10,10)
	#setOrder(obj) <- round(dim(obj@I)/10);
	#calculate moment invariants
	Moments(obj) <-NULL
	#save invariants (moments) as vector
	invariants[[i]] = c(obj@moments);
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



