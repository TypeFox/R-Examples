"AdaptPred" <-
function(pointsin,X,coeff,nbrs,remove,intercept,neighbours){

#does local adaptive prediction for the point remove based on N 
#points (chooses method of prediction and intercept);

details<-NULL
results<-list()
w<-list()

intercept<-FALSE  #does prediction schemes with no intercept

out1<-LinearPred(pointsin,X,coeff,nbrs,remove,intercept)
pred1<-out1$pred
w1<-out1$weights
w[[1]]<-w1
details[1]<-coeff[remove]-pred1

out1<-QuadPred(pointsin,X,coeff,nbrs,remove,intercept)
pred1<-out1$pred
w1<-out1$weights
w[[2]]<-w1
details[2]<-coeff[remove]-pred1

out1<-CubicPred(pointsin,X,coeff,nbrs,remove,intercept)
pred1<-out1$pred
w1<-out1$weights
w[[3]]<-w1
details[3]<-coeff[remove]-pred1

intercept<-TRUE

out1<-LinearPred(pointsin,X,coeff,nbrs,remove,intercept)
pred1<-out1$pred
w1<-out1$weights
w[[4]]<-w1
details[4]<-coeff[remove]-pred1

out1<-QuadPred(pointsin,X,coeff,nbrs,remove,intercept)
pred1<-out1$pred
w1<-out1$weights
w[[5]]<-w1
details[5]<-coeff[remove]-pred1

out1<-CubicPred(pointsin,X,coeff,nbrs,remove,intercept)
pred1<-out1$pred
w1<-out1$weights
w[[6]]<-w1
details[6]<-coeff[remove]-pred1

minindex<-order(abs(details))[1]
pred<-coeff[remove]-details[minindex]
coeff[remove]<-details[minindex]
int<-TRUE
scheme<-NULL
if(minindex<=3){
	int<-FALSE
}

if((minindex==1)|(minindex==4)){
	scheme<-"Linear"
}
if((minindex==2)|(minindex==5)){
	scheme<-"Quad"
}
if((minindex==3)|(minindex==6)){
	scheme<-"Cubic"
}

weights<-w[[minindex]]

return(list(weights=weights,pred=pred,coeff=coeff,int=int,scheme=scheme,details=details,minindex=minindex))

}
