"AdaptPredmp" <-
function(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g){

# does local adaptive prediction for the point remove based on N 
# points (chooses method of prediction and intercept)


details<-NULL
results<-list()
d<-matrix(0,1,length(coefflist[[remove]]))

p<-list()
w<-list()

intercept<-FALSE  #does prediction schemes with no intercept

out<-LinearPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)
p[[1]]<-out$pred
w[[1]]<-out$weights

out<-QuadPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)
p[[2]]<-out$pred
w[[2]]<-out$weights

out<-CubicPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)
p[[3]]<-out$pred
w[[3]]<-out$weights

intercept<-TRUE

out<-LinearPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)
p[[4]]<-out$pred
w[[4]]<-out$weights

out<-QuadPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)
p[[5]]<-out$pred
w[[5]]<-out$weights

out<-CubicPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)
p[[6]]<-out$pred
w[[6]]<-out$weights

pre<-matrix(0,1,6)
for (k in 1:6){
	pred<-p[[k]]
	if (length(pred)>1){
		md<-matrix(0,1,length(coefflist[[remove]]))
		pr<-matrix(0,1,length(coefflist[[remove]]))
		if (mpdet=="min"){
			for (i in 1:length(coefflist[[remove]])){
				pr[i]<-order(abs(coefflist[[remove]][i]-pred))[1]
				md[i]<-(coefflist[[remove]][i]-pred)[pr[i]]
			} 
		}
		else{
			for (i in 1:length(coefflist[[remove]])){
			md[i]<-mean(coefflist[[remove]][i]-pred)
			}
		}
	
		if (mpdet=="min"){
			sel<-order(abs(md))[1]
			details[k]<-md[sel]
			pre[k]<-pred[pr[sel]]
		}
		else{
			details[k]<-mean(md)
			pre[k]<-mean(pred)
		}


	}
	else{
		for (i in 1:length(coefflist[[remove]])){
			d[i]<-coefflist[[remove]][i]-pred
		}
		aved<-mean(d)
		mind<-min(d)
		pre[k]<-pred	
		if (mpdet=="min"){
			details[k]<-mind
		}
		else{
			details[k]<-aved
		}

	} 

} 

minindex<-order(abs(details))[1]
pred<-pre[minindex]
coefflist[[remove]]<-details[minindex]
coeff[remove]<-coefflist[[remove]]

int<-NULL
scheme<-NULL

if(minindex<=3){
	int<-FALSE
}
else{
	int<-TRUE
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

results[[1]]<-w[[minindex]]
results[[2]]<-pred
results[[3]]<-coeff
results[[4]]<-int
results[[5]]<-scheme
results[[6]]<-details
results[[7]]<-minindex

results
}
