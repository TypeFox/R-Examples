"findB.design" <-
function(m,alpha,e0,e1){
 r1<-floor( alpha*(m+1) )
 r0<-ceiling( (1-alpha)*(m+1) )
 
 Bvalue<-function(M,i,w,P=alpha){
	sqrt(i/M) * (w/i-P)/sqrt(P*(1-P)/i)  }	

theta<-matrix(NA,r0+r1,2)

## first do top of design
row.num<-0	
bmin<-1
for (a in 0:(r0-1)){
	row.num <-row.num+1
	keep.doing.b.loop<-TRUE
      # Note: there is much room for efficiency gains here: take out b loop, use vector for b, perhaps use break(), 
	for (b in bmin:r1){
		if (keep.doing.b.loop==TRUE && Bvalue(M=m,i=a+b,w=b,P=alpha) >= qnorm(1-e1) ){
			theta[row.num,1]<- b
			theta[row.num,2]<- a+b
			bmin<-b
			keep.doing.b.loop<-FALSE
		}
		if (keep.doing.b.loop==TRUE && b==r1){
			theta[row.num,1]<- b
			theta[row.num,2]<- a+b
		}
	}
}	

## now do right hand side of design
amin<-1
for (b in 0:(r1-1)){
	row.num <-row.num+1
	keep.doing.a.loop<-TRUE
	for (a in amin:r0){
		if (keep.doing.a.loop==TRUE && Bvalue(M=m,i=a+b,w=b,P=alpha) <= qnorm(e0) ){
			theta[row.num,1]<- b
			theta[row.num,2]<- a+b
			amin<-a
			keep.doing.a.loop<-FALSE
		}
		if (keep.doing.a.loop==TRUE && a==r0){
			theta[row.num,1]<- b
			theta[row.num,2]<- a+b
		}
	}
}	


	
out<-list(theta=theta)
out	
}

