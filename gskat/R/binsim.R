#Author: Xuefeng Wang
#Based on PARK et al. (1996), only work for positive R 


mybin<- function (p,R) {
	q=1-p
	k=nrow(R)
	#step 0
	a=matrix(NA,k,k)
	for (i in 1:k) {
		for (j in 1:k) {
			a[i,j]=log(1+R[i,j]*sqrt(q[i]*(p[i])^(-1)*q[j]*(p[j])^(-1)))
		}
	}
	l=0
	beta=0
	myS=list()

	while (sum(a)!=0 ){	

		#step 1
		l=l+1
		beta[l]=min(a[a>0])
		r=max((which(a==beta[l]))) %% k  ;  if (r==0){r=k}
		s=ceiling((max((which(a==beta[l]))) / k))

		if (a[r,r]==0 | a[s,s]==0 ){stop("feasibility check stop")}

		S=unique(c(r,s))
		for (i in 1:k) {
			myprod=1
			for (j in 1:length(S)) {
				myprod<- myprod* a[i,S[j]]
			}
			if ( myprod >0 ) {S=unique(c(S,i))}	
		}
		S=sort(S)

		a

		#step 2
		for (i in S) {
			for (j in S) {
				a[i,j]=a[i,j]-beta[l]
			}
		}
		a=round(a,10)
		myS[[l]]=S
		r
		s
		S
		beta[l]
	}

	#step 3
	X=rpois(l,lambda=beta) 
	I=matrix(0,ncol=l,nrow=k)    #index matrix
	for (i in 1:l) {I[myS[[i]],i]=1}
	Y=I%*%X
	Z<- as.numeric(Y==0)

  return(Z)
}


#Frechet bounds check
frechet<-function (p,r) {
  q=1-p
  n=nrow(r)
  lower<-matrix(NA,n,n); upper<-matrix(NA,n,n)
  check<-matrix(NA,n,n);
  for (j in 1:n) {
	for (k in 1:n) {
		lower[j,k]<-max(-sqrt((p[j]*p[k])/(q[j]*q[k])), -sqrt((q[j]*q[k])/(p[j]*p[k])))
		upper[j,k]<-min(sqrt((p[j]*q[k])/(p[k]*q[j])), sqrt((p[k]*q[j])/(p[j]*q[k])))
		check[j,k]<-(r[j,k]>=lower[j,k])& (r[j,k]<=upper[j,k])
		#if (check[j,k]==FALSE ){stop("Frechet bounds not meet")}
	}  
  }
  if (all(check)){return("OK")} else {return("NOK")}
}
