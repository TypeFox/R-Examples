gehan.test<-function(time,event,trt) {


#  number of times x clearly beats y - number of times y clearly beats x

	if( length(unique(trt)) != 2 ) stop("trt must have two levels")

	nvec<-table(trt)

	n<-sum(nvec)
	U<-matrix(0,nrow=n,ncol=n)
	D<-sign(outer(time,time,'-'))		

	ind11<-outer(event,event,'&')
	ind01<-outer(!event,event,'&')
	ind10<-outer(event,!event,'&')

	U[ind11]<-D[ind11]
	U[ind01&(D>=0)]<- 1
	U[ind10&(D<=0)]<- -1

	W<-apply(U,1,sum)

	statistic<-sum(W[trt==trt[1]])/sqrt(prod(nvec)/sum(nvec)*var(W))

	res<-list(statistic=statistic,p.value=2*pnorm(abs(statistic),lower.tail=FALSE))

	class(res)<-'rank.test'

	res

}
