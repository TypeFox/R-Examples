NBF <-
function(y,G,P,a,b,s2,nu){
	Z=data.matrix(G)
	P=data.matrix(P)
	X=Z%*%P
	ZZ=Z%*%t(Z)
	XX=X%*%t(X)

	BF=c()
	model0.log=logproby(y=y,XX=XX,ZZ=ZZ,aa=a,bb=b,s2=s2,nu=nu)

	i=1
	BF=foreach(i=1:ncol(X),.combine="cbind")%dopar%{
 	XX.new=XX-X[,i]%*%t(X[,i])
  	model1.log=logproby(y=y,XX=XX.new,ZZ=ZZ,aa=a,bb=b,s2=s2,nu=nu)
  	exp(model1.log-model0.log)
	}
	BF=as.numeric(BF)
	names(BF)=colnames(P)
BF
}
