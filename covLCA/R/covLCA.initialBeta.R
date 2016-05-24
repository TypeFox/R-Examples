covLCA.initialBeta <-
function(yy,RR=2,xx,maxit=50000,nrp=20) 
{

	stopifnot(require(mlogit))
	stopifnot(require(poLCA))
	namy=colnames(yy)
	colnames(xx)[2:(dim(xx)[2])]=namx=paste("x",1:(dim(xx)[2]-1),sep="")
	
	fm=as.formula(paste("cbind(",paste(namy,collapse=","),")~1"))

	dataRed=data.frame(yy,xx)
	res=poLCA(fm,dataRed,nclass=RR,nrep=nrp,calc.se=FALSE,maxiter=maxit,verbose=FALSE) 
	
	attr=res$predclass
	
	attr=as.factor(attr)
	xwinter=xx[,2:(dim(xx)[2])]
	data1=data.frame(attr,xwinter)
	colnames(data1)[2:dim(data1)[2]]=namx 
	
	fm2=as.formula(paste("attr~1|",paste(namx,collapse="+")))
	data2<-mlogit.data(data1,varying=NULL,choice="attr",shape="wide")
	mlogit.model<-mlogit(fm2,data=data2,reflevel=max(levels(attr))) # reference latent class : the last one
	retBeta=matrix(mlogit.model$coef,nrow=(RR-1)) #A: matrix where rows=LC, columns=covariates
	Beta=c(t(retBeta)) #A: vector with beta_jp, jp=11,...,1P, 21,...,2P, (J-1)1,...,(J-1)P
	
	return(Beta)
}
