covLCA.initialAlphaGamma <-
function(yy,zz,R,K.j,S2,...)
{

	library(mlogit)
	
	retAlpha=array(dim=c(dim(yy)[2],S2*(K.j[1]-1))) #Adapted to "S2=1"
	retGamma=matrix(nrow=dim(yy)[2],ncol=R*(K.j[1]-1))

	
	for (i in 1:dim(yy)[2])
	{

		data1=data.frame(as.factor(yy[,i]),zz)
		names(data1)[1]="y"
	
		names(data1)[2:dim(data1)[2]]=nam=paste("Z",1:(dim(data1)[2]-1),sep="")
		fm=as.formula(paste("y~1| ",paste(nam,collapse="+")))
	
		data2<-mlogit.data(data1, varying=NULL, choice="y", shape="wide")
		mlogit.model<- mlogit(fm, data = data2, reflevel=K.j[1])#,...)
		retAlpha[i,]=mlogit.model$coef[(K.j[1]):(length(mlogit.model$coef))] #A: the first (K.j[1]-1) parameters are the intercepts
		#A: row m contains alpha_mpk, where pk=11,12,...,1(K-1),21,22,...,2(K-1),...,P1,P2,...P(K-1)
		#retGamma[i,]=rep(mlogit.model$coef[1:(K.j[1]-1)],rep(3,K.j[1]-1))
		retGamma[i,]=rep(mlogit.model$coeff[1:(K.j[1]-1)],R) #A: modif : OK ?
		#A: row m contains a repetition of gamma_mk: m1,m2,...,m(K-1),m1,m2,...,m(K-1),...
	}
	return(list(Alpha=retAlpha,Gamma=retGamma))
}
