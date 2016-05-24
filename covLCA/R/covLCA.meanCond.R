covLCA.meanCond <-
function(a,g,z,J,K.j,R,S2,N)
{
	gama <- array(g,dim=c(J,K.j[1]-1,R)) #A: gamma_mkj
	
	alpha <- array(a,dim=c(J,K.j[1]-1,S2)) #A: alpha_mkp
	
	p=array(dim=c(J,K.j[1],R)) #A: p_mkj
	
	z.mean=apply(z,2,mean)
	for (m in 1:J) #A: for each manifest variable
	{
		#p[m,,]=exp(rbind(gama[m,,],0)+matrix(rep(z.mean%*%t(rbind(alpha[m,,],0)),R),nrow=K.j[1],ncol=R)) / matrix(rep(apply(exp(rbind(gama[m,,],0)+matrix(rep(z.mean%*%t(rbind(alpha[m,,],0)),R),nrow=K.j[1],ncol=R)),2,sum),rep(K.j[1],R)),nrow=K.j[1],ncol=R)
		if (dim(alpha)[3]!=1) {p[m,,]=exp(rbind(gama[m,,],0)+matrix(rep(z.mean%*%t(rbind(alpha[m,,],0)),R),nrow=K.j[1],ncol=R)) / matrix(rep(apply(exp(rbind(gama[m,,],0)+matrix(rep(z.mean%*%t(rbind(alpha[m,,],0)),R),nrow=K.j[1],ncol=R)),2,sum),rep(K.j[1],R)),nrow=K.j[1],ncol=R)}
		if (dim(alpha)[3]==1) {p[m,,]=exp(rbind(gama[m,,],0)+matrix(rep(z.mean*c(alpha[m,,],0),R),nrow=K.j[1],ncol=R)) / matrix(rep(apply(exp(rbind(gama[m,,],0)+matrix(rep(z.mean*c(alpha[m,,],0),R),nrow=K.j[1],ncol=R)),2,sum),rep(K.j[1],R)),nrow=K.j[1],ncol=R)}

	}
	return(p)
}
