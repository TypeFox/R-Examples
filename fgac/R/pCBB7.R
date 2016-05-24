"pCBB7" <-
function(theta,delta,s,t)
{n<-min(length(s),length(t));
	pc<-pcopula1(theta,delta,psiKS,phiBB7,ivpsiKS,ivphiBB7,s,t);
	for(i in 1:n){if(pc[i]=="NA"|pc[i]=="NaN"){pc[i]<-0}};
	resu<-pc
}

