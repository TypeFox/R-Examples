"pCBB4" <-
function(theta,delta,s,t)
{n<-min(length(s),length(t));
	pc<-pcopula2(theta,delta,psiKS,1,ivpsiKS,1,s,t);
	for(i in 1:n){if(pc[i]=="NA"|pc[i]=="NaN"){pc[i]<-0}};
	resu<-pc
}

