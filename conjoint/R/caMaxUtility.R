caMaxUtility<-function(sym,y,x)
{
	options(contrasts=c("contr.sum","contr.poly"))
	outdec<-options(OutDec="."); on.exit(options(outdec))
	options(OutDec=",")
	y<-m2v(y)
	m<-length(x)
	n<-nrow(x)
	S<-nrow(y)/n
	xnms<-names(x)
	Lj<-vector("numeric",m)
	for(j in 1:m) {Lj[j]<-nlevels(factor(x[[xnms[j]]]))}
	p<-sum(Lj)-m+1
	xtmp<-paste("factor(x$",xnms,sep="",paste(")"))
	xfrm<-paste(xtmp,collapse="+")
	usl<-partutils(xfrm,y,x,n,p,S)
	psc<-profsimcode(sym)
	Zsym<-as.matrix(psc)
	Usym<-Zsym%*%t(usl)
	r<-nrow(Zsym)
	maxutil<-maxutility(Usym,S,r)
	return(maxutil)
}