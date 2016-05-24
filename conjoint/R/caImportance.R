caImportance<-function(y,x)
{
	options(contrasts=c("contr.sum", "contr.poly"))
	outdec<-options(OutDec="."); on.exit(options(outdec))
	options(OutDec=",")
	y<-m2v(y)
	m<-length(x)
	n<-nrow(x)
	S<-nrow(y)/n
	xnms<-names(x)
	Lj<-vector("numeric",m)
	for (j in 1:m) {Lj[j]<-nlevels(factor(x[[xnms[j]]]))}
	p<-sum(Lj)-m+1
	xtmp<-paste("factor(x$",xnms,sep="",paste(")"))
	xfrm<-paste(xtmp,collapse="+")
	usl<-partutils(xfrm,y,x,n,p,S)
	imps<-matrix(0,S,m)
   	for(s in 1:S)
	{
	      u<-usl[s,]
	      ul<-utilities(u,Lj)
	      imp<-importance(ul,Lj)*100
	      imps[s,]<-imp
	}
	impS<-round(apply(imps,2,"mean"),2)
	return(impS)
}