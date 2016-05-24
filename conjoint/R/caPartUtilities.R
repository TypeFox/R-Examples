caPartUtilities<-function(y,x,z)
{
 	options(contrasts=c("contr.sum","contr.poly"))
	outdec<-options(OutDec="."); on.exit(options(outdec))
	options(OutDec=",")
	y<-m2v(y)
	m<-length(x)
	n<-nrow(x)
	S<-nrow(y)/n
	xnms<-names(x)
	Lj<-vector("numeric", m)
	for (j in 1:m) {Lj[j]<-nlevels(factor(x[[xnms[j]]]))}
	p<-sum(Lj)-m+1
	xtmp<-paste("factor(x$",xnms,sep="",paste(")"))
	xfrm<-paste(xtmp,collapse="+")
	uslall<-round(partutilities(xfrm,y,x,n,p,S,z),3)
	return(uslall)
}