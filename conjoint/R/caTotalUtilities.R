caTotalUtilities<-function(y,x)
{
 	options(contrasts=c("contr.sum","contr.poly"))
	outdec<-options(OutDec="."); on.exit(options(outdec))
	options(OutDec=",")
	y<-m2v(y)
	m<-length(x)
	n<-nrow(x)
	S<-nrow(y)/n
	xnms<-names(x)
	xtmp<-paste("factor(x$",xnms,sep="",paste(")"))
	xfrm<-paste(xtmp,collapse="+")
	Lj<-vector("numeric",m)
	for(j in 1:m) {Lj[j]<-nlevels(factor(x[[xnms[j]]]))}
	p<-sum(Lj)-m+1
	Usi<-round(totalutilities(xfrm,y,x,n,p,S),3)
	return(Usi)
}