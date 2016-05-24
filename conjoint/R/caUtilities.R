caUtilities<-function(y,x,z)
{
 	options(contrasts=c("contr.sum","contr.poly"))
	outdec<-options(OutDec="."); on.exit(options(outdec))
	options(OutDec=",")
	y<-m2v(y)
	m<-length(x)
	n<-nrow(x)
	S<-nrow(y)/n
	xnms<-names(x)
	ynms<-names(y)
	xtmp<-paste("factor(x$",xnms,sep="",paste(")"))
	xfrm<-paste(xtmp,collapse="+")
	yfrm<-paste("y$",ynms,sep="","~")
	frml<-as.formula(paste(yfrm,xfrm))
	Lj<-vector("numeric",m)
	for(j in 1:m) {Lj[j]<-nlevels(factor(x[[xnms[j]]]))}
	x<-as.data.frame(matexpand(m,n,S,x))
	camodel<-lm(frml)
	print(summary.lm(camodel))
	u<-as.matrix(camodel$coeff)  
	intercept<-u[1]
	ul<-utilities(u,Lj)
	utlsplot(ul,Lj,z,m,xnms)
	uli<-c(intercept,ul)
	return(uli)
}
