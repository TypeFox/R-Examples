NANOVA.test<-function(data,f1,f2,type=2,B=100, robustify=FALSE,equal.size=FALSE,eb=FALSE){
	if (robustify==FALSE){
		tm<-0
	}
	if (robustify==TRUE){
		tm<-0.2
	}
	F<-vector()
	F.null<-matrix()
	F<-F.stat(data,f1,f2,type=type,trim=tm,equal.size=equal.size,eb=eb)
	F.null<-F.stat.null(data, f1, f2, type=type,B=B, trim=tm,equal.size=equal.size,eb=eb)
	temp<-z.score(F,F.null)
	delta<-temp$z
	gene<-sort(delta,decreasing=TRUE,index=TRUE)$ix
	p<-vector()
	n<-length(F)
	for (i in 1:n){
		p[i]<-1-ecdf(F.null[i,])(F[i])
	}
	list(gene.order=gene, F=F, F.null=F.null, pvalue=p, delta=delta)
} 