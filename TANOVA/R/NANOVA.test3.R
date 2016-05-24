NANOVA.test3<-function(data,f1,f2,tp,type=2,B=100,robustify=FALSE,eb=FALSE){
	if (robustify==FALSE){
		tm<-0
	}
	if (robustify==TRUE){
		tm<-0.2
	}
	F<-vector()
	F.null<-matrix()
	F<-F.stat2(data,f1,f2,tp=tp,type=type,trim=tm,eb=eb)
	F.null<-F.stat.null2(data, f1, f2, tp=tp,type=type,B=B,trim=tm,eb=eb)
	delta<-z.score(F,F.null)$z
	gene<-sort(delta,decreasing=TRUE,index=TRUE)$ix
	p<-vector()
	n<-length(F)
	for (i in 1:n){
		p[i]<-1-ecdf(F.null[i,])(F[i])
	}
	list(gene.order=gene, F=F, F.null=F.null, pvalue=p, delta=delta)
}