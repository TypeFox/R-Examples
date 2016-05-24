fdr.table<-function(obj){
	gene.order<-obj$gene.order
	delta<-sort(obj$delta,decreasing=TRUE)
	F<-obj$F
	F.null<-obj$F.null
	z.null<-z.score(F,F.null)$z.null
	n<-dim(F.null)[1]
	B<-dim(F.null)[2]
	s<-seq(from=1,to=n-1,by=5)
	r<-length(s)
	temp<-vector()
	for (i in 1:B){
		temp[i]<-sum(delta>max(z.null[,i]))/n
	}
	pi0<-1-median(temp)
	p<-c(0.25,0.5,0.75,0.9)
	table<-matrix(nrow=r,ncol=length(p)+2)
	for (i in 1:r){
		c<-delta[s[i]+1]
		rej<-vector()
		for (j in 1:B){
			rej[j]<-sum(z.null[,j]>c)
		}
		table[i,1]<-s[i]
		table[i,2:5]<-quantile(pi0*rej/s[i],probs=p,na.rm=TRUE)
		table[i,6]<-mean(pi0*rej/s[i])
}
return (table)
}