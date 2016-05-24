gene.classifier3<-function(data, f1, f2, tp, B=100,FDR=0.05, qt=-1,robustify=FALSE,eb=FALSE){
	p<-length(f1)
	n<-dim(data)[1]
	
	s<-c(1:n)
	o2<-NANOVA.test3(data,f1=f1,f2=f2,tp=tp,type=2, B=B,robustify=robustify,eb=eb)
	fdr<-fdr.table(o2)
	c<-sig.number(fdr, FDR=FDR, qt=qt)
	ix<-o2$gene.order[1:c]
	s2<-ix[o2$pvalue[ix]<=0.05]
	p2<-o2$pvalue
	dt2<-o2$delta
	
	index<-s2
	
	o1<-NANOVA.test3(data[index,],f1=f1,f2=f2, tp=tp,type=1,B=B,robustify=robustify,eb=eb)
	fdr<-fdr.table(o1)
	c<-sig.number(fdr, FDR=FDR, qt=qt)
	ix<-o1$gene.order[1:c]
	ix<-ix[o1$pvalue[ix]<=0.05]
	s1<-s2[ix]
	p1<-o1$pvalue[ix]
	dt1<-o1$delta[ix]
	m1<-cbind(s1,p1,dt1)
	
	ix<-setdiff(s2,s1)
	
	index<-ix
	
	o3<-NANOVA.test3(data[index,],f1=f1,f2=f2,tp=tp,type=3,B=B,robustify=robustify,eb=eb)
	fdr<-fdr.table(o3)
	c<-sig.number(fdr, FDR=FDR, qt=qt)
	ixx<-o3$gene.order[1:c]
	ixx<-ixx[o3$pvalue[ixx]<=0.05]
	s3<-ix[ixx]
	p3<-o3$pvalue[ixx]
	dt3<-o3$delta[ixx]
	
	o4<-NANOVA.test3(data[index,],f1=f1,f2=f2, tp=tp,type=4,B=B,robustify=robustify,eb=eb)
	fdr<-fdr.table(o4)
	c<-sig.number(fdr, FDR=FDR, qt=qt)
	ixx<-o4$gene.order[1:c]
	ixx<-ixx[o4$pvalue[ixx]<=0.05]
	s4<-ix[ixx]
	p4<-o4$pvalue[ixx]
	dt4<-o4$delta[ixx]
	
	t<-setdiff(s3,intersect(s3,s4))
	if (length(t)>0){
		k<-s3%in%t
		m3<-cbind(s3,p3,dt3)[k,]
		if (is.matrix(m3)==FALSE){
			m3<-matrix(m3,nrow=1)
		}
	}
	if (length(t)==0){
		message('C3=0')
		m3<-matrix(c(0,0,0),nrow=1)
	}
	
	t<-setdiff(s4,intersect(s3,s4))
	if (length(t)>0){
		k<-s4%in%t
		m4<-cbind(s4,p4,dt4)[k,]
		if (is.matrix(m4)==FALSE){
			m4<-matrix(m4,nrow=1)
		}
	}
	if (length(t)==0){
		message('C4=0')
		m4<-matrix(c(0,0,0),nrow=1)
	}
	
	k<-setdiff(s2,union(union(m3[,1],m4[,1]),s1))
	if (length(k)>0){
		m2<-cbind(s,p2,dt2)[k,]
		if (is.matrix(m2)==FALSE){
			m2<-matrix(m2,nrow=1)
		}
	}
	if (length(k)==0){
		message('C2=0')
		m2<-matrix(c(0,0,0),nrow=1)
	}
	
	S5<-setdiff(s,union(union(union(m1[,1],m2[,1]),m3[,1]),m4[,1]))
	
	C1<-m1[,1]
	C2<-m2[,1]
	C3<-m3[,1]
	C4<-m4[,1]
	
	a1<-proj.dir2(data[C1,],f1=f1,f2=f2,tp=tp,type=1,robustify=robustify)
	
	a2<-proj.dir2(data[C2,],f1=f1,f2=f2,tp=tp,type=2,robustify=robustify)
	
	a3<-proj.dir2(data[C3,],f1=f1,f2=f2,tp=tp,type=3,robustify=robustify)
	
	a4<-proj.dir2(data[ix,],f1=f1,f2=f2,tp=tp,type=4,robustify=robustify)
	
	return(list(C1=m1[,1],C2=m2[,1], C3=m3[,1],C4=m4[,1], C1.pvalue=m1[,2],C2.pvalue=m2[,2], C3.pvalue=m3[,2],C4.pvalue=m4[,2],C1.delta=m1[,3],C2.delta=m2[,3],C3.delta=m3[,3],C4.delta=m4[,3],a1=a1,a2=a2,a3=a3,a4=a4))
}
