tanova<-function(data, f1, f2, tp, B=100, FDR=0.05,robustify=FALSE, equal.size=FALSE, qt=-1, longitudinal=TRUE, test.type=0,eb=FALSE,df=0){
	type<-test.type
	if (length(tp)==1){
		if (type==0){
			return(gene.classifier1(data,f1,f2,B=B,FDR=FDR,robustify=robustify, equal.size=equal.size,qt=qt,eb=eb))
		}
		if (type!=0){
			o<-NANOVA.test(data=data,f1=f1,f2=f2,robustify=robustify, equal.size=equal.size, type=type,B=B,eb=eb)
			fdr<-fdr.table(o)
			c<-sig.number(fdr, FDR=FDR, qt=qt)
			s<-o$gene.order[1:c]
			p<-o$pvalue
			dt<-o$delta
			return(list(genes=s, pvalue=p[s],delta=dt[s],obj=o))
		}
	}
	if (longitudinal==TRUE & length(tp)>1){
		temp<-data.form(data=data,f1=f1,f2=f2,tp=tp)
		d<-temp$d
		fc1<-temp$fc1
		fc2<-temp$fc2
		if (type==0){
			return(gene.classifier2(data=d,f1=fc1,f2=fc2,B=B,FDR=FDR,robustify=robustify, equal.size=equal.size,time.course=length(unique(tp)),qt=qt,eb=eb,df=df))
		}
		if (type!=0){
			o<-NANOVA.test2(data=d,f1=fc1,f2=fc2,time.course=length(unique(tp)),robustify=robustify, equal.size=equal.size, type=type,B=B,eb=eb,df=df)
			fdr<-fdr.table(o)
			c<-sig.number(fdr, FDR=FDR, qt=qt)
			s<-o$gene.order[1:c]
			p<-o$pvalue
			dt<-o$delta
			a<-proj.dir(data=d,f1=fc1,f2=fc2,time.course=length(unique(tp)),type=type,df=df)
			return(list(genes=s, pvalue=p[s],delta=dt[s],a=a[s,],obj=o,dir=a))
		}
	}
	if (longitudinal==FALSE & length(tp)>1){
		if (type==0){
			return (gene.classifier3(data=data,f1=f1,f2=f2,tp=tp,B=B,FDR=FDR,qt=qt,robustify=robustify,eb=eb))
		}
		if (type!=0){
			o<-NANOVA.test3(data=data,f1=f1,f2=f2,tp=tp,type=type,B=B,robustify=robustify,eb=eb)
			fdr<-fdr.table(o)
			c<-sig.number(fdr, FDR=FDR, qt=qt)
			s<-o$gene.order[1:c]
			p<-o$pvalue
			dt<-o$delta
			tm<-0
			if (robustify==TRUE){
				tm<-0.2
			}
			a<-proj.dir2(data=data,f1=f1,f2=f2,tp=tp,type=type,trim=tm)
			return(list(genes=s, pvalue=p[s],delta=dt[s],a=a[s,],obj=o,dir=a))
}
}
}
