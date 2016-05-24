MergeBreakPointsByChr <- function(chr,id,seg.B,seg.L,map.chr,thr=5){
	m<-chr
	ss1.chr<-seg.B
	ss2.chr<-seg.L
	ss1.chr<-ss1.chr[order(as.numeric(ss1.chr[,3])),]
	ss2.chr<-ss2.chr[order(as.numeric(ss2.chr[,3])),]
	tmp<-as.numeric(unlist(t(ss2.chr[,3:4])))
	cmp<-as.numeric(unlist(t(ss1.chr[,3:4])))
	cmp<-getUnifiedMap(tmp,cmp,map.chr,thr=thr)

	#Use logR segment as template and BAF as compare
	#Inter-logR segment regions are truly empty, while BAF not.

	flist<-c()
	j<-1
	Pc<-cmp[1]
	for(n in seq(1,length(tmp),2)){
		Pt1<-tmp[n]
		Pt2<-tmp[n+1]
		cur.i<-Pt1
		while(cur.i<Pt2){
			while(Pc<=cur.i){
				#flist=c(flist,cmp[j],cmp[j+1])						
				j<-j+1
				if(j>length(cmp)){break}
				Pc<-cmp[j]
			}
			if(j>length(cmp)){break}
			if(Pc>Pt2){
				flist<-c(flist,cur.i,Pt2)
				cur.i<-Pt2
				next
			}

			#start of a segment

			flist<-c(flist,cur.i,Pc)
			cur.i<-Pc
		}
		if(j>length(cmp)){break}
	}		
	Ns<-length(flist)/2
	new.ss<-t(matrix(flist,ncol=Ns))
	new.ss<-cbind(rep(m,Ns),new.ss)
	new.ss<-cbind(rep(id,Ns),new.ss)
	return(new.ss)
}