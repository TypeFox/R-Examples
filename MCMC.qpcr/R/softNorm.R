softNorm=function(data,controls){
	if (is.null(controls)) { 
		stop("Error: normalize=TRUE is specified wihtout specifying controls")
	}
	ccc=c();saa=c();ccs=c();sal=c()
	for (sa in unique(data$sample)){
		sale=0
		for (cc in controls) {
			if (length(data[data$gene==cc & data$sample==sa,"count"])>sale) {
				sale=length(data[data$gene==cc & data$sample==sa,"count"])
			}
		}
		sal=append(sal,sale)
		saa=append(saa,sa)
	}
	saac=c()
	for (s in 1:length(saa)){
		sa=saa[s]
		for (cc in controls) {
			cd=data[data$gene==cc & data$sample==sa,"count"]
			for (i in 1:sal[s]) {
				if(length(cd)>=i) { 
					ccs=append(ccs,cd[i])
				} else { 
					ccs=append(ccs,NA)
				}
				ccc=append(ccc,cc)
				saac=append(saac,sa)
			}
		}
	}
	
	CTR=data.frame(cbind("count"=ccs,"gene"=ccc,"sample"=saac))
	CTR$count=as.integer(as.character(CTR$count))
	sam=unstack(CTR,sample~gene)[,1]
	CTR=data.frame(unstack(CTR,count~gene))
	CTR$sample=sam
	
	if (length(controls)>1) {
		CTR$NORM=round(apply(CTR[,controls]+1,1,prod)^(1/length(controls)),0)
	}	else { CTR$NORM=CTR[,controls[1]] }		
	
	rest=c()
	for (sa in CTR$sample){
		dss=data[data$sample==sa,]
		dss=dss[1,]
		rest=rbind(rest,dss[3:length(dss)])
	}
	rest=cbind("count"=CTR$NORM,"gene"=rep("NORM"),rest)
	data=rbind(data,rest)
	data=data[!(data$gene %in% controls),]
	data$gene=factor(data$gene,unique(data$gene))
	data$gene=relevel(data$gene,ref="NORM")
	return(data)
}
