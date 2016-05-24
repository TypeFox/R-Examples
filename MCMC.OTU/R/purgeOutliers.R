purgeOutliers=function(data,count.columns,sampleZcut=(-2.5),otu.cut=0.001,zero.cut=0) {
# samples
	lsums=log(apply(data[,count.columns],1,sum))
	zscore=(lsums-mean(lsums))/sd(lsums)
	badz=zscore[zscore<sampleZcut]
	sams=grep("[sS]ample",names(data))
	names(data)[sams]="sample"
	badSam=data$sample[zscore<sampleZcut]
	print(paste("samples with counts below z-score",sampleZcut,":"))
	print(as.character(badSam))
	print("zscores:")
	print(badz)
	data=data[!(data$sample %in% badSam),]
# OTU sum counts
	goodOTU=names(data)[count.columns]
	if (otu.cut>0) {	
		tot=sum(data[,count.columns])
		otuCut=tot*otu.cut
		sums=apply(data[,count.columns],2,sum)
		badcount=sums[sums<otuCut]
		goodOTU=names(data)[count.columns][sums>=otuCut]
		print(paste("OTUs passing frequency cutoff ",otu.cut,":",length(goodOTU)))
	}
#	print(goodOTU)
# zeroes
	if (zero.cut>0) {	
		zf=apply(data[,goodOTU],2,function(x){ return(sum(x>0)/length(x)) } )
		print(paste("OTUs with counts in",zero.cut,"of samples:"))
		print(table(zf>=zero.cut))
		goodOTU=goodOTU[zf>=zero.cut]
	}
	conds=1:length(names(data))
	conds=conds[!(conds %in% count.columns)]
	cdat=data[,conds]
	cotu=data[,names(data) %in% goodOTU]
	return(cbind(cdat,cotu))
}
