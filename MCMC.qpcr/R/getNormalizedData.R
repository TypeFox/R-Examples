getNormalizedData=function(model,data,controls=NULL){
	dd=data
	if (!is.null(controls)) {
		dd=dd[!(dd$gene %in% controls),]
	}
	prr=predict(model,type="terms",marginal=~sample)
	dd$pred=prr[1:nrow(dd)]
	predicted=c()
	for (g in unique(dd$gene)) {
			if (g %in% controls) { next }
			gs=c();coo=c()
			for(s in unique(dd$sample)) {
				gs=append(gs,dd$pred[dd$gene==g & dd$sample==s][1]/log(2))
				coo=rbind(coo,dd[dd$gene==g & dd$sample==s,][1,4:(ncol(dd)-1)])
			}		
			predicted=cbind(predicted,gs)
	}
	predicted=data.frame(predicted)
	conditions=data.frame(coo)
	names(predicted)=unique(dd$gene)
	row.names(predicted)= unique(dd$sample)
	row.names(conditions)= unique(dd$sample)
	return(list("normData"=predicted,"conditions"=conditions))
}
