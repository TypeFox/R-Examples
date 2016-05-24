getNormalizedOTUdata <-
function(model,data,log.base=10){
#model=mm;data=gss
	dd=data
	dd$pred=predict(model,type="terms",marginal=~sample)
	predicted=c()
	for (g in unique(dd$otu)) {
			gs=c();coo=c()
			for(s in unique(dd$sample)) {
				gs=append(gs,dd$pred[dd$otu==g & dd$sample==s][1]/log(log.base))
				coo=rbind(coo,dd[dd$otu==g & dd$sample==s,][1,3:(ncol(dd)-1)])
			}		
			predicted=cbind(predicted,gs)
	}
	
	predicted=data.frame(predicted)
	conditions=data.frame(coo)
	names(predicted)=unique(dd$otu)

	if("summ" %in% unique(dd$otu)) {
		norm=predicted[,"summ"]-mean(predicted[,"summ"])
		predicted[,"summ"]=NULL
		predicted=predicted-norm
	}
	
	row.names(predicted)= unique(dd$sample)
	row.names(conditions)= unique(dd$sample)
	return(list("normData"=predicted,"conditions"=conditions))
}
