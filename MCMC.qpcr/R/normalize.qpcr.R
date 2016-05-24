normalize.qpcr=function(data,controls,center=T) {
	dln=c()
	for(s in levels(data$sample)){ 
		sub=subset(data,sample==s)
		cm=0
		for(g in controls) {
			cm=cm+mean(sub[sub$gene==g,"count"])
		}
		nor=round(cm/length(controls),2)
		sub$count=sub$count-nor
		dln=rbind(dln,sub)
	}
	# centering
	if (center==T){
		for(ge in levels(data$gene)){ 
			dln[dln$gene==ge,"count"]=dln[dln$gene==ge,"count"]-round(mean(dln[dln$gene==ge,"count"]),2)
		}
	}
	return(dln)
}
