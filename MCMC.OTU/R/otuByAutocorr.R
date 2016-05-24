otuByAutocorr=function(model,data,ac.cut=0.1){
	ac=autocorr(model$Sol)
	ocorr=c();otus=c()
	for(o in unique(data$otu)){
		if (o=="summ") { next}
		spots=grep(paste(o,":",sep=""),dimnames(ac)[[2]])
		corrs=c()
		for (s in spots) {
			corrs=append(corrs,ac[2:length(dimnames(ac)[[1]]),s,s])
		}
		ocorr=append(ocorr,mean(corrs))
		otus=append(otus,o)
	}
	res=data.frame(cbind("autocorr"=ocorr))
	res=cbind("otu"=otus,res)
	goods=res[res$autocorr<ac.cut,"otu"]
	return(as.character(goods))
}
