outlierSamples=function(model,data,z.cutoff=-2) {
	allmean=apply(model$Sol,2,mean)
	seff=allmean[grep("sample.",names(allmean))][c(1:length(levels(data$sample)))]
	seff.z=(seff-mean(seff))/sd(seff)
	return(sub("sample.","",names(which(seff.z<z.cutoff))))
}
