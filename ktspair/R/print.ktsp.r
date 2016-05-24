print.ktsp <- function(x,...){
	ktspobj <- x
	cat(c("k-TSP object with:",ktspobj$k, "TSPs\n"))
	cat(c("Pair:\t\tTSP Score\t\tIndices\n"))
	for(i in 1:length(ktspobj$ktspscore)){
		cat(c("TSP",i,":","\t",round(ktspobj$ktspscore[i],2),"\t\t\t",ktspobj$index[i,],"\n"))
	}
}

