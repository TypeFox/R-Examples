stratsample<-function(strata, counts){
	strata<-as.character(strata)
	n<-length(strata)
	rval <- integer(sum(counts))
	allrows<-1:n
	j<-0
	for(i in 1:length(counts))	{
		thisstrat<-names(counts)[i]
	   rval[j+(1:counts[i])]<-sample(allrows[strata==thisstrat],counts[i])
	   j<-j+counts[i]
	}
	rval
	}