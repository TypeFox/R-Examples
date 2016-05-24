MLD <-
function(samps){
	mld<-(1/length(samps))*sum(log(mean(samps)/samps))
	return(mld)
}
