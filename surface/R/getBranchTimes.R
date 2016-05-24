getBranchTimes <-
function(h){
	h<-as(h,"data.frame")
	brtimes<-rep(NA,dim(h)[1])
	for(i in 2:length(brtimes)){
		brtimes[i]<-h$times[as.numeric(as.character(h$ancestors[i]))]		
		}
brtimes
}
