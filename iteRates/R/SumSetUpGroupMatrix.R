SumSetUpGroupMatrix <-
function(SUGM){
	final<-length(SUGM$groups)
	group<-SUGM$groups[[final]]
	VecBin<-SUGM$VecBin
	return(list(group=group,VecBin=VecBin))
	}

