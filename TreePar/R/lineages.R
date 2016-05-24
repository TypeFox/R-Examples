lineages <-
function(x,t){
	lin<-c(length(x)+1)
	if (length(t)>1) {
		for (j in 2:length(t)){
			lin<-c(lin, (1+length(which(x[]>t[j]))))
		}
	}
	lin
	}

