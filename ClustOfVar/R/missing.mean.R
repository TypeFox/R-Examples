missing.mean <-
function(C1){
	moy <- mean(C1,na.rm=T)
	ind <- which(is.na(C1)==T)
	if(length(ind)>=1){C1[ind]<-moy
	}
	return(C1)
}

