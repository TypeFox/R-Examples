entropy.corrected <-function(H,ErrorHX,HXmax){
	
	H[H>(HXmax-ErrorHX)]<-HXmax
	H[H<ErrorHX]<-0
	return(H)
	}

