countTextures <- function(Imgs,medpol=TRUE,...) {

I<-Imgs

if(medpol==TRUE){
	for(k in 1:length(I)){
		I[[k]] <- medpolish(I[[k]])$resid
	}
}

indu <- 1:(length(I)) 	# list of images yet to be classified

Iclass <- rep(0,length(I))

counter <- 0

while (length(I)>1) {

	res <- NULL
	counter <- counter + 1
	testcount<-0

	for(i in 1:(length(I)-1)) {
		res[i] <- compareImages(I[[1]], I[[i+1]],...)
	}

	testcount<-testcount+length(res)

	Iclass[indu[c(1,which(res==1)+1)]] <- counter

	I <- I[which(res==0)+1] 

	indu <- indu[which(res==0)+1]
}

Iclass[indu] <- counter +1

return(Iclass)

}
