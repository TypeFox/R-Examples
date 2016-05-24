##	Utility Functions


#####################################
#		  Skewness			#
#####################################




skewness<-function(x, finite=TRUE){
	n=length(x)
	S=(1/n)*sum((x-mean(x))^3)/(((1/n)*sum((x-mean(x))^2))^1.5)
	if(finite==FALSE){
		S=S
	}else{
		S=S*(sqrt(n*(n-1)))/(n-2)
		}
	return(S)	
}



#####################################
#		  Kurtosis			#
#####################################



kurtosis<-function(x, finite){
	n=length(x)
	K=(1/n)*sum((x-mean(x))^4)/(((1/n)*sum((x-mean(x))^2))^2) - 3
	if(finite==FALSE){
		K=K
		}
	else{
		K=((n-1)*((n+1)*K - 3*(n-1))/((n-2)*(n-3))) +3
		}
	return(K)	
}



#####################################
#		Amplitudes			#
#####################################


amps<-function(x){

	dens<-density(x)
	dd<-data.frame(dens$x,dens$y)

	maxima.ind<-which(diff(sign(diff(dd[,2])))==-2)+1
	minima.ind<-which(diff(sign(diff(dd[,2])))==2)+1
	colnames(dd)<-c("Location (x)","Amplitude (y)")
	peaks<-as.matrix(dd[maxima.ind,])
	antimode<-as.matrix(dd[minima.ind,])
	ls<-list(peaks, antimode)
	names(ls)<-c("Peaks", "Antimode")
	return(ls)
}



#####################################
#	N-th Highest Value		#
#####################################



nth_highest<-function(x,k=1){
	x<-as.vector(x)
	N<-length(x)
	sort(x,partial=N-k+1)[N-k+1]
}
