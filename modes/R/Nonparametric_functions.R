##	Nonparametric Functions

#####################################
#	Bimodality Coefficient		#
#####################################



bimodality_coefficient<-function(x, finite=TRUE,...){
	if(finite==TRUE){
		G=skewness(x,finite)
		sample.excess.kurtosis=kurtosis(x,finite)
		K=sample.excess.kurtosis
		n=length(x)
		B=((G^2)+1)/(K+ ((3*((n-1)^2))/((n-2)*(n-3))))
	}
	else{
		G=skewness(x,FALSE)
		K=kurtosis(x,FALSE)
		B=((G^2)+1)/(K)
	}
	return(B)
}


#####################################
#	Bimodality Amplitude		#
#####################################

bimodality_amplitude<-function(x,fig,...){

	dens<-density(x)
	dd<-data.frame(dens$x,dens$y)

	maxima.ind<-which(diff(sign(diff(dd[,2])))==-2)+1
	minima.ind<-which(diff(sign(diff(dd[,2])))==2)+1
	maxima<-dd[,1][maxima.ind]
	minima<-dd[,1][minima.ind]

	min.maxima<-min(dd[,2][maxima.ind])
	antinode<-dd[,2][minima.ind]

	B=(min.maxima-antinode)/min.maxima

	if(length(B)>1){
	B=min(abs(B))
	warning("Distribution is too bumpy. Consider bootstrapping
	or alternative method to get more data in the tails.")
	}

	if (fig=="TRUE"){
	plot(dens, main="Density plot with Peaks and Antimodes")
	abline(v=maxima[1],col="blue", lwd=2)
	abline(v=maxima[2],col="blue", lwd=2)
	abline(v=minima,col="darkgreen", lwd=2)
	text((max(maxima)+.06*diff(range(dd[,1]))),1.009*max(range(dd[,2])), "Max 1", col = "red") 
	text((min(maxima)-.06*diff(range(dd[,1]))),1.009*max(range(dd[,2])), "Max 2", col = "red") 
	text((minima+.06*diff(range(dd[,1]))),0, "Antimode", col = "red") 
	}

	return(B)	
}



#####################################
#	    Bimodality Ratio		#
#####################################


bimodality_ratio<-function(x,list=FALSE,...){

temp<-function(x){

	ls<-amps(x)
	dens<-density(x)
	dd<-data.frame(dens$x,dens$y)
	maxima.ind<-which(diff(sign(diff(dd[,2])))==-2)+1
	minima.ind<-which(diff(sign(diff(dd[,2])))==2)+1
	maxima<-dd[,1][maxima.ind]
	minima<-dd[,1][minima.ind]


	maxima<-dd[,1][maxima.ind]
	min.maxima<-min(dd[,2][maxima.ind])

	r.ind<-which(ls[[1]][,1] == max(ls[[1]][,1]))
	r<-(ls[[1]][r.ind,])
	l<-(ls[[1]][-r.ind,])

	R=r[2]/l[2]

	return(R)
			}
	if (list==FALSE){
	return(temp(x))

	}else{
	R<-lapply(x,temp)
	return(R)
	}
}


#####################################
#	 	    Modes			#
#####################################

modes<-function(data,type=1,digits="NULL",nmore="NULL"){
	if (type=="1"){
		tabs<-rle(sort(as.integer(data)))
		tabs<-t(as.matrix(cbind(tabs$values,tabs$lengths)))
		nmode<-as.matrix(rbind(tabs[1,][tabs[2,]==max(tabs[2,])],
			tabs[2,][tabs[2,]==max(tabs[2,])]))

}
	if (type=="2"){
		tabs<-rle(sort(signif(data,digits)))
		tabs<-t(as.matrix(cbind(tabs$values,tabs$lengths)))
		nmode<-as.matrix(rbind(tabs[1,][tabs[2,]==max(tabs[2,])],
			tabs[2,][tabs[2,]==max(tabs[2,])]))

}	
	if(!nmore=="NULL" & is.numeric(nmore) & type=="1"|type=="2"){
		for (i in 1:nmore){
		nmode<-cbind(nmode,as.matrix(rbind(tabs[1,][tabs[,2]==max(
		tabs[,2][tabs[,2]<nmode[2,i]])],tabs[2,][tabs[,2]==max(
		tabs[,2][tabs[,2]<nmode[2,i]])])))

			}
		}

	if (type=="3"){
	tabs<-rle(sort(as.vector(as.factor(data))))
	tabs<-t(as.matrix(cbind(tabs$values,tabs$lengths)))
	nmode<-as.matrix(rbind(tabs[1,][tabs[2,]==max(tabs[2,])],
		tabs[2,][tabs[2,]==max(tabs[2,])]))

	if(!nmore=="NULL" & is.numeric(nmore) & type=="3"){

	j<-dim(nmode)[2]
	for (i in j+1:nmore){
	tmp<-nth_highest(tabs[2,],i)
	tmp2<-tabs[1,][tabs[2,]==tmp]
	nmode<-cbind(nmode,rbind(tmp2,tmp))
	}}}

	if (any(nmode[2,]==1)==TRUE){warning ("A single observation
	 is being observed as a mode.
	Double check the class or inspect the data.
	Alternatively, you may have specified 'nmore' too many times 
	for this data.")}

	rownames(nmode)<-c("Value","Length")
	nmode<-as.matrix(nmode[,!duplicated(nmode[1,])])
return(nmode)
}


