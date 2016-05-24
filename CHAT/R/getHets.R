getHets <- function(baf,paired.baf=NULL,thr=0.18){
	
	#Retrieve heterozygous SNP from BAF file
	#if paired.baf is given, the heterozygous markers are estimated using the normal tissue
	
	tag<-c()
	window.size<-100 #the size of sliding window
	k<-1
	while(k+window.size<=length(baf)){
		if(!is.null(paired.baf)){ #case with paired samples
			tmp<-paired.baf[k:(k+window.size-1)]
			hets<-k+which(tmp>=thr&tmp<=(1-thr))-1
			tag<-c(tag,sort(hets))
			k<-k+window.size
			next
		}
		tmp<-baf[k:(k+window.size-1)]
		hets<-k+which(tmp>=thr&tmp<=(1-thr))-1
		if(length(hets)<=window.size/10){	#region of LOH
			hets<-k+round(runif(round(window.size/3),1,window.size))-1
		}
		tag<-c(tag,sort(hets))
		k<-k+window.size
	}
	return(tag)
}