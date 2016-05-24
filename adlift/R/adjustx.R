"adjustx" <-
function(x,f,type="mean"){

#adjusts the grid (x) values and function values if necessary (due to repeat grid values)

#possible types: mean and jitter

groups<-list()

if (type=="jitter"){
	sepx<-NULL
	sepf<-f
	if (length(x)==length(unique(x))){	#do nothing
		sepx<-x
	}
	else{
		for (i in 1:length(unique(x))){
			q<-which(x==unique(x)[i])
			groups[[i]]<-q
			if (length(q)==1){ 		#don't change x 
				sepx[q]<-x[q]
			}
			else{
				eps<-1+(10^(-7))*(0:(length(q)-1))
				sepx[q]<-x[q]*eps
			}
		} 
	}
}

if (type=="mean"){
	sepx<-unique(x)
	sepf<-matrix(0,1,length(sepx))

	for (i in 1:length(sepx)){
		q<-which(x==sepx[i])
		groups[[i]]<-q
		sepf[i]<-mean(f[q])
	}
}

return(list(sepx=sepx, sepf=sepf,groups=groups))

}
