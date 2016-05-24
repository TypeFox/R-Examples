`logchoose` <-
function(n,D,last=FALSE) {
    tmp <- c(0,cumsum(log((n+1-(1:D))/(1:D))))
    if(last) tmp <- tmp[length(tmp)]
    return(tmp)
}
            
calcNL<-function(y,BL,right) {
	n<-length(y)
	bn<-length(BL)
	NL<-rep(0,bn)
	if (right) {
		k<-sum(y<=BL[1])
		NL[1]<-k
		ind<-2
		while (k<n) {
			k<-k+1
			rbd<-BL[ind]
			while (y[k]>rbd) {
				NL[ind]<-k-1
				ind<-ind+1
				rbd<-BL[ind]
				}
			NL[ind]<-k
		}
	}
	else {
		k<-sum(y<BL[1])
		NL[1]<-k
		ind<-2
		while (k<n) {
			k<-k+1
			rbd<-BL[ind]
			while (y[k]>=rbd) {
				NL[ind]<-k-1
				ind<-ind+1
				rbd<-BL[ind]
				}
			NL[ind]<-k
		}
	}
	if (ind<bn) NL[ind:bn]<-n
	NL
}
