#############################################
#   This code is subject to the license as stated in DESCRIPTION.
#   Using this software implies acceptance of the license terms:
#    - GPL 2
#
#   (C) by F. Hoffgaard, P. Weil, and K. Hamacher in 2009.
#
#  keul(AT)bio.tu-darmstadt.de
#
#
#  http://www.kay-hamacher.de
#############################################



freq1p<-function(aln,i=NULL){
	freq1p.func<-function(i,aln,mat=NULL){
		l<-length(aln[,i])
		freq<-summary(as.factor(aln[,i]),maxsum=l)
		if(!is.null(mat)){
			mat[names(freq),i]<-freq
			return(mat[,i])
		}else{
			return(freq)
		}
	}
	if(!is.null(i)){
		return(freq1p.func(i=i,aln=aln))
	}else{
		lett<-unique(as.vector(aln))
		n<-dim(aln)[2]
		mat<-matrix(0,nrow=length(lett),ncol=n)
		rownames(mat)<-lett
		res<-apply(as.array(1:n),1,freq1p.func,aln,mat)
		return(res)
	}
}
freq2p<-function(i,aln,j2=NULL,lett=NULL,cluster=NULL){
	freqs<-function(ind,aln,i,j,mat){
		lp<-length(aln[,i])
		nb<-paste(aln[,i],aln[,j[ind]],sep="")
		freq<-summary(as.factor(nb),maxsum=lp)
		#print(freq)
		#print(j[ind])
		mat[names(freq),ind]<-freq
		return(mat[,ind])
	}
	if(is.null(lett)){
			lett<-unique(as.vector(aln))
	}
	nlett<-length(lett)
	l<-dim(aln)[2]
        cores<-1
	if(!is.null(cluster)){
        	clusterEvalQ(cluster,library(BioPhysConnectoR))
		cores<-length(cluster)
	}
	func<-function(bb,bbvec){
		return(paste(bb,bbvec,sep=""))
	}

	if(is.null(j2)){
		j2<-i:l
		ln<-length(j2)
	}else{
		ln<-length(j2)
	}
	ind<- 1:ln
	ok<-apply(as.array(lett),1,func,lett)
	mat<-matrix(0,nlett^2,ln)
	rownames(mat)<-ok
	ind<-as.list(ind)
	### new code starts here ####
	##callingfunction for parLapply / lapply
	fun<-function(ind,aln,i,j2,mat){
		return(lapply(ind,freqs,aln,i,j2,mat))
	}
	##list formatting to get a list of lists with length cores
	lind<-length(ind)
  if(cores > lind){
    cores<-lind
  }
	overlap<-lind%%cores
	ll<-floor(lind/cores)
	nind<-list()
	it<-0
	for(k in 1:cores){
		it<-it+1
		if(overlap>=k){
			itn<-it+ll
		}else{
			itn<-it+ll-1
		}
		nind[[k]]<-ind[it:itn]
		it<-itn
	}
	if(!is.null(cluster)){
		bb<-parLapply(cluster,nind,fun,aln,i,j2,mat)
	}else{
		bb<-lapply(nind,fun,aln,i,j2,mat)
	}
	lb<-lbb<-cnam<-vector()
	for(k in 1:cores){
		lb<-c(lb,length(bb[[k]]))
		for (l in 1:lb[k]){
			lbb<-c(lbb,length(bb[[k]][[l]]))
			cnam<-rbind(cnam,names(bb[[k]][[l]]))
		}
	}
	t1<-unique(lb)
	t2<-unique(lbb)
	t3<-unique(cnam)
	if(length(t2)==1 && dim(t3)[2]==t2){
		resbb<-matrix(unlist(bb),ncol=lbb,byrow=TRUE)
		resbb<-t(resbb)
		rownames(resbb)<-t3
	## original code continues here ####
		return(resbb)
	}
}
