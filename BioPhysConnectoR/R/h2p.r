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


 
get.entropy2p<-function(aln,bool=FALSE,gapchar="NOGAPCHAR",verbose=FALSE){
	gc<-paste(gapchar,collapse="|")
	l<-dim(aln)[2]
	hfunc<-function(j,i,aln,bool,gc){
		pp<-paste(aln[,i],aln[,j],sep="")
		ss<-summary(as.factor(pp),maxsum=dim(aln)[1])
		
		if(bool){
			ss[grep(gc,names(ss))]<-0 #so fällt es in der berechnung raus...
			}
		
		sum<-sum(ss)
		ss<-ss/sum
		
		ss<-ss[which(ss>0)]
		if(length(ss)==0){
			jH<-0
			if(verbose){
				cat("No pair without a gap character at the positions",i,j,".\n")
				}
			}
		jH<- -1*sum(ss*log2(ss))
		return(jH)
		}
		hhfunc<-function(i,aln,bool,gc){
			j<-i:l
			return(c(rep(0,i-1),apply(as.array(j),1,hfunc,i,aln,bool,gc)))
			#return(apply(as.array(j),1,hfunc,i,aln,bool,gc))
			#um eine vollständige matrix zu erhalten, obige zeilen aus- und die letzte einkommentieren
		}
		H<-apply(as.array(1:l),1,hhfunc,aln,bool,gc)
		return(H)
	} 
