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


 
get.entropy<-function(aln,bool=FALSE,gapchar="NOGAPCHAR",verbose=FALSE){
	
	aln<-toupper(aln)
	
	entropy<-function(i,aln,bool,gc){
		ss<-summary(as.factor(aln[,i]))
		
		if(bool){
			ss[grep(gc,names(ss))]<-0
			}
		
		sum<-sum(ss)
		ss<-ss/sum
		
		ss<-ss[which(ss>0)]
		
		if(length(ss)==0){
			ret<-0
			if(verbose){
				cat("No Non-Gap Character.\n")
				}
			}else{
				ret<- -1*sum(ss*log2(ss))
				}
				
		return(ret)
		}
		
	gc<-paste(gapchar,collapse="|")
	l<-ncol(aln)
	indizes<-1:l
	H<-apply(as.array(indizes),1,entropy,aln,bool,gc)
	
	return(H)
	}
 
