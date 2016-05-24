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


build.interact<-function(cseq,mj1,mj2=mj1,d,alpha=82){

	if(!is.numeric(cseq)){
		if(is.character(cseq)){
			cseq<-aa2num(cseq,0,verbose=FALSE)
			}
			else{
				stop("Wrong Argument Type.\n")
				break
				}
		}

	m<-20;
	n<-length(cseq)
	nd<-length(d)
	interaction.mat<-matrix(data=0,ncol=n,nrow=n)
	out<-.C("buildInteract",cseq=as.integer(cseq),n=as.integer(n),mj1=as.double(mj1),			mj2=as.double(mj2),m=as.integer(m),d=as.integer(d),nd=as.integer(nd),interaction.mat=as.double(interaction.mat),PACKAGE="BioPhysConnectoR")

	interaction.mat<-matrix(data=out$interaction.mat,ncol=sum(d))
	chains<-cumsum(d)
	nchains<-length(chains)
	covind<-1:(n-1)
	if(nchains>1){
	  covind<-covind[-chains[-nchains]]
	}
	diag(interaction.mat[,-1])[covind]<-diag(interaction.mat[-1,])[covind]<-alpha
	diag(interaction.mat)<-0
	return(interaction.mat)
	}
