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


 
get.contact.list<-function(cm,d=NULL,single=TRUE,val=1){
	cl<-which(cm==val,arr.ind=TRUE)
	if(single){
		cl<-cl[which(cl[,1]<cl[,2]),]
	}
	if(!is.null(d)){
			oo<-which((cl[,1]+1==cl[,2])&(!cl[,1]%in%d))
			cl<-cl[-oo,]
			if(!single){
				oo<-which((cl[,1]-1==cl[,2])&(!cl[,1]%in%d))
				cl<-cl[-oo,]
			}
		}
	cl<-as.list(unname(as.data.frame(t(cl))))
	return(cl)
	} 
