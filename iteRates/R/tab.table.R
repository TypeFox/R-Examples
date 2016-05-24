tab.table <-
function(cst){
	sumcomp<-list()
	for(i in 1:length(cst)){
		ul<-unlist(cst[[i]])
		nms<-unique(names(ul))
		k<-as.numeric(ul[names(ul)==nms[1]])
		groups<-as.character(ul[names(ul)==nms[2]])
		P1<-as.numeric(ul[names(ul)==nms[3]])
		P2<-as.numeric(ul[names(ul)==nms[4]])
		lik<-as.numeric(ul[names(ul)==nms[5]])
		modID<-as.numeric(ul[names(ul)==nms[7]])#change from 6 to 7
		tn<-as.numeric(ul[names(ul)==nms[8]])#change from 7 to8
		sumcomp[[i]]<-list(k=k,groups=groups,P1=P1,P2=P2,lik=lik,modID=modID,tn=tn)
		}
	Sumcomp<-list()
	for(i in 1:length(sumcomp)){
		ka<-sumcomp[[i]]$k
		groupsa<-sumcomp[[i]]$groups
		P1a<-matrix(sumcomp[[i]]$P1,nrow=(length(sumcomp[[i]]$groups)),byrow=TRUE)
		P2a<-matrix(sumcomp[[i]]$P2,nrow=(length(sumcomp[[i]]$groups)),byrow=TRUE)
		lika<-matrix(sumcomp[[i]]$lik,nrow=(length(sumcomp[[i]]$groups)),byrow=TRUE)
		modIDa<-matrix(sumcomp[[i]]$modID,nrow=(length(sumcomp[[i]]$groups)),byrow=TRUE)
		tna<-matrix(sumcomp[[i]]$tn,nrow=(length(sumcomp[[i]]$groups)),byrow=TRUE)
		Sumcomp[[i]]<-list(k=ka,groups=groupsa,P1=P1a,P2=P2a,lik=lika,modID=modIDa,tn=tna)
		}
		
tempcolnam<-as.character(seq(1,2+(length(Sumcomp)*5)))

fintab.a<-list()
	for (i in 1:length(Sumcomp)){	
				rawsumtab<-data.frame(Sumcomp[[i]]$k,Sumcomp[[i]]$groups,Sumcomp[[i]]$P1,Sumcomp[[i]]$P2,Sumcomp[[i]]$lik,Sumcomp[[i]]$modID,Sumcomp[[i]]$tn)
				tempres<-rawsumtab[,-c(1,2)]
				colZm<-dim(tempres)[2]
				sets<-colZm/5
				colorder<-numeric()
				for (k in 1:sets){
					temp<-seq(k,colZm,by=sets)
					colorder<-c(colorder,temp)
					}
				colreorder<-c(1,2,colorder+2)
				ntempres<-rawsumtab[,colreorder]
						coltoadd<-length(tempcolnam)-dim(ntempres)[2]
						rowtoadd<-dim(ntempres)[1]
						filler<-array(dim=c(rowtoadd,coltoadd))
				ntempres<-cbind(ntempres,filler)
				colnames(ntempres)<-tempcolnam
			fintab.a[[i]]<-ntempres
		}
		
fintab.b<-fintab.a[[1]]		
	if ((length(fintab.a))>1){
		for (i in 2:length(fintab.a)){
			fintab.b<-rbind(fintab.b,fintab.a[[i]])
						}
					}	
	return(fintab.b)					
	}

