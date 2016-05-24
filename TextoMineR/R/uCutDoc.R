uCutDoc <-
function(base, num.text, idiom="en", N=1000,alfa=0.05){

if (!is.null(num.text)) {
 if (is.character(num.text)) 
 num.text<- which(colnames(base) %in% num.text)
  if (is.numeric(num.text)) 
  num.text<- num.text
 if(length(num.text)==1)
  num.text<-num.text 
 if(length(num.text)>1){
   for(i in 1:length(num.text)){
    if(i==1)
    text1<-base[,num.text[1]]
    else text1<-paste(text1,base[,num.text[i]],sep=".") 
   }
    base[,(ncol(base)+1)]<-text1
    num.text<-ncol(base)
  }
}

Doc<-DocWordTable(base, num.text, idiom)
DocTermR<-Doc$DocTerm
res.ca<-CA(DocTermR[apply(DocTermR,1,sum)>0,], graph=FALSE)

Homo.groups<-uHomo.groups(res.ca$row$coord,N,alfa)
    GroupComposition<- Homo.groups$groups
  nb_groups<-length(Homo.groups$groups)
	bloque<-vector("list",nb_groups)
		for(i in 1:nb_groups){
		bloque[[i]]<- rep(i,length(Homo.groups$groups[[i]]))
	
	}
Group<-as.vector(unlist(bloque))
Group<-as.factor(Group)
Relim<-as.numeric(rownames(DocTermR[apply(DocTermR,1,sum)<=0,]))

if(length(Relim)==0)
	SentGroup<-cbind.data.frame(Group,base[,num.text]) 
else SentGroup<-cbind.data.frame(Group,base[-Relim,num.text])
   colnames(SentGroup)<-c("Group", "Sentence")
  PW<-matrix(nrow=length(Homo.groups$groups),ncol=ncol(DocTermR))
	for(i in 1:nb_groups){
      	elem<-vector()
     		elem<-(Homo.groups$groups[[i]])
		if (length(elem)==1) {
               	 PW[i,]<-DocTermR[Homo.groups$groups[[i]],]
                   } else {
     			 PW[i,]<-apply(DocTermR[Homo.groups$groups[[i]],],2,sum)
		}
        }

	rownames(PW)<-paste("P",1:nb_groups,sep="")
	colnames(PW)<-colnames(DocTermR)
	res = list(GroupWords=PW,SentGroup=SentGroup, GrpComposition=GroupComposition, Num.groups=nb_groups)
      return(res)
  }
