CharDocWord <-
function(base=NULL,DocTerm,num.text=NULL,num.agg=NULL,maxDocs=10,list.details=FALSE)
{
if(!is.null(base)){
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
#### Fonctions auxiliaires
CALVAL<-function(x,y) 
{
  x<-as.matrix(x)
  y<-as.matrix(y)
  res<-sum((as.matrix(x))*(as.matrix(y)))
  return(res)
 }
# arguments:  
colgrou<-num.agg
ntext<-num.text
#  calculus of the v.test of the words for every group
res.descfreq<-descfreq(DocTerm, by.quali = base[,colgrou], proba=1)

        CharWord=vector(mode="list",length=length( res.descfreq))
        names(CharWord)<-names( res.descfreq)
        for(i in 1:length(res.descfreq)){
        bWord <- as.data.frame(res.descfreq[[i]])
        over <- subset(bWord, bWord$v.test > 0)
        Over<- row.names(over)
        infra <- subset(bWord, bWord$v.test < 0)
        Infra<- row.names(infra)
        OverInfra<-list(Over,Infra)
        names(OverInfra)<-c("Over_represented_word","Infra_represented_word")
        CharWord[[i]]<-OverInfra
        }
        CharWords<-CharWord

nmot<-length(row.names(res.descfreq[[1]]))

vlev<-levels(base[,colgrou]) 
nlev<-length(vlev)
motsval<-matrix(nrow=nmot,ncol=nlev)
 

 for (igru in 1:nlev)
    {
      matess<-res.descfreq[[igru]][,6]
      matessbis<-matess[order(row.names(res.descfreq[[igru]]))]
      motsval[,igru]<-matessbis
    }

  row.names(motsval)<-colnames(DocTerm)
  colnames(motsval)<-vlev

  tmotsval<-t(motsval)
  colnames(tmotsval)<-row.names(motsval)
  row.names(tmotsval)<-colnames(motsval)

   lisresult<-list()
   for (igru in 1:nlev)
   {
     tabgru<-data.frame()
     resp<-which(base[,colgrou]==vlev[igru])
     nrepgcurs<-length(resp)
     ntrep<-min(maxDocs,nrepgcurs)
	 
     if (ntrep > 1)
        {
        DocTermcurs<-matrix(nrow=nrepgcurs,ncol=nmot)
        print(dim(DocTermcurs))
        DocTermcurs<-(DocTerm[which(base[,colgrou]==vlev[igru]),])
        DocTermcurs<-as.matrix(DocTermcurs)
       y<-as.matrix(tmotsval[igru,])
       a<- apply(DocTermcurs,MARGIN=1,FUN=CALVAL,y)
       b<-apply(DocTermcurs,MARGIN=1,FUN=sum)
       repvaltest<-a
       repvaltest[which(b>0)]<-a[which(b>0)]/b[which(b>0)]  
       ordrep<-order(repvaltest,decreasing="TRUE")  
       for (i in 1:ntrep)
         {
         tabgru[i,1]<-row.names(base)[resp[ordrep[i]]]
         tabgru[i,2]<-repvaltest[ordrep[i]]
         tabgru[i,3]<-base[resp[ordrep[i]],ntext]
         } 
      } 
     if (ntrep == 1)
      {
         aa<-0
         bb<-0
         for (imot in 1:nmot)
         {
           aa<- aa + DocTerm[resp,imot]*motsval [imot,igru]
           bb<- bb + DocTerm[resp,imot]
         }
          rvaltest <-0
          if (bb>0) rvaltest<-aa/bb
          tabgru[1,1]<-row.names(base)[resp]
          tabgru[1,2]<-rvaltest
          tabgru[1,3]<-base[resp,ntext]
       }
      colnames(tabgru)<-c("n.rep", "val.moy", "document")
      lisresult[[igru]]<-tabgru
      }
    names(lisresult)<-levels(base[,colgrou])
    CharDoc<-lisresult
     if(list.details)
      CharWord<-res.descfreq
	 else CharWord<-CharWords
 }else{
     CharDoc<-NULL
     CharWord1<-descfreq(DocTerm, proba=0.05)
    if(list.details)
     CharWord<-CharWord1
     else{
      CharWord=vector(mode="list",length=length( CharWord1))
        names(CharWord)<-names( CharWord1)
        for(i in 1:length(CharWord1)){
        bWord <- as.data.frame(CharWord1[[i]])
        over <- subset(bWord, bWord$v.test > 0)
        Over<- row.names(over)
        infra <- subset(bWord, bWord$v.test < 0)
        Infra<- row.names(infra)
        OverInfra<-list(Over,Infra)
        names(OverInfra)<-c("Over_represented_word","Infra_represented_word")
        CharWord[[i]]<-OverInfra
        }
        CharWord<-CharWord
	}	
    }
     if(!is.null(CharDoc)){
      CharDocWord1<-vector(mode="list",length(CharDoc))
      names(CharDocWord1)<-names(CharDoc)
      CharDoc1<-vector(mode="list",length(CharDoc))
      names(CharDoc1)<-names(CharDoc)
       for(i in 1:length(CharDoc)){
       ChDoc<-as.data.frame(CharDoc[[i]])
       DocLb<-ChDoc[,1]
       DocCont<-ChDoc[,3]
       labelCont<-list(DocLb,DocCont)
       names(labelCont)<-c("Document_label", "Document_content")
       CharDoc1[[i]]<-labelCont
        X <-list(CharWord[[i]],CharDoc1[[i]])
            names(X)<-c("Words","Documents")
               CharDocWord1[[i]]<-X
         }
        CharDoc<-CharDoc1
        }else CharDocWord1<-NULL
     

    res<-list(CharDoc=CharDoc,CharWord=CharWord, CharDocWord=CharDocWord1)
    return(res)
  
 }
