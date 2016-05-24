TxCharClust <-function(base, res, num.text){
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
res.hcpc<-res
n<-length(res.hcpc$desc.ind$para)
m<-length(res.hcpc$desc.ind$dist)
resPara<-vector(mode="list",n)
resdist<-vector(mode="list",m)
DocumentPara<-vector(mode="list",n)
DocumentDist<-vector(mode="list",m)
for(i in 1:n){
	IndPara<-as.data.frame(res.hcpc$desc.ind$para[i])
       resp<-base[which(rownames(base)%in%rownames(IndPara)),c(num.text)]
       resPara[[i]]<-resp
       IndPara$Resp<-resp
       DocumentPara[[i]]<-rownames(IndPara)
      IndPara$Document<-as.vector(rownames(IndPara))
	IndPara$cluster<-rep(i,length(res.hcpc$desc.ind$para[[i]]))          
	colnames(IndPara)[1]<-"X"
	if(i==1) tP<-IndPara    
	else tP<-rbind(tP,IndPara)
}

for(i in 1:m){
       IndDist<-as.data.frame(res.hcpc$desc.ind$dist[i])
	respD<-base[which(rownames(base)%in%rownames(IndDist)),c(num.text)]
       resdist[[i]]<-respD
       IndDist$RespD<-respD
       DocumentDist[[i]]<-rownames(IndDist)
      IndDist$Document<-as.vector(rownames(IndDist))
	IndDist$cluster<-rep(i,length(res.hcpc$desc.ind$dist[[i]]))          
	colnames(IndDist)[1]<-"X"
      IndDist<-IndDist[,c(2,3,4)]
	if(i==1) tD<-IndDist    
	else tD<-rbind(tD,IndDist)
}
    data.clust<-res.hcpc$data.clust
    desc.var <- descfreq(data.clust[, -which(sapply(data.clust, 
        is.factor))], data.clust[, ncol(data.clust)], proba = 0.05)
     CharWord1 = desc.var
     nClust<-length(CharWord1)
     CharClust<-vector(mode="list",nClust)
     names(CharClust)<-paste("clust", 1:nClust, sep = "")
     OverCharWord = vector(mode = "list", length = length(CharWord1))
     names(OverCharWord) <- names(CharWord1)
     InfraCharWord = vector(mode = "list", length = length(CharWord1))
     names(InfraCharWord) <- names(CharWord1)
 for (i in 1:length(CharWord1)){
        bWord<-as.data.frame(CharWord1[[i]])
        over<-subset(bWord, bWord$v.test>0)
        infra<-subset(bWord, bWord$v.test<0)
        OverCharWord[[i]] <- row.names(over)
        InfraCharWord[[i]] <- row.names(infra)
           Para<-resPara[[i]]
           DocP<-DocumentPara[[i]]
           Dist<-resdist[[i]]
           DocD<-DocumentDist[[i]]
       X <-list(OverCharWord[[i]],InfraCharWord[[i]],DocP,Para,DocD,Dist)              
          names(X)<-c("Over_represented_word","Infra_represented_word",
              "Close_to_centroid_documents_(label)","Close_to_centroid_documents_(content)",
              "Far_from_other_clusters_document_(label)" ,"Far_from_other_clusters_document_(content)")
          CharClust[[i]]<-X
         }
    res=list(Char_word_doc=CharClust,CharWord_details=CharWord1)
 return(res) 
}
