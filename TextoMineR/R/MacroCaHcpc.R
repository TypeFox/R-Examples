MacroCaHcpc<-function(base,num.text, idiom="en",nb.clust =-1, Fmin=5, Dmin=5, Fmax=NULL,
equivalence=NULL, stop.word.user=NULL, stop.word.tm=FALSE,
  lmd=3, lmk=3,  ncp=5, row.sup = NULL, col.sup = NULL, quanti.sup=NULL,
   quali.sup = NULL, axes = c(1,2)) {

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

###Corpus
Corpus<-DocWordTable(base,num.text,idiom,lminword=1,Remov.Number=TRUE,lower=TRUE)
#summary(Corpus)
DocTerm<-Corpus$DocTerm   

# Selection of the words depending from two threshold minimum global frequency minimum nr of documents 
res.TxCA<-TxCA(DocVar=base,DocTerm, idiom=idiom, Fmin=Fmin,Dmin=Dmin,Fmax=Fmax,
     equivalence=equivalence,stop.word.user=stop.word.user,stop.word.tm=stop.word.tm, ncp=ncp, 
              row.sup=row.sup,col.sup= col.sup, quanti.sup,quali.sup=quali.sup,graph=FALSE)
res.ca<-res.TxCA$res.ca

Ch<-TRUE
i=1
while(Ch & i<=nrow(res.ca$eig)){
  if(res.ca$eig[i,1]<mean(res.ca$eig[,1])){
   ncp<-i-1
   Ch<-FALSE
  }
  else i<-i+1
}
 if(ncp==1){ 
 ncp=2
}

res.TxCA<-TxCA(DocVar=base,DocTerm, idiom=idiom, Fmin=Fmin,Dmin=Dmin,Fmax=Fmax,
     equivalence=equivalence,stop.word.user=stop.word.user,stop.word.tm=stop.word.tm, ncp=ncp, 
              row.sup=row.sup,col.sup= col.sup, quanti.sup,quali.sup=quali.sup,graph=FALSE)
#summary(res.TxCA)
DocTermR<-res.TxCA$DocTermR
Table<-res.TxCA$Table
res.ca<-res.TxCA$res.ca

# Classification 
res.hcpc<-HCPC(res.ca,cluster.CA="rows", nb.clust =nb.clust, order=TRUE)

res.TxCharClust<-TxCharClust(base,res.hcpc,num.text)
res.TxCharClust
plot(res.ca,invisible="row")
dev.new()
plot(res.ca,invisible="col")

res<-list(Corpus=Corpus, res.TxCA=res.TxCA, res.hcpc=res.hcpc, res.TxCharClust=res.TxCharClust, ncp=ncp)
class(res)<-c("MacroCaHcpc","list")
return(res)
}



