DocWordTable <-
function(base,num.text,idiom="en",lminword=1,lower=TRUE,Remov.Number=TRUE)
{
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
  corpus <- Corpus(DataframeSource(base[num.text]), readerControl=list(language=idiom))
  dtmCorpus <- corpus
   filt="(['?]|[[:punct:]]|[[:space:]]|[[:cntrl:]])+"
   	dtmCorpus <- tm_map(dtmCorpus, content_transformer(function(x) gsub(filt, " ", x)))
   if (Remov.Number == TRUE)dtmCorpus <- tm_map(dtmCorpus, removeNumbers)
   dtm <- DocumentTermMatrix(dtmCorpus, control=list(tolower=lower,
   wordLengths=c(lminword, Inf)))
   rownames(dtm)<-rownames(base)
   
   dtmM<-as.matrix(dtm)
   dimdtm<-dim(dtmM)
   Ndoc<-dimdtm[1]
   Nword<-dimdtm[2]
   Nfreqword<-apply(dtmM,MARGIN=2,FUN=sum)
   Nlength<-sum(Nfreqword)
   dtmA<-dtmM
   dtmA[dtmM[,]>0]<-1
   Ndocword<-apply(dtmA,MARGIN=2,FUN=sum)
   Nfreqwords<-as.data.frame(Nfreqword)
   Ndocwords<-as.data.frame(Ndocword)
   Table<-cbind(Nfreqwords,Ndocwords)
   TFreq<-with(Table,Table[order(Nfreqwords,Ndocwords,decreasing=TRUE),])
   colnames(TFreq)<-c("Frequency","N.Documents")
   res<-list(DocTerm=dtmM,Ndoc=Ndoc,Nlength=Nlength,Nword=Nword,Tfreq=TFreq, Nfreqword=Nfreqword, Ndocword=Ndocword )
   class(res)<-c("DocWordTable","list")
return(res)
 }
