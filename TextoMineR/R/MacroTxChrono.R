MacroTxChrono<-function(base, num.text, divide=TRUE, SentLength=100,  idiom ="en", Fmin=5, Dmin=1, 
Fmax=NULL, equivalence=NULL, stop.word.user=NULL,stop.word.tm=FALSE, lmk=3,lmd=3, ncp=5, 
 row.sup = NULL, col.sup = NULL, quanti.sup=NULL, quali.sup = NULL, graph = TRUE, axes = c(1,2),
 N=5000, alfa=0.15, CorChronoDim=0.10, HierWords = TRUE, SegRep=FALSE){
if(divide){
#Divide the text into sentences homogeneous
 phrase<- uSentences(base, num.text, SentLength) 
#Phrases are grouped into homogeneous parts.
         
      Homo.Groups<-uCutDoc(phrase,idiom, num.text=1, N, alfa)
      base<-Homo.Groups$SentGroup
      
}else{
#Phrases are grouped into homogeneous parts.
     Homo.Groups<-uCutDoc(base,idiom,num.text,N,alfa)
     base<-Homo.Groups$SentGroup
     phrase=NULL
}
  res.Homo.Groups<-list(composition=Homo.Groups$GrpComposition,nb.groups=Homo.Groups$Num.groups)

#Contextual variables
   DocVar<-DocVarTable(base, VarSel=1)
#DocumentosXword Table
   Corpus<-DocWordTable(base,num.text=2,idiom) 
   DocTerm<-Corpus$DocTerm

#CA Show trajectory of the parties

#CA lexica-Agreggated table
res.TxCA1<-TxCA(DocVar,DocTerm,num.agg=1,idiom=idiom,Fmin=Fmin, Dmin=Dmin,Fmax=Fmax,
  equivalence=equivalence,stop.word.user=stop.word.user,stop.word.tm=stop.word.tm, ncp=ncp, 
  row.sup=row.sup,col.sup= col.sup,graph =FALSE)
 Tagreg<-res.TxCA1$Tagreg
#Chrono variable
DocVarAgreg<-as.data.frame(rep(1:nrow(Tagreg)))
rownames(DocVarAgreg)<-rownames(Tagreg)
colnames(DocVarAgreg)<-"Chrono"

#define the number of axes
res.TxCA2<-TxCA(DocVarAgreg,Tagreg,idiom=idiom,Fmin=Fmin,Dmin=1,Fmax=Fmax,quanti.sup=1, graph =FALSE)
res.ca<-res.TxCA2$res.ca
Nc<-res.ca$quanti.sup$cos2
rownames(Nc)<-"Correlation"

if(Nc[1]< 0.05){ 
 print("Chronology is not obviously related to the first axis")
 print(Nc)
 stop
 }
Ch<-TRUE
i=1
while(Ch & i<=length(Nc)){
  if(Nc[i]< CorChronoDim){
   ncp<-i-1
   Ch<-FALSE
  }
  else i<-i+1
}
 if(ncp==1){ 
 ncp=2
}
#Trayectory
res.TxCA<-TxCA(DocVar,DocTerm, num.agg=1, idiom=idiom,Fmin=Fmin, Dmin=Dmin,Fmax=Fmax,
     equivalence=equivalence,stop.word.user=stop.word.user,stop.word.tm=stop.word.tm, ncp=ncp, 
     row.sup=row.sup,col.sup= col.sup, quanti.sup=quanti.sup,quali.sup=quali.sup, graph =FALSE)
Table<-res.TxCA$Table
res.ca<-res.TxCA$res.ca

dev.new()
res.meta1<-META.CA(res.ca,naxes=ncp,axe.x=axes[1],axe.y=axes[2],lmd=Inf,lmk, main="CA words")
dev.new()
plot(res.ca,invisible="col")
lines(res.ca$row$coord[,1],res.ca$row$coord[,2],lwd=2,col="black")
dev.new()
sel<-which( (res.ca$col$contrib[,1]> 3*mean(res.ca$col$contrib[,1]))
            | (res.ca$col$contrib[,2]>3*mean(res.ca$col$contrib[,2])) )
plot(res.ca,choix="CA",selectCol=sel,unselect=1)
lines(res.ca$row$coord[,1],res.ca$row$coord[,2],lwd=2,col="black")

# Hierarchical clustering
res.ca<-res.TxCA$res.ca
res.chcpc<-TxCHCPC(res.ca,cluster.CA = "rows", nb.clust = -1) 
    if (HierWords) {
        HierWord <- HierarchWords(res.ca, Table)
    } else HierWord = NULL
#Hegments in the hierarchy
if(SegRep){
	res.segment<-SegmentsRep(base, num.text=2, nxlon=5, nfreq=2)
	Tab.SegR<-res.segment$tab.seg
        num.agg <-1
        agg <- as.factor(base[, num.agg])
        dis.X <- tab.disjonctif(agg)
        Tagreg.LexSeg <- t(Tab.SegR) %*% dis.X
        TagSeg<- t(Tagreg.LexSeg)
	res.ca<-res.TxCA$res.ca
	HierSegment<-HierarchWords(res.ca, TagSeg)
}else HierSegment=NULL

#Vocabulary index
res.VocIndex<-VocIndex(base, num.text=2, Fmin=Fmin)

res<-list(SentenceList=base,Homo.Groups=res.Homo.Groups,Corpus=Corpus,Correlation=Nc,
       res.TxCA=res.TxCA,res.chcpc=res.chcpc,HierWord=HierWord,HierSegment=HierSegment,VocIndex=res.VocIndex)
class(res)<-c("MacroTxChrono","list")
return(res)
}









