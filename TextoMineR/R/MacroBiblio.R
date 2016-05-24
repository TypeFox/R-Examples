MacroBiblio <-function(base, num.text="Abstract", num.agg="Year", idiom ="en",lminword=3, Fmin=10, Dmin=5,Fmax=NULL,
             equivalence=NULL,stop.word.user=NULL, lmd=3, lmk=3, ncp=10, row.sup=NULL, col.sup=NULL, 
              graph=TRUE, axes=c(1,2), proba=0.01){


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

#Contextual variables
 VarSel<-DocVarTable(base, VarSel=c("Year","Journal" ))
 DocVar<-as.data.frame(VarSel[,c("Year","Year")])
 DocVar[,1]<-as.numeric(DocVar[,1])
 DocVar[,2]<-as.factor(DocVar[,2])
   summary(DocVar)

# Corpus and DocWordTable	
 Corpus<-DocWordTable(base, num.text, idiom , lminword,  
                      lower = TRUE, Remov.Number = TRUE)             
        DocTerm<-Corpus$DocTerm
        
# direct  CA for  Metakeys and Metadocs 
 res.TxCA.Dir<-TxCA(DocVar,DocTerm, Fmin, Dmin, Fmax, idiom=idiom, equivalence, num.agg = NULL,stop.word.user,
                  lmd,lmk,ncp, row.sup=NULL,col.sup ,stop.word.tm=TRUE,graph=FALSE, axes,quanti.sup =1,quali.sup=2)
 			  
	 res.caD<-res.TxCA.Dir$res.ca  
       DocTermR<-res.TxCA.Dir$DocTermR
       res.meta<-res.TxCA.Dir$res.meta
       Glossary<-res.TxCA.Dir$Glossary
       DimWords<- res.meta$DimWord[which(res.meta$DimWord[,1]>1),]
      
# Aggregate CA
 if(num.agg=="Year"){
 DocTerm<-Corpus$DocTerm
 res.TxCA.Agreg<-TxCA(DocVar,DocTerm, Fmin, Dmin,Fmax, idiom=idiom, equivalence, num.agg=1,stop.word.user,
                   lmd,lmk,ncp, row.sup,col.sup ,stop.word.tm=TRUE,graph=FALSE,axes,quanti.sup=NULL,quali.sup=NULL)
	
 }else{
    DocTerm<-Corpus$DocTerm
   res.TxCA.Agreg<-TxCA(DocVar=base,DocTerm, Fmin, Dmin,Fmax, idiom=idiom, equivalence, num.agg=num.agg,stop.word.user,
                   lmd,lmk,ncp, row.sup,col.sup ,stop.word.tm=TRUE,graph=FALSE,axes,quanti.sup=NULL,quali.sup=NULL)
 }

	 res.caAg<-res.TxCA.Agreg$res.ca
       Table<-res.TxCA.Agreg$Table

#characteristic words	
     CharWord.Compl<-descfreq(Table, by.quali = NULL, proba)
	CharWord=vector(mode="list",length=length( CharWord.Compl))
	names(CharWord)<-names( CharWord.Compl)
      for(i in 1:length( CharWord.Compl)){
            bWord <- as.data.frame( CharWord.Compl[[i]])
            over <- subset(bWord, bWord$v.test > 0)
            Over <- row.names(over)
            infra <- subset(bWord, bWord$v.test < 0)
            Infra <- row.names(infra)
            OverInfra <- list(Over, Infra)
            names(OverInfra) <- c("Over_represented_word", "Infra_represented_word")
            CharWord[[i]] <- OverInfra
        }
        CharWords<-CharWord
        res.CharWord<-list( CharWord.Compl= CharWord.Compl, CharWords=CharWords)
   
# Constrained hierarchical clustering 
 res.chcpc<-TxCHCPC(res.caAg,cluster.CA = "rows",nb.clust=-1, graph=FALSE)

#cronological evolution AFMTC
  DocTermR<-res.TxCA.Dir$DocTermR
  DocTermRM<-DocTermR[apply(DocTermR,1,sum)>0,]
  DocVarR<-as.data.frame(DocVar[which(rownames(DocVar)%in%rownames(DocTermRM)),])
  MDocWord<-cbind.data.frame(DocTermRM,DocVarR)

  res.mfact<-MFA(MDocWord, group=c(dim(DocTermRM)[2],1,1), type =c("f","s","n"),
                name.group=c("FREQ","YEAR","YEAR_CAT"), num.group.sup = 3,graph=FALSE)    

# most contributive word in the DIM1	
   WordDim1<-as.data.frame((res.mfact$freq$coord[,1]))
   WordDim1$Words<-as.data.frame(row.names(WordDim1))
   WordDim<-cbind(WordDim1[,2],round(WordDim1[,1],3))
   colnames(WordDim)<-c("Word", "Coord.Dim1")
   WordDim1Comp<-with( WordDim, WordDim[order(Coord.Dim1,Word),])
   sel<-which(res.mfact$freq$contrib[,1]>lmk*mean(res.mfact$freq$contrib[,1]))
   WordMoreContrib<- WordDim1Comp[which( WordDim1Comp[,1]%in%names(sel)),]  
 res.WordDim1<-list(WordMoreContrib= WordMoreContrib, WordDim1Comp=WordDim1Comp)
   
#pioneers articles 
  SelYear<-as.numeric(rownames(subset(res.mfact$quali.var.sup$coord,
         res.mfact$quali.var.sup$coord[,1]>0)))
  SelYear<-SelYear[1:(length(SelYear)-2)]

  MDocWordA<-MDocWord[which(MDocWord[,ncol(MDocWord)-1]%in%SelYear),]
  MDocWordA<-MDocWordA[,apply(MDocWordA[,1:(ncol(MDocWordA)-2)],2,sum)>0]
	mfact.res<-MFA(MDocWordA, group=c((ncol(MDocWordA)-2),1,1), type =c("f","s","n"),
                name.group=c("FREQ","YEAR","YEAR_CAT"), num.group.sup = 3,graph=FALSE)  

###  distance between the partial points
	G1<-mfact.res$ind$coord.partiel[seq(1,(ncol(MDocWordA)-1),2),]
      G1<-G1[which(apply(G1,1,sum)!=0),]
	G2<-mfact.res$ind$coord.partiel[seq(2,(ncol(MDocWordA)),2),]
	G2<-G2[which(apply(G2,1,sum)!=0),]
      G<-G1-G2

 if(graph){

#graph of AFMTC
      lim=0.8
         plot(1,1,pch=16,col="white",xlim=c(-0,4),ylim=c(-4,5),xlab=paste("DIM",axes[1],"(",round(mfact.res$eig[axes[1],2],2),"%)",sep=""),
		ylab=paste("DIM",axes[2],"(",round(mfact.res$eig[axes[2],2],2),"%)",sep=""),main="Superimposed Representation (partial axes)",cex=0.75)
      points(G2[,1],G2[,2],pch=16,col="white",cex=0.75)
	points(G1[,1],G1[,2],pch=16,col="white",cex=0.75)
	points(mfact.res$ind$coord[,1],mfact.res$ind$coord[,2],pch=16,col="white",cex=0.75)
	segments(G1[,1],G1[,2],mfact.res$ind$coord[,1],mfact.res$ind$coord[,2],lty=2,col="white",cex=0.75)
	segments(mfact.res$ind$coord[,1],mfact.res$ind$coord[,2],G2[,1],G2[,2],lty=1,col="white")
	abline(h=0,v=0,lty=2)
	points(G1[which(G[,1]>lim),1],G1[which(G[,1]>lim),2],pch=16,col="red",cex=0.75)
	points(G2[which(G[,1]>lim),1],G2[which(G[,1]>lim),2],pch=16,col="green",cex=0.75)
	points(mfact.res$ind$coord[which(G[,1]>lim),1],mfact.res$ind$coord[which(G[,1]>lim),2],pch=16,col="black",cex=0.75)

	text(mfact.res$ind$coord[which(G[,1]>lim),1],mfact.res$ind$coord[which(G[,1]>lim),2],rownames(mfact.res$ind$coord[which(G[,1]>lim),]),col="black",cex=0.75,font=2,pos=c(3,4,4,3,1,3,1))
		
	segments(G1[which(G[,1]>lim),1],G1[which(G[,1]>lim),2],mfact.res$ind$coord[which(G[,1]>lim),1],mfact.res$ind$coord[which(G[,1]>lim),2],lty=2,col="red",cex=0.75)
	segments(mfact.res$ind$coord[which(G[,1]>lim),1],mfact.res$ind$coord[which(G[,1]>lim),2],G2[which(G[,1]>lim),1],G2[which(G[,1]>lim),2],lty=1,col="green")
	legend("topleft",c("Vocabulary","Year"),lty=c(2,1),text.col=c("red","green"),col=c("red","green"))

   plot.MFA(res.mfact, choix = "ind", invisible = "ind",  habillage = "group", 
                axes = axes, new.plot = TRUE)
    points(res.mfact$quali.var.sup$coord[,1],res.mfact$quali.var.sup$coord[,2]
                ,type="l",cex=0.75,pch=19,col="blue",lwd=2,lty=2)
     dev.new()
     sel1<-which( (res.mfact$freq$contrib[,1]> lmk*mean(res.mfact$freq$contrib[,1])) | 
                   (res.mfact$freq$contrib[,2]>lmk*mean(res.mfact$freq$contrib[,2])) )
      par(cex=0.8)
	plot.MFA(res.mfact,choix="freq",invisible="row",select=sel1,axes,
			unselect=1,palette=palette(c("black","black","black")),col.hab=c("green","blue"),title="Words" )

    plot.MFA(res.mfact, choix = "var", habillage = "group", 
                axes = axes, new.plot = TRUE, shadowtext = TRUE) 
    plot.MFA(res.mfact, choix = "group", axes = axes, new.plot = TRUE) 

     ##Graph of TxCHCPC
     res.chcpc<-TxCHCPC(res.caAg,cluster.CA = "rows",nb.clust=-1)	
     #Graph aggregate AC
      res.metaA<- META.CA(res.caAg, naxes = ncp, axe.x = axes[1],  axe.y = axes[2], lmd = -Inf, lmk, main = "CA documents/words")
         dev.new()
         plot(res.caAg, invisible = c("col", "col.sup","quali.sup"), axes = axes,  title = "CA documents")
	   points(res.caAg$row$coord[,1],res.caAg$row$coord[,2],type="l",cex=0.75,pch=19,col="black",lwd=2,lty=2)
        dev.new()
	  #Graph direct AC 
         res.meta <- META.CA(res.caD, naxes = ncp, axe.x = axes[1],axe.y = axes[2], lmd, lmk)  
         dev.new()        
         plot(res.caD, invisible = c("col", "col.sup","quali.sup"), axes = axes,  title = "CA documents")
         dev.new()
         res.meta1<- META.CA(res.caD, naxes = ncp, axe.x = axes[1],  axe.y = axes[2], lmd = Inf, lmk, main = "CA words")
   
   ###pioneers articles
        pioneers<-base[row.names(base)%in%names(which(mfact.res$ind$coord[which(G[,1]>lim),1]>0)),]               
        pioneers<- pioneers[,colnames(pioneers)%in%c("Title", "Year","Journal")]      
} 
  else pioneers<-NULL  
res<-list(Corpus=Corpus,Glossary=Glossary,DocTermR=DocTermR,Tagreg=Table, Metakeys.Metadocs=res.meta, res.CA=res.TxCA.Dir, res.CA.Agreg=res.TxCA.Agreg,CharWord=res.CharWord,
        res.CHCPC=res.chcpc, res.MFACT=res.mfact,OrdWord=res.WordDim1,pioneers=pioneers)
 class(res)<-c("MacroBiblio","list")
return(res)
}
