TxMFACT <-
function(MDocWord, group, type = rep("s",length(group)), col.sup=NULL, ind.sup = NULL, 
           ncp = 5, name.group = NULL, num.group.sup = NULL, 
          graph = TRUE, weight.col.mfa = NULL, row.w = NULL,axes = c(1,2),tab.comp=NULL)
 {
    if (!is.null(col.sup)){
         if(is.character(col.sup))
           col.sup<-which(colnames(MDocWord)%in%col.sup)
          if(is.numeric(col.sup))
	     col.sup<-col.sup  
   MSup<-MDocWord[,col.sup]
   MDocWordR<-MDocWord[,-col.sup]
   MDocWord<-cbind.data.frame(MDocWordR,MSup)
       } 


    if (!is.null(ind.sup)){
         if(is.character(ind.sup))
           ind.sup<-which(rownames(MDocWord)%in%ind.sup)
          if(is.numeric(ind.sup))
	     ind.sup<-ind.sup
       } 

  mfact<- MFA (MDocWord, group, type, ind.sup, ncp, name.group, num.group.sup, 
              graph, weight.col.mfa, row.w,axes,tab.comp)
	if(graph){
      dev.new()
	sel1<-which( (mfact$freq$contrib[,1]> 2*mean(mfact$freq$contrib[,1])) | 
                   (mfact$freq$contrib[,2]>2*mean(mfact$freq$contrib[,2])) )

      par(cex=0.7)
	plot.MFA(mfact,choix="freq",invisible="row",select=sel1,axes,
			unselect=1,palette=palette(c("black","black","black")),col.hab=c("green","blue"),title="Words" )
	
	if (!is.null(col.sup)){
            selSup<-colnames(MSup)
            sel2<-which(rownames(mfact$freq.sup$cos2)%in%selSup)	
	     sel<-c(sel1,sel2)
       dev.new()
	par(cex=0.7)
	plot.MFA(mfact,choix="freq",invisible="row",select=sel,axes,
		unselect=1,palette=palette(c("black","black","black")),col.hab=c("green","blue"),title="Words")
	}
	dev.new()
      plot.MFA(mfact,choix="ind",axes, title="Documents")
}
  res<-mfact
 class(res)<-c("TxMFACT","MFA","list")
  return(res)
}
