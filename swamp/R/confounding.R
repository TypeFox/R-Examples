confounding <-
function(o,method="chisq",workspace=2e7,smallest=-20,diagonal.zero=F,label=colnames(o),note=T,notecol="black",notecex=1,breaks=50,col=c(heat.colors(48),"white"),key=T,cexRow=1,cexCol=1,margins=c(7,7),colsep=NULL,rowsep=NULL,sepcolor="black",sepwidth=c(0.05,0.05)){

      if(smallest>=0){stop("smallest has to be less than 0")}
      classes<-unlist(lapply(unclass(o),class))
      if(any(classes=="character")){stop("o contains characters")}
      nrlevels<-unlist(lapply(unclass(o),function(x)length(levels(x))))
      if(any(nrlevels==1)){stop("o contains factors with only one level")}
    matp<-matrix(nrow=ncol(o),ncol=ncol(o),dimnames=list(colnames(o),colnames(o)))
    sums<-matp
    test<-matp
    for(i in 1:ncol(o)){
      for(j in 1:ncol(o)){
                       xi<-table(o[!is.na(o[,j]),i])
                       li<-length(xi[xi>0])
                       xj<-table(o[!is.na(o[,i]),j])
                       lj<-length(xj[xj>0])
                  if(any(li<2,lj<2)){warning(paste("Annotations o$",colnames(o)[i]," and o$",colnames(o)[j]," in combination dont have 2 different values each. Output is set to NA.",sep=""))
                       matp[i,j]<-NA
                       sums[i,j]<-NA
                       test[i,j]<-"none"                       
                       }
                  else{                      
            if(classes[i]=="factor"&classes[j]=="factor"){                                                                        
                if(method=="fisher"){
                matp[i,j]<-fisher.test(table(o[,i],o[,j]),workspace=workspace)$p.value
                sums[i,j]<-sum(table(o[,i],o[,j]))
                test[i,j]<-"fisher.test" }
                    if(method=="chisq"){                    
                    tab<-table(o[,i],o[,j])
                    which(rowSums(tab)==0)
                    tab<-tab[which(rowSums(tab)!=0),(which(colSums(tab)!=0))] #get rid of the factors with no entry (fisher test does this by default)                    
                    matp[i,j]<-chisq.test(tab)$p.value
                    sums[i,j]<-sum(tab)
                    test[i,j]<-"chisq.test"
                    }}
            if(classes[i]%in%c("factor","numeric","integer")&classes[j]%in%c("numeric","integer")){
                s<-summary(lm(o[,j]~o[,i]))
                matp[i,j]<-pf(s$fstatistic[1],s$fstatistic[2],s$fstatistic[3],lower.tail=FALSE)
                sums[i,j]<-nrow(o)-length(s$na.action)
                test[i,j]<-"lm"}
            if(classes[j]%in%c("factor","numeric","integer")&classes[i]%in%c("numeric","integer")){
                s<-summary(lm(o[,i]~o[,j]))
                matp[i,j]<-pf(s$fstatistic[1],s$fstatistic[2],s$fstatistic[3],lower.tail=FALSE)
                sums[i,j]<-nrow(o)-length(s$na.action)
                test[i,j]<-"lm"}
                  }
                }}   
   require(gplots)
   if(diagonal.zero==T){matp[-(which(upper.tri(matp)|lower.tri(matp)))]<-0}
   linp10<-log10(matp)
   linp10<-replace(linp10,linp10<=smallest,smallest)
   heatmap.2(linp10,Colv=F,Rowv=F,dendrogram="none",trace="none",symbreaks=F,symkey=F,breaks=breaks,key=key,col=col,cexRow=cexRow,cexCol=cexCol,colsep=colsep,rowsep=rowsep,sepcolor=sepcolor,sepwidth=sepwidth,
            main="",labCol=label,labRow=label,xlab="",margins=margins,cellnote=if(note==T){signif(matp,1)}else{matrix(ncol=ncol(matp),nrow=nrow(matp))},notecol=notecol,notecex=notecex)
            
      return(list(p.values=matp,n=sums,test.function=test,classes=classes))
      }
