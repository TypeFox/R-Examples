## Julia Bischof
## 15-10-2015

#library(ineq)

clones.giniIndex<-function(clone.size=NULL, PDF=NULL){
  
  if(length(clone.size)==0){
    stop("--> Clone size vector is required")
  }
  
  gini<-ineq(sort(clone.size,decreasing = F), type="Gini")
  
  if(length(PDF)>0){
    pdf(paste(PDF,"_Lorenz-curve.pdf",sep=""),width=7,height=7, pointsize = 16)
    par(mar=c(5,5,4,3))
    plot(Lc(as.numeric(sort(clone.size,decreasing = F)),
            as.numeric(sort(clone.size)/sum(clone.size,na.rm=T))),col="darkred",
         main=c("Lorenz curve",paste("(Gini Index = ",round(gini,2),")",sep="")))
    dev.off()
  }
  return(gini)
}
