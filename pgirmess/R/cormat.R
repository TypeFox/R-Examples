"cormat" <-
function(donnees,method="spearman",sep=FALSE){
  
    noms<-names(donnees)
    cells<-matrix(data=rep(0,length(noms)^2),nrow=length(noms),ncol=length(noms))
    cellsp<-cells
    for (i in 1:length(noms)) {cells[i,i]<-1; cellsp[i,i]<-0}
    debut<-1
        for (i in debut:length(noms)){
            for (j in (i+1):length(noms)){
                    if(i!=length(noms)){
                    res<-cor.test(donnees[,i],donnees[,j],method=method)
                    if (i==1) met<-res$method
                    cells[i,j]<-res$estimate;cells[j,i]<-res$estimate
                    cellsp[i,j]<-res$p.value;cellsp[j,i]<-res$p.value
                    }
            }
        debut<-debut+1
        }
        cells<-as.data.frame(cells)
        names(cells)<-names(donnees)
        row.names(cells)<-noms
        cellsp<-as.data.frame(cellsp)
        names(cellsp)<-names(donnees)
        row.names(cellsp)<-noms
        cells<-round(cells,2)
        cellsp<-round(cellsp,3)
        if (sep)return(list(method=met,coef.estimates=cells,p.values=cellsp))
        else {
            mat<-cells;mat[lower.tri(mat)]<-cellsp[lower.tri(cellsp)]
            return(list(method=met,prob.cor=mat))
        }
}
