ora<-function(geneList,GPSrepo){
    idmap<-get(data(idmap,envir=as.environment(parent.frame())))
    fr<-GPSrepo$origRepo[[3]]
       sp1<-GPSrepo$pathwaydescriptions
          if(length(intersect(geneList,GPSrepo$origRepo[[2]]))==0){
              t1<-which.max(sapply(idmap,function(x)length(intersect(x,GPSrepo$origRepo[[2]]))))
            t2<-which.max(sapply(idmap,function(x)length(intersect(x,geneList))))
            geneList<-idmap[which(idmap[,t2]%in%geneList),t1]
            print(paste("Mapped identifiers from" ,colnames(idmap)[t2]," to ",colnames(idmap)[t1],"..."))
             }
    g1m<-match(geneList,GPSrepo$origRepo[[2]])   
    frO<-fr[fr[,2]%in%g1m,]
    npwys<-table(as.character(GPSrepo$origRepo[[1]][(frO[,1])]))
    nn<-cbind(sp1[match(names(npwys),sp1[,1]),],as.numeric(as.vector(npwys)))
    PPwys<-table(as.character(GPSrepo$origRepo[[1]][(fr[,1])]))
    ps<-phyper(npwys-1,PPwys[match(names(npwys),names(PPwys))],
               length(GPSrepo$origRepo[[2]])-PPwys[match(names(npwys),names(PPwys))],length(intersect(geneList,GPSrepo$origRepo[[2]])),lower.tail=F)
    ps<-signif(ps,digits = 4)
    res<-data.frame(nn,PPwys[match(names(npwys),names(PPwys))],as.numeric(ps),as.numeric(p.adjust(ps,method='bonfer')))
    colnames(res)<-c("pathwyid","description","success","pathwaySize","pvalues","Bonfer")
    res<- res[with(res, order(pvalues)),]
   return(res[which(as.numeric(res[,ncol(res)])<0.05),])
}
