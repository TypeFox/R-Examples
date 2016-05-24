
thintrack<-function(spdf,mindist=100){
  st<-FALSE
  cds<-coordinates(spdf)
  cds<-cds[!duplicated(cds),]
  nbb<-1
  cds2<-NULL
   while(!st){
    cds2<-rbind(cds2,cds[nbb,])
    mycirc<-polycirc(mindist,unlist(cds[nbb,]))
    idx<-inout(cds,mycirc)    
    cds<-cds[!idx,,drop=FALSE]
    if (nrow(cds)>1) {
        nnb<-knearneigh(rbind(cds2[nrow(cds2),],cds),k=1)
        nbb<-nnb$nn[1,1]-1
    } else st<-TRUE
  }
  cds2<-data.frame(cds2)
  names(cds2)<-c("x","y")
  coordinates(cds2)<-~x+y
  proj4string(cds2)<-proj4string(spdf)
  cds2
}