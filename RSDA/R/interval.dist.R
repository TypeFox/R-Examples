interval.dist <-
function(sym.data,distance=c('hausdorff','centers','interscal'),p=2) {
  distance<-match.arg(distance)
  if((distance=='hausdorff')||(distance=='centers')) {
    idn<-all(sym.data$sym.var.types=='$I')
    if(idn==FALSE)
      stop("The two variables have to be interval type")         
    nn<-sym.data$N
    mm<-sym.data$M
    dist.matrix<-matrix(0,nn,nn)
    for(i in 1:nn) {
      for(j in 1:i) {
        dist.matrix[i,j]<-interval.dist.tobj(sym.obj(sym.data,i),sym.obj(sym.data,j),
                                             distance)
        dist.matrix[j,i]<-dist.matrix[i,j]
      }      
    }
    return(as.dist(dist.matrix))
  }  
  if((distance=='interscal')) {
    idn<-all(sym.data$sym.var.types=='$I')
    if(idn==FALSE)
      stop("The two variables have to be interval type")         
    nn<-sym.data$N
    mm<-sym.data$M
    dist.matrix<-matrix(0,nn,2*nn)
    min.matrix<-matrix(0,nn,nn)
    max.matrix<-matrix(0,nn,nn)
    pos<-1
    for(j in 1:nn) {      
      for(i in 1:nn) {
          dist.matrix[i,pos]<-interval.dist.tobj(sym.obj(sym.data,i),sym.obj(sym.data,j),distance)[1]
          dist.matrix[i,pos+1]<-interval.dist.tobj(sym.obj(sym.data,i),sym.obj(sym.data,j),distance)[2]
          min.matrix[i,j]<-dist.matrix[i,pos]
          max.matrix[i,j]<-dist.matrix[i,pos+1]
      }     
      pos<-pos+2
    }
    return(list(interval.dist=dist.matrix,min.matrix=min.matrix,max.matrix=max.matrix))
  }  
}
