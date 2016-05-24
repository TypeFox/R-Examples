sym.histogram.pca <-
function(sym.data,method=c('histogram','classic')) {
  method<-match.arg(method)  
  if(method=='histogram') {
    dam<-downarrow.matrix(sym.data)
    ram<-rightarrow.matrix(sym.data)
    cpc<-PCA(ram,graph=FALSE)
    res<-sym.interval.pca(dam,'centers')
    q<-min(res$Sym.Components$M,dim(cpc$ind$coord)[2])
  }
  k<-max(sym.data$sym.var.length)
  if(k==1) {
    return(res)    
  }  
  else {
    pos<-1
    for(i in 1:sym.data$N) {
      for(s in 1:k) {
        colm<-2
        for(j in 1:q) {
          res$Sym.Components$meta[pos,colm]<-
            res$Sym.Components$meta[pos,colm]+cpc$ind$coord[i,j]
          res$Sym.Components$meta[pos,colm+1]<-
            res$Sym.Components$meta[pos,colm+1]+cpc$ind$coord[i,j]
          colm<-colm+3
        }
        pos<-pos+1
      }
    }
    dam$meta<-res$Sym.Components$meta
    return(dam)      
  }
}
