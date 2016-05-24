HET<-
function(data){
  genos<-as.matrix(data)
  n.genos<-ncol(genos)
  m<-nrow(genos)
  #Make seperate matrices, describing number of A and B alleles
  mat.a<-matrix(0,nrow=m,ncol=n.genos)  
  mat.b<-matrix(0,nrow=m,ncol=n.genos)
  mat.a[genos==1]<-1
  mat.a[genos==0]<-0.5
  mat.b[genos==-1]<-1
  mat.b[genos==0]<-0.5
  #Tally the number of A and B alleles at each locus
  a<-rowSums(mat.a)
  b<-rowSums(mat.b)
  #Calculate mean HET
  return(mean(1-(a^2+b^2)/(a+b)^2,na.rm=TRUE))
}

