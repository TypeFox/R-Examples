delList <- function(vector,maxdist){
  bucket<-list()
  if(maxdist==1){
    for(word in vector){
      bucket[[word]]<-deletion1(word)
    }
  }
  else if(maxdist==2){
    for(word in vector){
      bucket[[word]]<-deletion2(word)
    }
  }
  else if(maxdist==3){
    for(word in vector){
      bucket[[word]]<-deletion3(word)
    }
  }
  else{stop("argument \"maxdist\" out of range (1:3)")}
  bucket
}
