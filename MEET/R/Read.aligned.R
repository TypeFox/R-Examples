Read.aligned <-function(file){
  DNA<-read.table(file, colClasses="character")
  DNA<-as.vector(DNA[[1]])
  sequence<-sapply(DNA,function(x){strsplit(x,"")}, USE.NAMES=FALSE)
  m<-matrix(,length(sequence)*0.5,length(sequence[[2]]))
	for(i in c(1:(length(sequence)*0.5))){
    m[i,]<-sequence[[2*i]]
  }

  m
}

