
close.bracket<-function(pattern, x){

  possible.end<-nchar(x)+1
  start<-gregexpr(pattern, x)[[1]]
  n<-length(start)
  openB<-gregexpr("\\(", x)[[1]]
  closeB<-gregexpr("\\)", x)[[1]]
  end<-1:n

  for(i in 1:n){
    nopen<-1
    pos<-openB[match(TRUE, openB>start[i])]+1
    while(nopen!=0){
      if(pos%in%openB){nopen<-nopen+1}
      if(pos%in%closeB){nopen<-nopen-1}
      pos<-pos+1
      if(pos>possible.end){stop("path formuala invalid")}
    }
    end[i]<-pos-1
  }
  cbind(start,end)
}



