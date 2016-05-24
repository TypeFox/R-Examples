crearFullDT <- function(input,maxdist){
  V1 <- c()
  if(maxdist < 1 || maxdist > 3){stop("argument \"maxdist\" out of range (1:3)")}
  else{
  DT <- crearDT()
  if(!file.exists(input)[1]){
    for(item in input){
      v1<-matricear(item,maxdist)
      DT<-rbindlist(list(v1,DT))
    }
  }
  else if (is.null(input)) {return()}
  else{
    content <- readLines(input, warn = FALSE)
    for(item in content){
      v1<-matricear(item,maxdist)
      DT<-rbindlist(list(v1,DT))
    }
  }
  setkey(DT,V1)
  if(maxdist==1){
    class(DT)<-append(class(DT), "dist1_DT")
  }
  else if(maxdist==2){
    class(DT)<-append(class(DT), "dist2_DT")
  }
  else if(maxdist==3){
    class(DT)<-append(class(DT), "dist3_DT")
  }
  }
  return(DT)
}
