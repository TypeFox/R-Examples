show.map<-function(){
  #data(ALINE.map,envir=environment())
  show.map<-cbind(apply(data.frame(ALINE.map$U.Val),MARGIN=1,FUN=intToUtf8),ALINE.map,deparse.level=0)
  colnames(show.map)[1]<-"IPA"
  return(show.map)
}