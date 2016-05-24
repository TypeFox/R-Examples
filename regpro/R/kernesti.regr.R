kernesti.regr<-function(arg,x,y,h=1,kernel="gauss",g=NULL,gernel="gauss",
vect=FALSE)
{
w<-kernesti.weights(arg,x,h=h,kernel=kernel,vect=vect,g=g,gernel=gernel)
y<-matrix(y,length(y),1)

if (!vect){
  notna<-(!is.na(y))
  w<-matrix(w,length(w),1)    
  w<-notna*w
  w<-w/sum(w)
  mu<-w*y
  est<-sum(mu,na.rm=TRUE)
}
else{
  est<-t(w)%*%y  #t(y)%*%w
}

return(est)
}



