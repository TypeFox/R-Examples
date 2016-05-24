# x integer
# r vector of 0/1/2
int2tern<-function(x){
  r<-vector(mode="integer",length=0)
  if (x<=2) r<-x
  else{
  while(x>2){
    r<-c(x%%3,r)
    x<-floor(x/3)
  }
  r<-c(1,r)
  }
  r
}
