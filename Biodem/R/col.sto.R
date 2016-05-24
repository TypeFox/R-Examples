"col.sto" <- function(x){
  y<-apply(x,2,sum)
  x1<-t(t(x)/y)
  x1
}
