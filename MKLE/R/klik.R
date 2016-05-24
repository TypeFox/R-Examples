"klik" <- function(delta,data,kde,grid,min){
  dshift<-data+delta
  vec<-(dshift-min)/grid
  index<-round(vec)+1
  sum(log(kde$y[index]+(kde$y[index+1]-kde$y[index])/grid*(dshift-kde$x[index])))
}
