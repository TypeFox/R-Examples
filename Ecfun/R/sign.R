sign <- function(x, zero=0){
  s <- (-1*(x<0)+zero*(x==0)+(x>0))
  s
}