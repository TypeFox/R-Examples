NormalizeRUVRand <- function(Y, ctl, k=NULL,lambda=NULL,plotk=TRUE){
  output<-RUVRand(Y=Y, ctl=ctl,lambda=lambda, k=k)
  return(structure(output, class="normdata"))  
}
  