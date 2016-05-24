Dens <-
function(x,int,...){
  sapply(x,int)*exp(-CumInt(x,int=int,...))
}

