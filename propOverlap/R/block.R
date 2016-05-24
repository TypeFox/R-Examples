block <- function( ... ){
  argv = list( ... )
  i = 0
  for(a in argv){
    m = as.matrix(a)
    if(i == 0)
      rmat = m
    else
    {
      nr = dim(m)[1]
      nc = dim(m)[2]
      aa = cbind(matrix(0,nr,dim(rmat)[2]),m)
      rmat = cbind(rmat,matrix(0,dim(rmat)[1],nc))
      rmat = rbind(rmat,aa)
    }
    i = i+1
  }
  rmat    
}