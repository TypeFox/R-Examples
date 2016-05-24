add_es <- function(eg,eg2,current_rank,orgn,ff=0,method=c("esm","isvd")){
  if(method=="esm"){
    out = add_eig(eg, eg2)}
  else{
    if (is.null(eg$m)) {
      #without orgn data is assumed as zero-mean
      m = dim(eg$u)[1]
      print(m)
    } 
    else
    {
      m = eg$m
    }
    if (missing("current_rank")) {
      #full rank
      current_rank =  2
    }
    B = eg2
    out = add_svd(eg,B,m,current_rank,orgn,ff = ff)
  }
  return(out)
}