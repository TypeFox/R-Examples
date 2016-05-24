
mesh <- function(x,dim){  
  nrows <- length(x)^(dim)
  out <- array(dim=nrows*dim)
  idx = 1;
  for(d in 1:dim){
    r = length(x)^(dim-d) 
    while(idx < d*nrows){
      for(v in x){
        out[idx:(idx+r-1)] = v
        idx = idx + r;
      }
    }
  }
  dim(out) <- c(nrows,dim)
  return(out)
}