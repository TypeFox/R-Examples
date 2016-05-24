Quant <-
function(p,int,tolerance=.Machine$double.eps,...){
  F <- function(x)1-exp(-CumInt(x,int=int,...))
  fn <- function(p){
    lo <- up <- k <- 0;
    while(F(up)<p){up <- up + 2^k; k <- k+1}
    while(abs(lo-up)>tolerance){
      mid <- (lo+up)/2;
      if(F(mid)>=p){up <- mid}else{lo <- mid}
    }
    return((lo+up)/2)
  }
  sapply(p,fn)
}

