# starting program ncp

ncp <- function(df1, df2, alpha, beta)
{
  f.min <- 0
  f.max <- 1000000000
  bRrepeat <- TRUE
  while(bRrepeat) {
    f.mid <- 0.5 * (f.min + f.max)
    rpower <- 1-pf(qf(1-alpha,df1,df2),df1,df2,f.mid)
    if(rpower < 1 - beta) {
      f.min <- f.mid
    }
    else {
      f.max <- f.mid
    }
    if(abs(1 - beta - rpower) < 1e-008) {
      bRrepeat <- FALSE
    }
  }
  return (f.mid)
} 


