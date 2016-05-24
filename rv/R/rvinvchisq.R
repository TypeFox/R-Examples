

rvinvchisq <- function (n=1, df, scale=1) {
  return(scale / (rvchisq(n=n, df=df) / df))
}


