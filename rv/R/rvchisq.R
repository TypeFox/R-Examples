

rvchisq <- function (n=1, df, ncp = 0) 
{
  if (missing(ncp)) {
    rvvapply(stats:::rchisq, n.=n, df=df)
  } else {
    rvvapply(stats:::rchisq, n.=n, df=df, ncp=ncp)
  }
}


