m2rasch.func <- function(nitem, ncat, nrec, myX, alphas, betas, nq, iprint){

  ## prep return variables
  g<-nitem*(ncat-1)
  for (i in 1:(nitem-1))
  {
    for (j in 0:(i-1))
    {
      for (k in 1:(ncat-1))
      {
        for (l in 1:(ncat-1))
        {
          g<-g+1
        }
      }
    }
  }
  
  samplemout=rep(0,g)
  m2statout=0
  dfout=0
  
  ## TO DO: add package arument here?
  out <- .C("Rm2rasch",
    as.integer(nitem), as.integer(ncat), as.integer(nrec), as.double(myX),
    as.double(alphas), as.double(betas), samplemout=as.double(samplemout),
    m2statout=as.double(m2statout), dfout=as.double(dfout), as.integer(nq))  
  list(samplem=out$samplemout, m2stat=out$m2statout, df=out$dfout)  
}
