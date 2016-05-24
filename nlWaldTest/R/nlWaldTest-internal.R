.getint <-
function (cf, parm,ses, level = 0.95, df = NULL,...) 
{
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  if (is.null(df))
    fac <- qnorm(a) 
  else
    fac <- qt(a, df)
  fac <- c(0,fac)
  fff = format(100*a, 3)
  pct=paste(fff,"%")
  pct=c("value",pct)
  ci <- array(NA, dim = c(length(parm), 3L), dimnames = list(parm,pct))
  ci[] <- cf + ses %o% fac
  ci
}
.smartsu <-
function(pat,repl,x)
{
  mx=paste("?",x,"?",sep="")
  bign=nchar(mx)
  n=nchar(pat)
  if (grepl(pat,x))
  {
    poss= gregexpr(pat,mx)[[1]]
    poss=as.numeric(poss)
    npos=length(poss)
    xvect=c(substr(mx,1,poss[1]-1))
    for (i in 1L:npos)
    {
      txti1=substr(mx,poss[i]-1,poss[i]+n-1)
      txti2=substr(mx,poss[i],poss[i]+n)
      tx=substr(mx,poss[i]-1,poss[i]-1)
      boo=grepl(tx,"0123456789",fixed=T)
      if ((txti1!=make.names(txti1)&&!boo)&&txti2!=make.names(txti2))
        xvect=paste(xvect,repl,sep="")
      else
        xvect=paste(xvect,pat,sep="") 
      if (i==npos)
        xvect=paste(xvect,substr(mx,poss[i]+n,bign),sep="")
      else
        if (poss[i+1]>(poss[i]+n))
          xvect=paste(xvect,substr(mx,poss[i]+n,poss[i+1]-1),sep="")  
    }
    bign=nchar(xvect)
    return(substr(xvect,2,bign-1))
  }
  else
    return(x)
}
.smartsub <- Vectorize(.smartsu,vectorize.args = "x",USE.NAMES =F)

