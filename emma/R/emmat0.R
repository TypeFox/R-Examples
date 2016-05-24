emmat0  <- function(in.name,nlev,lower,upper,out.name,nd,fn1=NULL,fn2=NULL,fn3=NULL,fn4=NULL)
{

  ## load packages
  library(earth)
  library(clusterSim)

  xspace <- generate(in.name,nlev,lower,upper)
  yspace <- data.frame(matrix(NA,nrow=nrow(xspace),ncol=length(out.name)))
  colnames(yspace) <- out.name

#  xpop <- initialise(xspace,nd)
  sam <- sample(seq(1:nrow(xspace)),nd)
  xpop <- xspace[sam,]

  tested <- sam

  print("PERFORM THE FOLLOWING EXPERIMENTS ( t = 0 )")
  print(xpop)

  if (!is.null(fn1)) ypop <- fn1(xpop)
  if (!is.null(fn2)) ypop <- cbind(ypop,fn2(xpop))
  if (!is.null(fn3)) ypop <- cbind(ypop,fn3(xpop))
  if (!is.null(fn4)) ypop <- cbind(ypop,fn4(xpop))
  if (is.null(fn1))
  { 
    ypop <- NULL
    print("ADD THE MEASURED RESPONSE VALUES TO ypop")
  }

  structure(list(xpop=xpop,ypop=ypop,xspace=xspace,yspace=yspace,nd=nd,tested=tested,time=0),class="emmat0")
}
