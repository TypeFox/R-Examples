#--------------------------------------------------------------------------------------------
ars<-function(n=1,f,fprima,x=c(-4,1,4),ns=100,m=3,emax=64,lb=FALSE,ub=FALSE,xlb=0,xub=0,...)
{
  mysample<-rep(0,n)
  iwv<-rep(0,ns+7)
  rwv<-rep(0,6*(ns+1)+9)
  hx<-f(x,...)
  hpx<-fprima(x,...)
  initial<-.C("initial_",as.integer(ns),as.integer(m),
                  as.double(emax),as.double(x),
                  as.double(hx),as.double(hpx),
                  as.integer(lb),as.double(xlb),
                  as.integer(ub),as.double(xub),
                  ifault=as.integer(0),
                  iwv=as.integer(iwv),
                  rwv=as.double(rwv))
  if(initial$ifault==0)
  {
    h<-function(x) f(x,...)
    hprima<-function(x) fprima(x,...)
    for(i in 1:n)
    {
       sample<-.C("sample_",as.integer(initial$iwv),as.double(initial$rwv),h,hprima,new.env(),
                  beta=as.double(0),ifault=as.integer(0))
       if(sample$ifault==0)
       {
          if(i<ns)
          {
            x<-c(x,sample$beta)                          #Add the sampled point to the support and update h and hprima
            h<-function(x) f(x,...)
            hprima<-function(x) fprima(x,...)
          }
          mysample[i]<-sample$beta  
       }
       else
       {
         cat("\nError in sobroutine sample_...") 
         cat("\nifault=",sample$ifault,"\n")
       }
    }
  }
  else
  {
      cat("\nError in sobroutine initial_...") 
      cat("\nifault=",initial$ifault,"\n")      
  }
  return(mysample)
}

#--------------------------------------------------------------------------------------------
.onAttach <- function(library, pkg)
{
  Rv <- R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "3.1.2"))
    stop("This package requires R 3.1.2 or later")
  assign(".ars.home", file.path(library, pkg),
         pos=match("package:ars", search()))
  ars.version <- "0.5 (2014-12-03)"
  assign(".ars.version", ars.version, pos=match("package:ars", search()))
  if(interactive())
  {
    packageStartupMessage(paste("Package 'ars', ", ars.version, ". ",sep=""),appendLF=TRUE)
    packageStartupMessage("Type 'help(ars)' for summary information",appendLF=TRUE)
  }
  invisible()
}
