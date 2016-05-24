gif <- function(data,gifset)
{
  famsize<-dim(data)[1]
  giflen<-length(gifset)
  gifval<-0
  z<-.C("gif",data=as.integer(t(data)),famsize=as.integer(famsize),
        gifset=as.integer(array(gifset)),giflen=as.integer(giflen),
        gifval=as.double(gifval),PACKAGE="gap")

  list(gifval=z$gifval)
}
