neighbors <- function(x,y,...){
  nbh=neighborhood(x,...)
  levs=rep(1:nbh$nu,nbh$size.nbh)
  nbh.list <- split(y[nbh$neighbors],levs)
  list(nbh=nbh,list=nbh.list)
}
