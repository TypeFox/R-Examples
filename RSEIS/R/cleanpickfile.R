'cleanpickfile'<-
  function(P)
{
  ############  clean up a pick file that has
  ############     old or non-conforming phase picks

  ########  if the pick has no station name it must be bpgus.  SEE EMPTYPickfile
  YPX=P$STAS
  w = which(YPX$name == "" | is.na(YPX$name) )

  if(length(w)<1) return(P)
  
  
  Y = deleteWPX(YPX, w)

  
  P$STAS =  as.list(Y)
  
  return(P)
  

}
