exportNetLogo <- function(object,objectname,....)
{
  if (is.character(object))
  {
    object <- paste0('"',object,'"')
  }
  if (is.vector(object))
  {
    res=paste(paste0("[",paste(object, collapse = " "),"]"),collapse="\n")
  } else if (is.matrix(object))
  {
    res=paste("[\n",paste(apply(object,1,function(s)paste0("[",paste(s, collapse = " "),"]")),collapse="\n"),"\n]")
  } else stop("Object not supported")
  write.table(res,file=paste0(objectname,".txt"),row.names=FALSE, col.names=FALSE,quote=FALSE)
  return(res)
}


# zo maken dat de weiadj en thresholds als aparte objecten eruit komen die als txt-file opgeslagen kunnen worden.
# Ook zo maken dat deze functie individueel opgeroepen kan worden. Bijv om de symptoomnamen in file te krijgen.

