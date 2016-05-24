fixpaste <-
function(x,fix)
# fixpaste() - paste elements of vector x onto a fixed character
{
  xfix <- rep(0,length(x))
  for(i in 1: length(x)){
      xfix[i] <- paste(x[i],":",fix,sep="", collapse=NULL)
  }
  return(xfix)
}
