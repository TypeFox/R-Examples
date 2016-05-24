`fits.isoreg` <-
function(iso, x0){
  o = iso$o
  if(is.null(o)) o=1:length(x)
  x = iso$x[o]
  y = iso$yf
  ind = cut(x0, breaks=x, labels=FALSE, include.l=TRUE)
  fits = sapply(seq(along=x0),
    function(i){
      j = ind[i]
      y[j] + (x0[i] - x[j])*(y[j+1] - y[j])/(x[j+1] - x[j])
      }
    )
  return(fits)
  }

