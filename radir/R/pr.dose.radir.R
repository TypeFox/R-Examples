pr.dose.radir <-
function(object, lod = 0, upd = object[[2]][length(object[[2]])])
{
  if (class(object)!="dose.radir") stop("Wrong object")
  
  cd  <- approxfun(object[[2]], object[[1]], rule=2)
  if(lod < 0) lod <- 0
  if(upd > object[[2]][length(object[[2]])]) upd <- object[[2]][length(object[[2]])] 
  res <- integrate(cd, lower=lod, upper=upd)$value
  if (res > 1) return(1)
  return(res)
}