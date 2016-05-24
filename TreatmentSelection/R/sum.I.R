sum.I <-
function(yy,FUN,Yi){

  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}

  pos <- rank(c(yy,Yi),ties.method='max')[1:length(yy)]-rank(yy,ties.method='max')

  if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos
  
  return(pos)
}
