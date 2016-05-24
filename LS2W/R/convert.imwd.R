convertimwd <-
function (imwd, ...) 
{
  if (imwd$type != "station")
  stop("Object to convert must be of type \"station\" ")

  J <- imwd$nlevels

  n <- c(2^J, 2^J)   
  dummy <- matrix(0,n[1],n[2])
  tmpwst2D <- wst2D(dummy, filter.number = imwd$filter$filter.number, family 
=imwd$filter$family)
  tmpwst2D$date <- imwd$date

  arrvec <- getarrvec(J+1, sort = FALSE)

  for (lev in (J - 1):0) {	 
  o<-arrvec[,J-lev] 						
  packmat<-packetj(imwd,lev,o) 					
  tmpwst2D$wst2D[lev+1,,]<-packmat
}

 tmpwst2D
}

