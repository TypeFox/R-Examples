"statTn" <-
 function(uv, cop=NULL, para=NULL, p=2, proot=FALSE, ...) {
   proot  <- ifelse(proot, 1/p, 1)
   parcop <- COP(uv[,1], uv[,2], cop=cop, para=para, ...)
   # The empirical is called 2nd because it can be CPU intensive
   # and if the parameteric call is going to error, let us have it
   # do it first and not waste users time.
   empcop <- EMPIRcopdf(para=uv, ...)$empcop
   (sum((parcop - empcop)^p))^proot
}

