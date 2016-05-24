symmetricKnots <- function(knotsFromFitobject, localLow, localUp, refpoint, lrthreshold) {
  #### reflexion if necessary, of the knotsFromFitobject, not of the keypts; note that reflected points are possibly redundant
  fitobjectSym <- constrSymOrReflWrapper(points=knotsFromFitobject, ,
                                         logLcheck=lrthreshold, refpoint=refpoint, 
                                         lower=localLow, upper=localUp) ## was reflexion of knots through mirrorface, now symmetry through refpoint...
  return(fitobjectSym)
} ## end def addSymmetricKnots
