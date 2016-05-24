get.slepians<-function(npoints=900, nwin=5, npi=3)
{
  ######   return nwin npi slepian tapers for an npoints vector
  if(missing(npoints)) npoints=900
  if(missing(nwin)) nwin=5
  if(missing(npi)) npi=3

  tapers = vector(length=nwin*npoints)
  ary = .C("CALL_slepian",  PACKAGE = "RSEIS",
    as.integer(npoints), as.integer(nwin), as.double(npi), as.double(tapers) )

  kout = ary[[4]]

  mout =  matrix(        kout, ncol=nwin)

  invisible(mout) 
}
