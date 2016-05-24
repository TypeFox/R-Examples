`restrictY` <-
function(dd,yy,rr,level,verbose=0) {
  if (sum(yy^2, na.rm = TRUE) == 0) return(yy)
  switch(level,
  "nominal"=return(nominalY(dd,yy,rr,verbose=verbose)),
  "ordinal"=return(ordinalY(dd,yy,rr,verbose=verbose)),
  "numerical"=return(numericalY(dd,yy,rr,verbose=verbose)),
  "polynomial"=return(polynomialY(dd,yy,rr,verbose=verbose)))
}

