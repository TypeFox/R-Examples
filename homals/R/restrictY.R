`restrictY` <-
function(d,y,r,level,verbose=0) {
  if (sum(y^2, na.rm = TRUE) == 0) return(y)
  switch(level,
  "nominal"=return(nominalY(d,y,r,verbose=verbose)),
  "ordinal"=return(ordinalY(d,y,r,verbose=verbose)),
  "numerical"=return(numericalY(d,y,r,verbose=verbose)),
  "polynomial"=return(polynomialY(d,y,r,verbose=verbose)))
}

