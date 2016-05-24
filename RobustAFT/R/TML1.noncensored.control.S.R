TML1.noncensored.control.S <- function(tlo=0.0001, mxf=50, mxs=50, ntm=50, tls=1e-6, h=100) {
  if (!is.numeric(tlo) || tlo <= 0)
        stop("value of 'tlo' must be > 0")
  if (!is.numeric(mxf) || mxf < 2)
        stop("value of 'mxf' must be an integer > 1")
  if (!is.numeric(mxs) || mxs < 1)
        stop("value of 'mxs' must be an integer > 0")
  if (!is.numeric(tls) || tls <= 0)
        stop("value of 'tls' must be > 0")
  if (!is.numeric(h) || h < 6)
        stop("'h' must be an integer > 5 ")
 list(tlo=tlo,mxf=mxf,mxs=mxs,ntm=as.numeric(ntm),tls=tls,h=h)}


