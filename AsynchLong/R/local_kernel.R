local_kernel <- function(t, h, kType){

  if( kType == "epan" ) {

    kt <- epanechnikov(t/h)

  } else if( kType == "uniform" ) {

    kt <- uniform(t/h)

  } else if( kType == "gauss" ) {

    kt <- gauss(t/h)

  } else {

    stop("unsupported kernel")

  }

  return(kt/h)
}

epanechnikov <- function(t){

  tst <- (-1.0 <= t) & (t <= 1.0 )

  kt <- 0.75*( 1.0 - t*t )
  kt[ !tst ] <- 0.0

  return(kt)
}

uniform <- function(t){

  tst <- (-1.0 <= t) & (t <= 1.0 )
  kt <- t
  kt[  tst ] <- 0.5
  kt[ !tst ] <- 0.0

  return(kt)
}

gauss <- function(t){

  kt <- exp(-t*t*0.5)/sqrt(2.0*pi)

  return(kt)

}
