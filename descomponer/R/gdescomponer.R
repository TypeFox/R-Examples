gdescomponer <- function(y,freq,type,year,q) {
  # Author: Francisco Parra Rodriguez
  # http://econometria.wordpress.com/2013/08/21/estimation-of-time-varying-regression-coefficients/ 
  serie <- descomponer (y,freq,type)
  TdsT <- c(serie$datos$TDST)
  Td <- c(serie$datos$TD)
  sT <- c(serie$datos$ST)
  TDST <- ts(TdsT,frequency=freq,start = c(year,q))
  TD <- ts(Td,frequency=freq,start = c(year,q))
  ST <- ts(sT,frequency=freq,start = c(year,q))
  par(mfrow=c(3,1))
  plot (TDST)
  plot (TD)
  plot (ST)
}