gammaIntegral <- function(total, params, offset=50, minx=1) {

  # require('gsl') # for hyperg_1F1 in the conditional expectation of the signal

  ## this bit is the most obvious "farm me out to C++" piece of all...
  if(length(total) > 1) return(sapply(total, gammaIntegral, params=params))
  
  g = params[1]
  a = params[2]
  d = params[3]
  b = params[4]
  bg.mean = d * b 
  bg.sd = sqrt( d * b * b )
  ## cat('total =',total,'... bg.mean =',bg.mean,'... bg.sd =',bg.sd,"\n")
  if(total > ( bg.mean + ( 3 * bg.sd ) )) {
    return(pmax(total - bg.mean, minx))
  } else {
    ## FIXME: need to write this as a function object for C++ to integrate it
    res = try(
      integrate( 
        function(x) {
          (exp(x*((1/b)-(1/a)))*(total**(1-g-d))*((total-x)**(d-1))*(x**(g-1)))*
          (1/(beta(g,d)*hyperg_1F1(g, g+d, total*((1/b)-(1/a)), strict=F)))*x 
        }, 
        0, total
      )$value # else will return a list with value, abs.error, subdivisions, ...
    ) # i.e., integrate(PrSignalGivenTotal, /* from */ 0, /* to */ total);
    if(class(res) == 'try-error') {
      return(pmax(total-bg.mean, minx)+offset)
    } else {
      return(pmax(res, minx)+offset)
    }
  }

}
