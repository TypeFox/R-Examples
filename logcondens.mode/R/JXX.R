
J10.new <- function (x, y) 
{
    m <- length(x)
    z <- exp(x)
    d <- y - x
    ed <- exp(d)
    ## Note the function is asymmetric in x,y; it's "fine" if e^y is Inf as
    ## long as e^x is not (will get Inf returned, of course)
    LARGE.POS <- z==Inf  ## "POS" refers to x being (large) positive; Not good
    LARGE.NEG <- ed==Inf & !LARGE.POS ## "NEG" refers to x being (large) negative; OK
    MED <- (abs(d) > 0.01) & !LARGE.POS & !LARGE.NEG
    SMALL <- !(MED | LARGE.POS | LARGE.NEG)
    II <- (1:m)[LARGE.NEG | MED]
    z[II] <-  (exp(y[II]) - z[II] - d[II]*z[II])/(d[II]^2)
    II <- (1:m)[SMALL]
    z[II] <- z[II] * (1/2 + d[II] * (1/6 + d[II] * (1/24 + d[II] * 
        (1/120 + d[II]/720))))#by defn II, z[II] hasn't been modified yet
    II <- (1:m)[LARGE.POS]
    z[II] <- z[II] * (exp(d[II]) - 1 - d[II])/(d[II]^2) ##logcondens method of computing
    {if (any(is.nan(z))) {
      save(file="JXX-debug.rsav", x, y) ## returnign NaNs is basically an
                                        ## error, worth save() overhead
      warning("J10 is returning NaNs; see JXX-debug.rsav for values causing error.")
    }
    else if ((sum(LARGE.POS) > 0) || any(exp(y)==Inf)) ## 
      #print("Warning: J10 is computing Infs") ##any time Inf appears, could be creating errors ..
      ##save(file="JXX-debug.rsav", x, y) ##not worth the save() overhead;
      ##Inf may be reasonable/expected behavior
      warning("Warning: J10 is computing Infs") ##any time Inf appears, could
                                                ##be creating errors ..
   }
    return(z)
}


J10 <- function (x, y) 
{
    m <- length(x)
    z <- exp(x)
    d <- y - x
    LARGE <- (abs(d) > 0.01)
    SMALL <- !LARGE
    II <- (1:m)[LARGE]
    z[II] <-  (exp(y[II]) - z[II] - d[II]*z[II])/(d[II]^2)    
    II <- (1:m)[SMALL]
    z[II] <- z[II] * (1/2 + d[II] * (1/6 + d[II] * (1/24 + d[II] * 
        (1/120 + d[II]/720))))#by defn II, z[II] hasn't been modified yet
    return(z)
}


## J10.new <- function (x, y) 
## {
##     m <- length(x)
##     z <- exp(x)
##     d <- y - x
##     ed <- exp(d)
##     LARGE <- ed==Inf
##     MED <- (abs(d) > 0.01) & !LARGE ##this case necessary? 
##     SMALL <- !(MED | LARGE)
##     II <- (1:m)[LARGE]
##     z[II] <-  (exp(y[II]) - z[II] - d[II]*z[II])/(d[II]^2)
##     II <- (1:m)[MED]
##     z[II] <- z[II] * (ed[II] - 1 - d[II])/(d[II]^2)
##     II <- (1:m)[SMALL]
##     z[II] <- z[II] * (1/2 + d[II] * (1/6 + d[II] * (1/24 + d[II] * 
##         (1/120 + d[II]/720))))
##     return(z)
## }


J20 <- function (x, y) 
{
    m <- length(x)
    z <- exp(x)
    d <- y - x
    LARGE <- (abs(d)>0.02)
    SMALL <- !LARGE
    II <- (1:m)[LARGE]
    z[II] <- 2* (exp(y[II]) - z[II] - z[II]*d[II] - z[II]*d[II]^2/2 ) / (d[II]^3)
    ##II <- (1:m)[MED]
    ##z[II] <- 2 * z[II] * (exp(d[II]) - 1 - d[II] - d[II]^2/2)/(d[II]^3)
    II <- (1:m)[SMALL]
    z[II] <- z[II] * (1/3 + d[II] * (1/12 + d[II] * (1/60 + d[II] * 
        (1/360 + d[II]/2520))))
    return(z)
}


J20.new <- function (x, y) 
{

  m <- length(x)
  z <- exp(x)
  d <- y - x
  ed <- exp(d)
  ## Note the function is asymmetric in x,y; it's "fine" if e^y is Inf as
  ## long as e^x is not (will get Inf returned, of course)
  LARGE.POS <- z==Inf  ## "POS" refers to x being (large) positive; Not good
  LARGE.NEG <- ed==Inf & !LARGE.POS ## "NEG" refers to x being (large) negative; OK
  MED <- (abs(d) > 0.01) & !LARGE.POS & !LARGE.NEG
  SMALL <- !(MED | LARGE.POS | LARGE.NEG)
  II <- (1:m)[LARGE.NEG | MED]
  z[II] <- 2* (exp(y[II]) - z[II] - z[II]*d[II] - z[II]*d[II]^2/2 ) / (d[II]^3)
  II <- (1:m)[SMALL]
  z[II] <- z[II] * (1/3 + d[II] * (1/12 + d[II] * (1/60 + d[II] * 
                                                   (1/360 + d[II]/2520))))
  II <- (1:m)[LARGE.POS]
  z[II] <- 2 * z[II] * (exp(d[II]) - 1 - d[II] - d[II]^2/2)/(d[II]^3) ##logcondens method of computing
  {if (any(is.nan(z))) {
    save(file="JXX-debug.rsav", x, y) ## returnign NaNs is basically an
    ## error, worth save() overhead
    warning("J20 is returning NaNs; see JXX-debug.rsav for values causing error.")
  }
  else if ((sum(LARGE.POS) > 0) || any(exp(y)==Inf)) ## 
                                        #print("Warning: J10 is computing Infs") ##any time Inf appears, could be creating errors ..
    ##save(file="JXX-debug.rsav", x, y) ##not worth the save() overhead;
    ##Inf may be reasonable/expected behavior
    warning("Warning: J20 is computing Infs") ##any time Inf appears, could
   ##be creating errors ..
 }
  return(z)
}


## ##This version ~10%  slower it seems, using uniformly drawn x,y on 0,-100
## J20 <- function (x, y) 
## {
##     m <- length(x)
##     z <- exp(x)
##     d <- y - x
##     ed <- exp(d)
##     LARGE <- ed==Inf
##     MED <- (abs(d) > 0.02) & !LARGE
##     SMALL <- !(MED | LARGE)
##     II <- (1:m)[LARGE]
##     z[II] <- 2* (exp(y[II]) - z[II] - z[II]*d[II] - z[II]*d[II]^2/2 ) / (d[II]^3)
##     II <- (1:m)[MED]
##     z[II] <- 2 * z[II] * (exp(d[II]) - 1 - d[II] - d[II]^2/2)/(d[II]^3)
##     II <- (1:m)[SMALL]
##     z[II] <- z[II] * (1/3 + d[II] * (1/12 + d[II] * (1/60 + d[II] * 
##         (1/360 + d[II]/2520))))
##     return(z)
## }

J11.new <-   function (x, y) 
{
  m <- length(x)
  z <- exp(x)
  d <- y - x
  ed <- exp(d)
  ## Note the function is asymmetric in x,y; it's "fine" if e^y is Inf as
  ## long as e^x is not (will get Inf returned, of course)
  LARGE.POS <- z==Inf  ## "POS" refers to x being (large) positive; Not good
  LARGE.NEG <- ed==Inf & !LARGE.POS ## "NEG" refers to x being (large) negative; OK
  MED <- (abs(d) > 0.01) & !LARGE.POS & !LARGE.NEG
  SMALL <- !(MED | LARGE.POS | LARGE.NEG)
  II <- (1:m)[LARGE.NEG | MED]
  z[II] <- (d[II]*(exp(y[II])+z[II]) - 2*(exp(y[II])-z[II])) / (d[II]^3)

  II <- (1:m)[SMALL]
  z[II] <- z[II] * (1/6 + d[II] * (1/12 + d[II] * (1/40 + d[II] * 
                                                   (1/180 + d[II]/1008))))
  
  II <- (1:m)[LARGE.POS]
  z[II] <- z[II] * (d[II] * (exp(d[II]) + 1) - 2 * (exp(d[II]) - 
                                                    1))/(d[II]^3)
  
  {if (any(is.nan(z))) {
    save(file="JXX-debug.rsav", x, y) ## returnign NaNs is basically an
    ## error, worth save() overhead
    warning("J11 is returning NaNs; see JXX-debug.rsav for values causing error.")
  }
  else if ((sum(LARGE.POS) > 0) || any(exp(y)==Inf)) ## 
                                        #print("Warning: J10 is computing Infs") ##any time Inf appears, could be creating errors ..
    ##save(file="JXX-debug.rsav", x, y) ##not worth the save() overhead;
    ##Inf may be reasonable/expected behavior
    warning("Warning: J11 is computing Infs") ##any time Inf appears, could
   ##be creating errors ..
 }
  return(z)
}


J11 <-   function (x, y) 
{
  m <- length(x)
  z <- exp(x)
  d <- y - x
  ##ed <- exp(d)
  ##LARGE <- ed==Inf
  ##MED <- (abs(d) > 0.02) & !LARGE
  ##SMALL <- !(MED | LARGE)
  LARGE <- (abs(d)>.02)
  SMALL <- !LARGE
  II <- (1:m)[LARGE]
  z[II] <- (d[II]*(exp(y[II])+z[II]) - 2*(exp(y[II])-z[II])) / (d[II]^3)
  ##II <- (1:m)[MED]
  ##z[II] <- z[II] * (d[II] * (exp(d[II]) + 1) -
  ##                  2 * (exp(d[II]) -  1))/(d[II]^3)
  II <- (1:m)[SMALL]
  z[II] <- z[II] * (1/6 + d[II] * (1/12 + d[II] * (1/40 + d[II] * 
                                                   (1/180 + d[II]/1008))))
  return(z)
}

## J11 <-   function (x, y) 
## {
##   m <- length(x)
##   z <- exp(x)
##   d <- y - x
##   ed <- exp(d)
##   LARGE <- ed==Inf
##   MED <- (abs(d) > 0.02) & !LARGE
##   SMALL <- !(MED | LARGE)
##   II <- (1:m)[LARGE]
##   z[II] <- (d[II]*(exp(y[II])+z[II]) - 2*(exp(y[II])-z[II])) / (d[II]^3)
##   II <- (1:m)[MED]
##   z[II] <- z[II] * (d[II] * (exp(d[II]) + 1) -
##                    2 * (exp(d[II]) -  1))/(d[II]^3)
##   II <- (1:m)[SMALL]
##   z[II] <- z[II] * (1/6 + d[II] * (1/12 + d[II] * (1/40 + d[II] * 
##                                                    (1/180 + d[II]/1008))))
##   return(z)
## }



J00.new <- function (x, y, v = 1) 
{
  m <- length(x)
  z <- exp(x)
  d <- y - x
  ed <- exp(d)
  ## Note the function is asymmetric in x,y; it's "fine" if e^y is Inf as
  ## long as e^x is not (will get Inf returned, of course)
  LARGE.POS <- z==Inf  ## "POS" refers to x being (large) positive; Not good
  LARGE.NEG <- ed==Inf & !LARGE.POS ## "NEG" refers to x being (large) negative; OK
  MED <- (abs(d) > 0.01) & !LARGE.POS & !LARGE.NEG
  SMALL <- !(MED | LARGE.POS | LARGE.NEG)
  II <- (1:m)[LARGE.NEG | MED]
  z[II] <- (exp(v*y[II]+(1-v)*x[II]) - z[II]) / d[II]
  II <- (1:m)[SMALL]
  z[II] <- z[II] * (v + d[II] * (v/2 + d[II] * (v/6 + d[II] * 
                                                (v/24 + d[II] * v/120))))
  II <- (1:m)[LARGE.POS]
  z[II] <- z[II] * (exp(v * d[II]) - 1)/d[II]
  {if (any(is.nan(z))) {
    save(file="JXX-debug.rsav", x, y) ## returnign NaNs is basically an
    ## error, worth save() overhead
    warning("J00 is returning NaNs; see JXX-debug.rsav for values causing error.")
  }
  else if ((sum(LARGE.POS) > 0) || any(exp(y)==Inf)) ## 
                                        #print("Warning: J10 is computing Infs") ##any time Inf appears, could be creating errors ..
    ##save(file="JXX-debug.rsav", x, y) ##not worth the save() overhead;
    ##Inf may be reasonable/expected behavior
    warning("Warning: J00 is computing Infs") ##any time Inf appears, could
   ##be creating errors ..
 }
  return(z)
}



J00 <- function (x, y, v = 1) 
{
  m <- length(x)
  z <- exp(x)
  d <- y - x
  ed <- exp(v*d)
  LARGE <- (abs(d)>0.005)
  SMALL <- !LARGE
  II <- (1:m)[LARGE]
  z[II] <- (exp(v*y[II]+(1-v)*x[II]) - z[II]) / d[II]
  II <- (1:m)[SMALL]
  z[II] <- z[II] * (v + d[II] * (v/2 + d[II] * (v/6 + d[II] * 
                                                (v/24 + d[II] * v/120))))
  return(z)
}

