# Copyright (c) 2012 Michel Crucifix <michel.crucifix@uclouvain.be>

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject

# the following conditions:

# The above copyright notice and this permission notice shall be
# incluudedin all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND INFRINGEMENT
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR

# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# When using this package for actual applications, always
# cite the authors of the original insolation solutions 
# Berger, Loutre and/or Laskar, see details in man pages


# ------------------------------------------------------------------
# R Code developed for R version 2.15.2 (2012-10-26) -- "Trick or Treat"
# ------------------------------------------------------------------ 

 ## calculates astronomical elements according to the solutions given below
 

.BER78 <- new.env()
.BER90 <- new.env()
.LA04  <- new.env()

astro <- function(t,solution=ber78,degree=FALSE) {solution(t,degree)}

 ## Calculates climate orbital elements according to the algorithm given in A. Berger (1978)
 # Berger, A. L. (1978).  Long-term variations of daily insolation and Quaternary climatic changes, 
 # J. Atmos. Sci., 35(), 2362-2367

 # This solution is valid for + / - 1e6 years. 
 # uses MBCS_BER78 provided by A. Berger, available on
 # ftp://ftp.astr.ucl.ac.be/pub/berger/berger78/INSOL.IN

## attach table to ber78 function
## Input :  t = time expressed in yr after 1950.0 (reference epoch)


.BER78 <- new.env()
.BER90 <- new.env()
.LA04  <- new.env()

ber78 <- function(t,degree=FALSE)
{
  if (!exists(".loaded", envir=.BER78))
  { 
     message('load BER78data')
     data('BER78', envir=.BER78)
     assign('.loaded',  TRUE,  envir = .BER78 )
  }

  
  psibar<- 50.439273/60./60. * pi/180 
  estar <- 23.320556
  e0    <- 0.028707
  zeta  <- 3.392506 * pi/180.
  twopi <- 2*pi
 
  sectorad <- pi/(180*60.*60.)

  M <- .BER78$Table4$Amp
  g <- .BER78$Table4$Rate*sectorad
  b <- .BER78$Table4$Phase*pi/180
  F <- .BER78$Table5$Amp*sectorad
  fp <- .BER78$Table5$Rate*sectorad
  d <- .BER78$Table5$Phase*pi/180.
  A <- .BER78$Table1$Amp/60./60.
  f <- .BER78$Table1$Rate*sectorad
  phi<- .BER78$Table1$Phase*pi/180.

  ## Obliquity

  eps <- estar + sum(A*cos(f*t+phi))
  epsp <- estar + sum(A*sin(f*t+phi))

  esinpi <- sum(M*sin(g*t+b))
  ecospi <- sum(M*cos(g*t+b))
  psi    <- psibar*t + zeta + sum(F*sin(fp*t +d))

  e <- sqrt(esinpi^2+ecospi^2)
  Pi <-atan(esinpi/ecospi)+pi*(ecospi<0)
  eps <- eps * pi/180.
  epsp <- epsp * pi/180.
  varpi <- (Pi+psi+pi) %% (twopi)

  if (degree) {rad2deg <- 180/pi
               eps <- eps*rad2deg
               varpi <- varpi*rad2deg}

  c(eps=eps,ecc=e,varpi=varpi,epsp=epsp)
}
 ## Calculates climate orbital elements according to the algorithm given in A. Berger (1978)
 # Berger, A. L. (1978).  Long-term variations of daily insolation and Quaternary climatic changes, 
 # J. Atmos. Sci., 35(), 2362-2367

 # This solution is valid for + / - 3e6 years. 

## attach table to ber90 function
## Input :  t = time expressed in yr after 1950.0 (reference epoch)

ber90 <- function(t,degree=FALSE)
{
  if (!exists(".loaded", envir=.BER90))
  {
     message("load BER90 data")
     data('BER90', envir=.BER90)
     assign('.loaded',  TRUE,  envir = .BER90 )
  }
 
  psibar<- 50.41726176/60./60. * pi/180 
  estar <- 23.33340950
  zeta  <- 1.60075265 * pi/180.
  twopi <- 2*pi
 
  sectorad <- pi/(180*60.*60.)

  M <- .BER90$Table4$Amp
  g <- .BER90$Table4$Rate*sectorad
  b <- .BER90$Table4$Phase*pi/180
  F <- .BER90$Table5$Amp*sectorad
  fp <- .BER90$Table5$Rate*sectorad
  d <- .BER90$Table5$Phase*pi/180.
  A <- .BER90$Table1$Amp/60./60.
  f <- .BER90$Table1$Rate*sectorad
  phi<- .BER90$Table1$Phase*pi/180.

  ## Obliquity

  eps <- estar + sum(A*cos(f*t+phi))
  epsp <- estar + sum(A*sin(f*t+phi))

  esinpi <- sum(M*sin(g*t+b))
  ecospi <- sum(M*cos(g*t+b))
  psi    <- psibar*t + zeta + sum(F*sin(fp*t +d))

  e <- sqrt(esinpi^2+ecospi^2)
  Pi <-atan(esinpi/ecospi)+pi*(ecospi<0)
  eps <- eps * pi/180.
  epsp <- epsp * pi/180.
  varpi <- (Pi+psi+pi) %% (twopi)

  if (degree) {rad2deg <- 180/pi
               eps <- eps*rad2deg
               varpi <- varpi*rad2deg}

  c(eps=eps,ecc=e,varpi=varpi,epsp=epsp)
}

 # This solution is valid for 50e6 years

## Input :  t = time expressed in yr after 1950.0 (reference epoch)

la04 <- function(t,degree=FALSE)
{
  if (!exists(".loaded", envir=.LA04))
  {
     data('LA04', envir=.LA04)
     message("LA04 data loaded")
     assign('.loaded',  TRUE,  envir = .LA04 )
  }
 
  tka = t/1000.
  if (tka>0) 
  {
   local(
   {
    F <-  floor(tka)
    ORB <<- .LA04$la04future[F+1, ]
    if (! (tka == F)) {
      D  <- tka - floor(tka)
      diff <- .LA04$la04future[F+2, ] - ORB
      # note : if the diff in varpi is greater than pi,
      # this probably means that we have skipped 2*pi, 
      # so we need to correct accordingly
      if (diff$varpi > pi) diff$varpi = diff$varpi - 2*pi
      if (diff$varpi < -pi) diff$varpi = diff$varpi + 2*pi
      #
      ORB <<- ORB + D*diff
    }
   })}
   else
   {
   local(
   {
    F <-  floor(tka)
    ORB <<- .LA04$la04past[-F+1, ]
    if (! (tka == F)) {
      D  <- tka - F

      diff <- .LA04$la04past[-F, ] - ORB
      # note : if the diff in varpi is greater than pi,
      # this probably means that we have skipped 2*pi, 
      # so we need to correct accordingly
      if (diff$varpi > pi) diff$varpi = diff$varpi - 2*pi
      if (diff$varpi < -pi) diff$varpi = diff$varpi + 2*pi
      #
      ORB <<- ORB + D*diff
    }
   })}
  
  if (degree) {rad2deg <- 180/pi
               ORB['eps'] <- ORB['eps']*rad2deg
               ORB['varpi'] <- ORB['varpi']*rad2deg}


   # must return a array (0.92 -> 0.93)
   names <- c('eps','ecc','varpi')
   OUT = as.numeric(ORB[names]) ; names(OUT) <- names
   OUT

   }

precession <- function(t,solution=ber78)
##  as astro, but returns only precession parameter e sin (varpi)
{ 
  O <- astro(t,solution, degree=FALSE)
  O['ecc'] * sin (O['varpi'])
}

obliquity <- function(t,solution=ber78,degree=FALSE)
##  as astro, but returns only obliquity
{ astro(t,solution, degree=degree)['eps'] }
