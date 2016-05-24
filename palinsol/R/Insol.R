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

 SIDERAL_YEAR = 365.25636 
 TROPIC_YEAR  = 365.24219876
 YEAR = TROPIC_YEAR 
 ## delibarate difference with 
 ## Berger et al.(2010), who consider the SIDERAL_YEAR

 
 # Returns daily mean incoming solar insolation after Berger (1978) 

 # eps  : obliquity (radians)
 # varpi: longitude of perihelion (radians)
 # e : eccentricity 
 # long : true solar longitude 
 #        (radians; = pi/2 for summer solstice )
 #                     pi   for autumn equinox  )
 #                     3* pi/2 for winter solstice )
 #                     0 for spring equinox )
 # lat : latitude
 # orbit is a list containing eps, varpi and e (see Ber78)
 # S0 : Total solar irradiance (default : 1365 W/m^2)
 # returns : daily mean incoming solar radiation at the top of the atmosphere
 #           in the same unit as S0. 
 Insol <- function (orbit,long=pi/2, lat=65*pi/180,S0=1365)
 {
  varpi <- NULL
  eps <- NULL
  ecc <- NULL
  for (i in names(orbit)) assign(i,orbit[[i]])
  nu <- long - varpi
  rho <- (1-ecc^2)/(1+ecc*cos(nu))
  sindelta <- sin(eps)*sin(long)
  cosdelta <- sqrt(1-sindelta^2)
  sinlatsindelta = sin(lat)*sindelta
  coslatcosdelta = cos(lat)*cosdelta
  cosH0 <- min(max(-1,-sinlatsindelta/coslatcosdelta),1)
  sinH0 <- sqrt(1-cosH0^2)
  H0 <- acos(cosH0)
  S0/(pi*rho^2)*(H0*sinlatsindelta+coslatcosdelta*sinH0)
  }
 
 ## time increment corresponding a tsl increment
 .dtdnu <- function (orbit,long=pi/2)
 {
  nu <- long - orbit['varpi']
  ecc <- orbit['ecc']
  xec <- ecc*ecc
  rho <- (1-xec)/(1+ecc*cos(nu))
  .dtdnu <- rho^2/sqrt(1.-xec)
 }
 # Provides an insolation times series
 # astrosol = astronomical solution (defaults to Ber78)
 # times = times in yr epoch 1950.0
 # ...   = any argument passed to Insol
 #InsolWrapper <- function(times=times,astrosol=ber78,...)
 #  sapply (times,function(tt)  Insol(astrosol(tt),...) )

 ## caloric_insolation
 ## integrated insolation over the 180 days receiving above median insolation
 calins <- function (orbit,lat=65*pi/180,...)
   {
    ins   <- sapply(seq(1:360)*pi/180, function(x) Insol(orbit,long=x, lat=lat,...))
    dt    <- sapply(seq(1:360)*pi/180, function(x) .dtdnu (orbit,long=x))
   
    is    <- sort(ins,decreasing=TRUE,index.return=TRUE)$ix
    cs    <- cumsum(dt[is])
    is    <- which(cs <= 180)           ## the 180 days whose cumulative length is half total
                                        ## year length, picking days by decreasing order of 
                                        ## insolation. 
    XCORR = 86.4 *  YEAR / 360
    sum(ins[is]*dt[is])* XCORR    ## result in kJ
   }


 ## integrated insolation over the 360 days receiving insolation above a threshold
 thrins <- function (lat=65*pi/180,orbit,threshold=400,...)
   {
    ins   <- sapply(seq(1:360)*pi/180, function(x) Insol(orbit,long=x, lat=lat,...))
    dt    <- sapply(seq(1:360)*pi/180, function(x) .dtdnu (orbit,long=x, ...))
    
    is    <- which(ins >= threshold)
    XCORR = 86.4 *  YEAR / 360
    sum(ins[is]*dt[is])* XCORR   ## result in kJ
   }



 ## time-integrated between two true solar longitude bounds
Insol_l1l2 <- function (orbit,l1=0,l2=2*pi,lat=65*pi/180,avg=FALSE,ell=TRUE,...)
   {
    # parameters: orbit : supplied by orbit calculator; e.g. : ber78 or ber90
    # l1 and l2 : longitudes bonds in radiants. Defaults to annual average.
    # discretize longitude intreval in N intervals
    # avg : supplies an average insolation
    # ell : defaults to TRUE, use elliptic integrals for calculation (much faster)
    #       rather than trapeze rule integral.  Currently incompatible
    #       with avg=TRUE (this can be fixed later)

     if (ell)
      ## use elliptic integrals
      {
          if (l1 < l2) 
          { 
          INT = W(lat, orbit['eps'], orbit['ecc'],l2,...) - 
                W(lat, orbit['eps'], orbit['ecc'], l1,...)} 
          else
          {
           INT = W(lat, orbit['eps'], orbit['ecc'], 2*pi,...) -
                 W(lat, orbit['eps'], orbit['ecc'], l1,...) +
                 W(lat, orbit['eps'], orbit['ecc'], l2,...)
           }
           if (avg) { 
                DT <- l2day(orbit,l2) - l2day(orbit,l1)
                if (DT <= 0) DT = DT+360.
                XCORR = 86.4 *  YEAR / 360
                INT = INT / (DT*XCORR) ## result in W/m2
                }
      }
      else
      ## integration using trapeze rule
      {
        Dl = ((l2-l1) %% (2*pi))
        if (Dl == 0) Dl=2*pi
        N =  1*ceiling(Dl* 180/pi)
        dl = Dl/N
        L  = l1+(0:N)*dl
        ins   <- sapply(L, function(x) Insol(orbit=orbit,long=x, lat=lat,...))
        dt    <- sapply(L, function(x) .dtdnu (orbit=orbit,long=x)) * 180./pi
   

        is    <- ins*dt
        XCORR = 86.4 *  YEAR / 360
        INT = (sum(is[2:N]) +  0.5 *is[1] +  0.5 *is[N+1]) * dl * XCORR  ## result in kJ
        
        if (avg) {
         DT =  (sum(dt[2:N]) + 0.5 * dt[1] + 0.5 * dt[N+1]) * dl * XCORR
         INT = INT / DT ## result in W/m2
         }
      }
      INT
    }

day2l  <- function (orbit,day)
   { 
    ## converts day to longitude.
    ## source : Berger 78, from Brower and Clemence.
    ## day using a 360-d calendar
    ## here :  day of spring equinox is in fact conventionally 80
    
    ecc = orbit['ecc']
    varpi= orbit['varpi']
    # definitions for speed-up
    xee= ecc*ecc
    xec = xee * ecc
    xse= sqrt(1.-xee)
    # true anomaly of vernal equinox point 
    xlp= - varpi 
    # mean anomaly  of the vernal equinox point
    lm =  xlp - 2.*((ecc/2 + xec/8)*(1+xse)*sin(xlp) - 
          xee/4*(0.5+xse)*sin(2*xlp) +
          xec/8*(1/3+xse)*sin(3*xlp) )

    # mean anomaly  of the point considered
    M = (day- 80) * pi/ 180. + lm 

    # TRUE anomaly of the point considered
    V = M + (2*ecc-xec/4.)*sin(M) +
          5/4*xee*sin(2*M) +
          13/12*xec*sin(3*M)

    # TRUE longitude of the point considered

    L = ( V + varpi ) %% (2*pi)
    L
   }

l2day <- function (orbit,l)
    ## converts true solar longitude to day
    ## source :  Brouwer and Clemence
   { 
    ecc = orbit['ecc']
    varpi= orbit['varpi']
    # definitions for speed-up
    xee= ecc*ecc
    xec = xee * ecc
    xse= sqrt(1.-xee)
    # true anomaly
    xlp= - varpi 
    # mean anomaly  of the vernal equinox point
    lm =  xlp - 2.*((ecc/2 + xec/8)*(1+xse)*sin(xlp) - 
          xee/4*(0.5+xse)*sin(2*xlp) +
          xec/8*(1/3+xse)*sin(3*xlp) )

    # true anomaly of the point considered
    V = (l + xlp) %% (2*pi)

    # mean anomaly  of the point considered
    M =  V - 2.*((ecc/2 + xec/8)*(1+xse)*sin(V) - 
          xee/4*(0.5+xse)*sin(2*V) +
          xec/8*(1/3+xse)*sin(3*V) )

    # anomaly in deg. elapsed between vernal equinox point and point

    DAY = ( 80 + (M - lm)  * 360.0/(2*pi) ) %% (360.0)
    DAY

   }


date_of_perihelion <- function(orbit)
 {
    ecc = orbit['ecc']
    varpi= orbit['varpi']
    # definitions for speed-up
    xee= ecc*ecc
    xec = xee * ecc
    xse= sqrt(1.-xee)
    # true anomaly
    xlp= - varpi 
    # mean anomaly  of the vernal equinox point
    lm =  xlp - 2.*((ecc/2 + xec/8)*(1+xse)*sin(xlp) - 
          xee/4*(0.5+xse)*sin(2*xlp) +
          xec/8*(1/3+xse)*sin(3*xlp) )

    # mean anomaly  of the point considered
    M =  0.

    # anomaly in deg. elapsed between vernal equinox point and perihelion passage

    DAY = ( 80 + (M - lm)  * 360.0/(2*pi) ) %% (360.0)
    names(DAY) <- 'day'
    DAY

   }



Insol_d1d2 <- function (orbit,d1,d2,lat=65*pi/180,avg=FALSE,...)
## as insol but given days rather than longitudes

   {
      l1 = day2l(orbit,d1)
      l2 = day2l(orbit,d2)
      Insol_l1l2(orbit,lat=lat,l1,l2,avg=avg,...)
    }  


# if (isTrue(getOption('debug')) && interactive())
# { t <- seq(-1e6,0,by=1e3)
#   F <- InsolWrapper(t)

Milankovitch <- function(orbit, S0=1365, lat=seq(-pi/2, pi/2, l=73), long=seq(0, 2*pi, l=145), deg=TRUE)

{
# returns a 'Milankovitch graph' of incoming solar irradiance
CONVERT=ifelse(deg, 180/pi, 1)

long = long[-length(long)]

M <- outer(long,lat,Vectorize(Insol,c("long","lat")),orbit=orbit,  S0=S0)
attr(M,  "class") <- 'Milankovitch'
attr(M,  "lat") <- lat * CONVERT
attr(M,  "long") <- long * CONVERT
attr(M,  "deg") <- deg
M
}

plot.Milankovitch <- function(x, months=TRUE, plot_function=image,...)
{
   long = attr(x, "long")
   if( ! attr(x,  "deg")) 
    {stop ('plot Milankovitch only when longitude in degrees')
    }
if (months)
 { 
   Col  = c(which(long >= (360-80)) , which(long < (360-80)))
   Col = c(Col, Col[1])
   Long = c(long, long[1]+360)
   MM = x[ Col,]
   plot_function(Long, attr(x, "lat"), MM, axes=FALSE,xlab='Month',ylab='Latitude',xaxs='i',yaxs='i',...)
    if(!exists("legend.only"))  {
   axis (1, at=seq(0,12)*30, labels=rep('',13))
   axis (3, at=seq(0,12)*30, labels=rep('',13))
   axis (1, at=seq(0,11)*30+15, labels=c('J','F','M','A','M','J','J','A','S','O','N','D'), tick=FALSE)}
 }
 else
 {   plot_function(attr(x, "long"), attr(x, "lat"), x, axes=FALSE, 
      xlab='True solar Longitude',ylab='Latitude',xaxs='i',yaxs='i',...)
      if(!exists("legend.only")) 
       {axis (1, at=seq(0,360,30))
       axis (3, at=seq(0,360,30))
       }
}
 if(!exists("legend.only")) 
 {
 axis(2, at=seq(-90,90,30), labels=c('90S','60S', '30S','Eq.','30N','60N','90S'))
 axis(4, at=seq(-90,90,30), labels=rep('',7))
 }
}


