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

## The present file is an implementation of 
## Andre Berger, Marie-France Loutre, and Qiuzhen Yin, 
## Total irradiation during any time interval of the year using elliptic integrals,  
## Quaternary Science Reviews, 29, 1968 - 1982  2010


# ------------------------------------------------------------------
# R Code developed for R version 2.15.2 (2012-10-26) -- "Trick or Treat"
# ------------------------------------------------------------------ 



SIDERAL_YEAR = 365.25636 
 TROPIC_YEAR  = 365.24219876
 YEAR = SIDERAL_YEAR


W  <- function (phi, eps, ecc, lambda,S0=1365,n=3)
   {
 
    # require(gsl) ## gnu scientific library ported by Robin Hankin. Thank you Robin !!! 
    pi2=pi/2
    H00   = pi2
    seps = sin(eps)
    ceps = cos(eps)
    sphi = sin(phi)
    cphi = cos(phi)
    tphi = sphi/cphi

    eq34 <- function (lambda)
    {
    slambda = sin(lambda)
    clambda = cos(lambda)
    k   = seps/cphi
    sesl   = seps*slambda
    tdelta = sesl / sqrt(1-sesl*sesl)
    H0    = acos(-tphi*tdelta)
    Flk   = ellint_F(lambda,k)
    Elk   = ellint_E(lambda,k)

    sphi*seps*(H00-clambda*H0) + cphi*Elk +
        sphi*tphi*Flk - sphi*tphi*ceps*ceps*
        ellint_P(lambda, k, -seps*seps)
    }
#

     eq40  <- function(lambda)
    {
    ## the max(-1, min(1, 
    ## is to account for a numerical artefact when lambda = lambda1,2,3,4

    slambda = sin(lambda)
    clambda = cos(lambda)

    k =  seps/cphi
    k1 = cphi/seps
    psi  = asin(max(-1,min(1,k*slambda)))

    if (clambda < 0) psi = pi-psi
     
    sesl   = seps*slambda
    tdelta = sesl / sqrt(1-sesl*sesl)
    H0    = acos(max(-1,min(1,-tphi*tdelta)))
    Fpk   = ellint_F(psi,k1)
    Epk   = ellint_E(psi,k1)
    Pipk  = ellint_P(psi,k1, -cphi*cphi)

    ( sphi*seps*(H00-clambda*H0) + seps*Epk +
        ceps* ceps/seps * Fpk - sphi*sphi*ceps*ceps/seps*
        Pipk)
    }
#

    eq38  <- function(lambda) { - pi * sphi*seps*cos(lambda) }


    T   = YEAR * 0.086400 * 1000
    xes = sqrt(1-ecc*ecc)
    W0  = S0*T/(2*pi*pi*xes)

    if (phi >= (pi2-eps) | phi <= -(pi2-eps) )
    {
    ## above polar circle 
    lambda1 = (asin(cphi/seps) )
    lambda2 = pi - lambda1
    lambda3 = pi + lambda1
    lambda4 = 2*pi - lambda1
  
    WW=0

    if (lambda > 0)         WW = WW + eq40(min(lambda,lambda1))
    if (lambda > lambda2)   WW = WW + eq40(min(lambda,lambda3)) - eq40(lambda2)
    if (lambda > lambda4)   WW = WW + eq40(lambda) - eq40(lambda4)

    if (phi >= (pi2-eps) ) {  ## northern hemisphere
      if (lambda > lambda1)   WW = WW + eq38(min(lambda,lambda2)) - eq38(lambda1)
      } else 
    { ## Southern hemisphere
      if (lambda > lambda3)   WW = WW + eq38(min(lambda,lambda4)) - eq38(lambda3)
    }

    WW = W0*WW 
    } else ## outside polar circle
    { WW = W0 * eq34(lambda) }
   WW
   }

