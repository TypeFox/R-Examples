## Copyright (C) 2001 Paulo Neis
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
## 

## N-ellip 0.2.1
##usage: [Zz, Zp, Zg] = ellip(n, Rp, Rs, Wp, stype,'s')
##
## Generate an Elliptic or Cauer filter (discrete and contnuious).
## 
## [b,a] = ellip(n, Rp, Rs, Wp)
##  low pass filter with order n, cutoff pi*Wp radians, Rp decibels 
##  of ripple in the passband and a stopband Rs decibels down.
##
## [b,a] = ellip(n, Rp, Rs, Wp, 'high')
##  high pass filter with cutoff pi*Wp...
##
## [b,a] = ellip(n, Rp, Rs, [Wl, Wh])
##  band pass filter with band pass edges pi*Wl and pi*Wh ...
##
## [b,a] = ellip(n, Rp, Rs, [Wl, Wh], 'stop')
##  band reject filter with edges pi*Wl and pi*Wh, ...
##
## [z,p,g] = ellip(...)
##  return filter as zero-pole-gain.
##
## [...] = ellip(...,'s')
##     return a Laplace space filter, W can be larger than 1.
## 
## [a,b,c,d] = ellip(...)
##  return  state-space matrices 
##
## References: 
##
## - Oppenheim, Alan V., Discrete Time Signal Processing, Hardcover, 1999.
## - Parente Ribeiro, E., Notas de aula da disciplina TE498 -  Processamento 
##   Digital de Sinais, UFPR, 2001/2002.
## - Kienzle, Paul, functions from Octave-Forge, 1999 (http://octave.sf.net).
##
## Author: Paulo Neis <p_neis@yahoo.com.br>
## Modified: Doug Stewart Feb. 2003


ellip <- function(n, ...) UseMethod("ellip")

ellip.FilterOfOrder <- function(n, Rp = n$Rp, Rs = n$Rs, W = n$Wc, type = n$type, ...)
  ellip(n$n, Rp, Rs, W, type, ...)

ellip.default <- function(n, Rp, Rs, W, type = c("low", "high", "stop", "pass"), plane = c("z", "s"), ...) {

  type <- match.arg(type)
  plane <- match.arg(plane)

  ## interpret the input parameters
  if (!(length(n)==1 && n == round(n) && n > 0))
    stop("filter order n must be a positive integer")

  stop <- type == "stop" || type == "high"
  digital <- plane == "z"

  if (length(W) != 1 && length(W) != 2)
    stop("frequency must be given as w0 or c(w0, w1)")

  if (digital && !all(W >= 0 & W <= 1))
    stop("critical frequencies must be in (0 1)")
  else if (!digital && !all(W >= 0))
    stop("critical frequencies must be in (0 inf)")

  if (Rp < 0)
    stop("passband ripple must be positive decibels")
  if (Rs < 0)
    stop("stopband ripple must be positive decibels")

  
  ## Prewarp to the band edges to s plane
  if (digital) {
    T <- 2       # sampling frequency of 2 Hz
    W <- 2 / T*tan(pi * W / T)
  }

  ##Generate s-plane poles, zeros and gain
  ZPG <- ncauer(Rp, Rs, n)

  ## s-plane frequency transform
  ZPG <- sftrans(ZPG, W = W, stop = stop)

  ## Use bilinear transform to convert poles to the z plane
  if (digital)
     ZPG <- bilinear(ZPG, T = T)

  as.Arma(ZPG)
}

ncauer <- function(Rp, Rs, n)  { #  [zer, pol, T0]

  ellip_ws <- function(n, rp, rs) {
    ellip_ws_min <- function(kl) {
      int <- ellipke(c(kl, 1-kl))$k
      ql <- int[1]
      q <- int[2]
      abs((ql/q) - x)
    }
    kl0 <- ((10^(0.1*rp) - 1) / (10^(0.1*rs) - 1))
    k0 <- 1 - kl0
    int <- ellipke(c(kl0, k0))$k
    ql0 <- int[1]
    q0 <- int[2]
    x <- n * ql0 / q0
#    kl <- optim(ellip_ws_min, method = "L-BFGS-B", lower = eps, upper = 1-eps)
    kl <- optimize(ellip_ws_min, 
        interval = c(.Machine$double.eps, 1-.Machine$double.eps))$minimum
    ws <- sqrt(1/kl)
    ws
  }

  ## Cutoff frequency = 1:
  wp <- 1
  
  ## Stop band edge ws:
  ws <- ellip_ws(n, Rp, Rs)
  
  k  <- wp / ws
  k1 <- sqrt(1-k^2)
  q0 <- (1/2)*((1-sqrt(k1))/(1+sqrt(k1)))
  q <-  q0 + 2*q0^5 + 15*q0^9 + 150*q0^13

  ####Filter order maybe this, but not used now:
  ##D <-  (10^(0.1*Rs)-1)/(10^(0.1*Rp)-1)
  ##n<-ceil(log10(16*D)/log10(1/q))

  l <- (1/(2*n))*log((10^(0.05*Rp)+1)/(10^(0.05*Rp)-1))
  sig01 <- 0
  sig02 <- 0
  for (m in 0:30) {
    sig01 <- sig01+(-1)^m * q^(m*(m+1)) * sinh((2*m+1)*l)
  }
  for (m in 1:30) {
    sig02 <- sig02+(-1)^m * q^(m^2) * cosh(2*m*l)
  } #
  sig0 <- abs((2*q^(1/4)*sig01)/(1+2*sig02))

  w <- sqrt((1+k*sig0^2)*(1+sig0^2/k))
  
  r <- (n - (n %% 2))/2
  
  wi <- matrix(0, 1,r)
  for (ii in 1:r) {
    mu <- ii - (1-(n %% 2))/2
    soma1 <- 0
    for (m in 0:30) {
      soma1 <- soma1 + 2*q^(1/4) * ((-1)^m * q^(m*(m+1)) * sin(((2*m+1)*pi*mu)/n))
    }
    soma2 <- 0
    for (m in 1:30) {
      soma2 <- soma2 + 2*((-1)^m * q^(m^2) * cos((2*m*pi*mu)/n))
    }
    wi[ii] <- soma1 / (1+soma2)
  }
  
  Vi <- sqrt((1-(k*(wi^2)))*(1-(wi^2)/k))
  A0i <- 1 / wi^2
  sqrA0i <- 1 / wi
  B0i <- ((sig0*Vi)^2 + (w*wi)^2) / ((1+sig0^2*wi^2)^2)
  ## not used: ## B1i <- (2 * sig0 * Vi) / (1 + sig0^2 * wi^2)

  ##Gain T0:
  if (n %% 2) # odd
    T0 <- sig0 * prod(B0i / A0i) * sqrt(ws)
  else
    T0 <- 10^(-0.05*Rp) * prod(B0i / A0i)
  
  ##zeros:
  zer <- c(1i*sqrA0i, -1i*sqrA0i)

  ##poles:
  pol <- c((-2*sig0*Vi + 2*1i*wi*w) / (2*(1 + sig0^2*wi^2)), (-2*sig0*Vi - 2*1i*wi*w) / (2*(1 + sig0^2*wi^2)))
  
  ##If n odd, there is a real pole  -sig0:
  if (n %% 2) # odd
    pol <- c(pol, -sig0)

  ##
  pole <- sqrt(ws) * pol
  zero <- sqrt(ws) * zer

  Zpg(zero = zero, pole = pole, gain = T0)
} 
