## Copyright (C) 2004 Julius O. Smith III
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  The GNU
## General Public License has more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, write to the Free
## Software Foundation, 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
##
## Plot, demos, and help info copied and adapted from 
## grpdelay.m, Copyright (C) 2000 Paul Kienzle

## Compute the group delay of a filter.
##
## [g, w] = grpdelay(b)
##   returns the group delay g of the FIR filter with coefficients b.
##   The response is evaluated at 512 angular frequencies between 0 and
##   pi. w is a vector containing the 512 frequencies.
##   The group delay is in units of samples.  It can be converted
##   to seconds by multiplying by the sampling period (or dividing by
##   the sampling rate fs).
##
## [g, w] = grpdelay(b,a)
##   returns the group delay of the rational IIR filter whose numerator
##   has coefficients b and denominator coefficients a.
##
## [g, w] = grpdelay(b,a,n)
##   returns the group delay evaluated at n angular frequencies.  For fastest
##   computation n should factor into a small number of small primes.
##
## [g, w] = grpdelay(b,a,n,'whole')
##   evaluates the group delay at n frequencies between 0 and 2*pi.
##
## [g, f] = grpdelay(b,a,n,Fs)
##   evaluates the group delay at n frequencies between 0 and Fs/2.
##
## [g, f] = grpdelay(b,a,n,'whole',Fs)
##   evaluates the group delay at n frequencies between 0 and Fs.
##
## [g, w] = grpdelay(b,a,w)
##   evaluates the group delay at frequencies w (radians per sample).
##
## [g, f] = grpdelay(b,a,f,Fs)
##   evaluates the group delay at frequencies f (in Hz).
##
## grpdelay(...)
##   plots the group delay vs. frequency.
##
## If the denominator of the computation becomes too small, the group delay
## is set to zero.  (The group delay approaches infinity when
## there are poles or zeros very close to the unit circle in the z plane.)
##
## Theory: group delay, g(w) = -d/dw [arg{H(e^jw)}],  is the rate of change of
## phase with respect to frequency.  It can be computed as:
##
##               d/dw H(e^-jw)
##        g(w) = -------------
##                 H(e^-jw)
##
## where
##         H(z) = B(z)/A(z) = sum(b_k z^k)/sum(a_k z^k).
##
## By the quotient rule,
##                    A(z) d/dw B(z) - B(z) d/dw A(z)
##        d/dw H(z) = -------------------------------
##                               A(z) A(z)
## Substituting into the expression above yields:
##                A dB - B dA 
##        g(w) =  ----------- = dB/B - dA/A
##                    A B
##
## Note that,
##        d/dw B(e^-jw) = sum(k b_k e^-jwk)
##        d/dw A(e^-jw) = sum(k a_k e^-jwk)
## which is just the FFT of the coefficients multiplied by a ramp.
##
## As a further optimization when nfft>>length(a), the IIR filter (b,a)
## is converted to the FIR filter conv(b,fliplr(conj(a))).
## For further details, see 
## http://ccrma.stanford.edu/~jos/filters/Numerical_Computation_Group_Delay.html


grpdelay <- function(filt, ...) UseMethod("grpdelay")

print.grpdelay <- function(x, ...){
    cat("- Group delay (gd) calculated at", x$ns, "points.\n")
    cat("- Frequencies (w) given in", if(x$HzFlag) "*Hz*." else "*radians*.", "\n")
    temp <- data.frame(do.call("cbind", x[c("gd", "w")]))
    if(nrow(temp) > 8L){
        print(head(temp, n=4L), row.names=FALSE, ...)
        cat(" .......  .......\n")
    } else print(temp, row.names=FALSE, ...)
    invisible(x)
}

plot.grpdelay <- function(x, xlab = if(x$HzFlag) 'Hz' else 'radian/sample', 
    ylab = 'Group delay (samples)', type = "l", ...){
  plot(x$w[1:x$ns], x$gd[1:x$ns], 
    xlab = xlab, ylab = ylab, type = type, ...)
}

grpdelay.Arma <- function(filt, ...) # IIR
  grpdelay(filt$b, filt$a, ...)

grpdelay.Ma <- function(filt, ...) # FIR
  grpdelay(filt$b, 1, ...)

grpdelay.Zpg <- function(filt, ...) # Zero-pole-gain ARMA
  grpdelay(as.Arma(filt), ...)

grpdelay.default <- function(filt, a = 1, n = 512, whole = FALSE, Fs = NULL, ...)   {
  nfft <- n
  b <- filt
  if (whole == "whole" || whole) {
    whole <- TRUE
  } else {
    whole <- FALSE
    nfft <- 2*nfft
  }
  if (is.null(Fs)) {
    HzFlag <- FALSE
    Fs <- 1
  } else {
    HzFlag <- TRUE
  }    
  
  w <- Fs * (0:(nfft-1)) / nfft
  if (!HzFlag)
    w <- w * 2 * pi

  oa <- length(a)-1             # order of a(z)
  if (oa < 0) {
    a <- 1
    oa <- 0
  }
  ob <- length(b)-1             # order of b(z)
  if (ob < 0) {
    b <- 1
    ob <- 0
  }                         
  oc <- oa + ob                 # order of c(z)
  
  c <- conv(b, rev(Conj(a)))       # c(z) <- b(z)*conj(a)(1/z)*z^(-oa)
#  c <- conv(b, rev(a))       # c(z) <- b(z)*conj(a)(1/z)*z^(-oa)
  cr <- c * (0:oc)                 # cr(z) <- derivative of c wrt 1/z 
  num <- fft(postpad(cr, nfft))
  den <- fft(postpad(c, nfft))
#  minmag <- 10*eps
#  polebins <- which(abs(den) < minmag)
  polebins <- which(abs(den) < 2 * .Machine$double.eps)
  if(any(polebins)){
    warning('grpdelay: setting group delay to 0 at singularity')
    num[polebins] <- 0
    den[polebins] <- 1
  }
    
  gd <- Re(num / den) - oa

  if (!whole) {
    ns <- nfft/2     # Matlab convention ... should be nfft/2 + 1
    gd <- gd[1:ns]
    w <- w[1:ns]
  } else {
    ns <- nfft # used in plot below
  } 
  res <- list(gd = gd, w = w, ns = ns, HzFlag = HzFlag)
  class(res) = "grpdelay"
  res
} 

# ------------------------ DEMOS -----------------------

###grpdelay(c(1,0.9),1,512,'whole',1)

###b = poly(c(1/0.9*exp(1i*pi*0.2), 0.9*exp(1i*pi*0.6)))
###a = poly(c(0.9*exp(-1i*pi*0.6), 1/0.9*exp(-1i*pi*0.2)))
####! title ('Two Zeros and Two Poles')
###grpdelay(b,a,512,'whole',1)

#!demo % 1
#! %--------------------------------------------------------------
#! % From Oppenheim and Schafer, a single zero of radius r=0.9 at
#! % angle pi should have a group delay of about -9 at 1 and 1/2
#! % at zero and 2*pi.
#! %--------------------------------------------------------------
#! title ('Zero at z = -0.9')
###grpdelay(c(1,0.9),1,512,'whole',1)
#! hold on
#! xlabel('Normalized Frequency (cycles/sample)')
#! stem([0, 0.5, 1],[0.5, -9, 0.5],'*b;target;')
#! hold off
#! 
#!demo % 2
#! %--------------------------------------------------------------
#! % confirm the group delays approximately meet the targets
#! % don't worry that it is not exact, as I have not entered
#! % the exact targets.
#! %--------------------------------------------------------------
###b = poly(c(1/0.9*exp(1i*pi*0.2), 0.9*exp(1i*pi*0.6)))
###a = poly(c(0.9*exp(-1i*pi*0.6), 1/0.9*exp(-1i*pi*0.2)))
####! title ('Two Zeros and Two Poles')
###grpdelay(b,a,512,'whole',1)
#! hold on
#! xlabel('Normalized Frequency (cycles/sample)')
#! stem([0.1, 0.3, 0.7, 0.9], [9, -9, 9, -9],'*b;target;')
#! hold off

#!demo % 3
#! %--------------------------------------------------------------
#! % fir lowpass order 40 with cutoff at w=0.3 and details of
#! % the transition band [.3, .5]
#! %--------------------------------------------------------------
#! subplot(211)
#! Fs = 8000;     % sampling rate
#! Fc = 0.3*Fs/2; % lowpass cut-off frequency
#! nb = 40
#! b = fir1(nb,2*Fc/Fs); % matlab freq normalization: 1=Fs/2 
#! [H,f] = freqz(b,1,[],1)
#! [gd,f] = grpdelay(b,1,[],1)
#! title(sprintf('b = fir1(%d,2*%d/%d);',nb,Fc,Fs))
#! xlabel('Normalized Frequency (cycles/sample)')
#! ylabel('Amplitude Response (dB)')
#! grid('on')
#! plot(f,20*log10(abs(H)))
#! subplot(212)
#! del = nb/2; % should equal this
#! title(sprintf('Group Delay in Pass-Band (Expect %d samples)',del))
#! ylabel('Group Delay (samples)')
#! axis([0, 0.2, del-1, del+1])
#! plot(f,gd)
#! axis(); oneplot()

#!demo % 4
#! %--------------------------------------------------------------
#! % IIR bandstop filter has delays at [1000, 3000]
#! %--------------------------------------------------------------
#! Fs = 8000
#! [b, a] = cheby1(3, 3, 2*[1000, 3000]/Fs, 'stop')
#! [H,f] = freqz(b,a,[],Fs)
#! [gd,f] = grpdelay(b,a,[],Fs)
#! subplot(211)
#! title('[b,a] = cheby1(3, 3, 2*[1000, 3000]/Fs, \'stop\');')
#! xlabel('Frequency (Hz)')
#! ylabel('Amplitude Response')
#! grid('on')
#! plot(f,abs(H))
#! subplot(212)
#! title('[gd,f] = grpdelay(b,a,[],Fs);')
#! ylabel('Group Delay (samples)')
#! plot(f,gd)
#! oneplot()

# ------------------------ TESTS -----------------------

# TEST 0A:a
#! [gd,w] = grpdelay([0,1])
#! [gd,w] = grpdelay([0,1],1)

#!test % 0A
#! [gd,w] = grpdelay([0,1],1,4)
#! assert(gd,[1;1;1;1])
#! assert(w,pi/4*[0:3]',10*eps)

#!test % 0B
#! [gd,w] = grpdelay([0,1],1,4,'whole')
#! assert(gd,[1;1;1;1])
#! assert(w,pi/2*[0:3]',10*eps)

#!test % 0C
#! [gd,f] = grpdelay([0,1],1,4,0.5)
#! assert(gd,[1;1;1;1])
#! assert(f,1/16*[0:3]',10*eps)

#!test % 0D
#! [gd,w] = grpdelay([0,1],1,4,'whole',1)
#! assert(gd,[1;1;1;1])
#! assert(w,1/4*[0:3]',10*eps)

#!test % 0E
#! [gd,f] = grpdelay([1, -0.9i],[],4,'whole',1)
#! gd0 = 0.447513812154696; gdm1 =0.473684210526316
#! assert(gd,[gd0;-9;gd0;gdm1],20*eps)
#! assert(f,1/4*[0:3]',10*eps)

#!test % 1:
#! gd= grpdelay(1,[1,.9],4)
#! assert(gd, [-0.47368;-0.46918;-0.44751;-0.32316],1e-5)

#!test % 2:
#! gd = grpdelay([1,2],[1,0.5,.9],4)
#! assert(gd,[-0.29167;-0.24218;0.53077;0.40658],1e-5)

#!test % 3
#! b1=[1,2];a1f=[0.25,0.5,1];a1=fliplr(a1f)
#! % gd1=grpdelay(b1,a1,4)
#! gd=grpdelay(conv(b1,a1f),1,4)
#! assert(gd, [0.095238;0.239175;0.953846;1.759360],1e-5)

#!test
#! Fs = 8000
#! [b, a] = cheby1(3, 3, 2*[1000, 3000]/Fs, 'stop')
#! [h, w] = grpdelay(b, a, 256, 'half', Fs)
#! [h2, w2] = grpdelay(b, a, 512, 'whole', Fs)
#! assert (size(h), size(w))
#! assert (length(h), 256)
#! assert (size(h2), size(w2))
#! assert (length(h2), 512)
#! assert (h, h2(1:256))
#! assert (w, w2(1:256))
