#' R Wireless Communications Package
#'
#' A communications simulation package supporting common modulations formats
#' including BPSK, QPSK, 8-PSK, 16-PSK, 16-QAM and 64-QAM. The package includes
#' AWGN noise generation, and raised cosine and square root raised cosine
#' pulse shaping. It also includes functions for plotting
#' constellation diagrams, density plots, stem plots and eye diagrams. The rwirelesscom
#' package includes the followng functions:
#' \itemize{
#' \item fNo(),
#' \item fbpskmod(), fbpskdemod(),
#' \item f8pskmod(), f8pskdemod(),
#' \item f16pskmod(), f16pskdemod(),
#' \item f16qammod(), f16qamdemod(),
#' \item f64qammod(), f64qamdemod,
#' \item rcosine(), sqrtrcosine(), sinc(),
#' \item iqscatterplot(), iqdensityplot(),
#' \item stemplot(), eyediagram()
#' }
#' Together these functions
#' enable the evaluation of bit error and symbol error rates, evalutation of pulse shaping and inter-symbol
#' interferance and support visualization of the respective signals and noise.
#'
#' @docType package
#' @name rwirelesscom
#' @import
#'   ggplot2
#'   stats
#'   Rcpp
NULL

#' AWGN
#'
#' Generates a vector of normally distributed noise samples with mean of zero and noise spectral density (No/2), a.k.a. AWGN.
#' @param N - number of noise samples
#' @param No - noise spectral density
#' @param type - "real" or "complex" defaults to real
#' @family rwirelesscom functions
#' @return returns a vector of distributed noise samples of length N, mean of zero and variance of No/2
#' @examples
#' n <- fNo(N=10,No=10)
#' bits <- sample(0:1,10, replace=TRUE)
#' s=fbpskmod(bits)
#' r=s+n
#'
#' n <- fNo(N=20,No=10,type="complex")
#' bits <- sample(0:1,20, replace=TRUE)
#' s=fqpskmod(bits)
#' r=s+n
#' @export
fNo <- function(N,No,type="real") {
  if (type=="real") n = sqrt(No/2)*rnorm(N, 0, 1)
  else if(type=="complex")  n=sqrt(No/2)*complex(real=rnorm(N,0,1),imaginary=rnorm(N,0,1))
  else n=c(rep(0,N))
  return(n)
}

#' BPSK Modulator
#'
#' Receives a vector of bits, each with value 0 or 1, and outputs a
#' vector with values 1 and -1, respectively.
#' @param bits - vector of bits (0's and 1's)
#' @param Ns - N samples per symbol (default, Ns = 1)
#' @param p - a vector defining the pulse shape of the transmitted waveform (default, p = 1)
#' @family rwirelesscom functions
#' @return Returns a BPSK modulated vector, each element taking on values of 1 or -1. If Ns > 1
#' then the returned signal is shaped with pulse shape, p.
#' @examples
#' bits <- sample(0:1,10, replace=TRUE)
#' fbpskmod(bits)
#' @export
fbpskmod <- function(bits,Ns=1,p=1) {
  s <- sapply(bits, function(x) if (x==0) s=-1 else s=x)
  if (Ns >= 1) {
    Nsymbols=length(s)
    x <-as.vector(rbind(s,matrix(c(0),Ns-1,Nsymbols)))
    s= convolve(x,p,type="open")
  }
  return(s)
}

#' BPSK Demodulator
#'
#' Receives a vector of real values, corresponding to a
#' BPSK modulated signal transmitted through a communications channel
#' (e.g., signal plus noise). An input value < 1 is mapped to an
#' output value of 0, otherwise to a value of 1.
#' @param r - received signal vector
#' @return returns a vector of 1's and 0's corresponding to BPSK demodulation of the input vector
#' @examples
#' Eb=1
#' Nbits=10
#' bits <- sample(0:1,Nbits, replace=TRUE)
#' s <- fbpskmod(bits)
#' EbNodB=8
#' No = Eb/(10^(EbNodB/10))
#' n <- fNo(Nbits,No)
#' r <- s+n
#' bitsr <- fbpskdemod(r)
#' biterrs<-bits[bitsr!=bits]
#' Pberr=length(biterrs)/length(bits)
#' @family rwirelesscom functions
#' @export
fbpskdemod <- function(r) {
  r <- sapply(r,function(x) if (x>=0) r=1 else r=0)
  return(r)
}

#' QPSK Modulator
#'
#' Receives a vector of bits (1's and 0's). The 1's and 0's are
#' mapped to in-phase (real) and quadrature (imaginary) components.
#' Correspondingly, a bit of 1 is mapped to +1/sqrt(2), otherwise to -1/sqrt(2)
#' according to the following mapping.
#'  \tabular{cc}{
#' input \tab output \cr
#' 00 \tab  (-1 - 1i) / sqrt(2)  \cr
#' 01 \tab  (-1 + 1i) / sqrt(2) \cr
#' 10 \tab  (+1 - 1i) / sqrt(2) \cr
#' 11 \tab  (+1 + 1i) / sqrt(2)
#' }
#' @param bits - received vector of bits (0's and 1's).
#' @param Ns - N samples per symbol (default, Ns = 1)
#' @param p - a vector defining the pulse shape of the transmitted waveform (default, p = 1)
#' @return Returns a complex vector of QPSK symbols. If Ns > 1
#' then the returned signal is shaped with pulse shape, p.
#' @examples
#' M=4
#' Nsymbols=10
#' Nbits=log2(M)*Nsymbols
#' bits <- sample(0:1,Nbits, replace=TRUE)
#' s <- fqpskmod(bits)
#' @family rwireless
#' @export
fqpskmod <- function(bits,Ns=1,p=1) {
  bi <- bits[seq(1,length(bits),2)]
  bq <- bits[seq(2,length(bits),2)]
  si <- (1/sqrt(2))*sapply(bi, function(x) if (x==0) r=-1 else r=x)
  sq <- (1/sqrt(2))*sapply(bq, function(x) if (x==0) r=-1 else r=x)
  s <- complex(real=si, imaginary=sq)
  if (Ns >= 1) {
    Nsymbols=length(s)
    x <-as.vector(rbind(s,matrix(c(0),Ns-1,Nsymbols)))
    s= convolve(x,p,type="open")
  }
  return(s)
}

#' QPSK Demodulator
#'
#' Receives a vector of complex values, r, corresponding to a
#' QPSK modulated signal transmitted through a communications channel
#' (e.g., signal plus noise). The received signal, r, is mapped to its in-phase (real part) and
#' quadrature parts (imaginary part) and demodulated, such that two binary bits are output
#' for each value of r. If the in-phase part is > 0 then the corresponding output bit value is 1,
#' otherwise 0. Similarly, if the quadrature part (imaginary) > 0 then the corresponding bit value
#' is 1, otherwise 0.
#' @param r - received signal plus noise.
#' @family rwirelesscom functions
#' @return returns a vector of 1's and 0's, 2 bits per input element (i.e., QPSK symbol)
#' @examples
#' M=4
#' Es=1
#' Nsymbols=10
#' Nbits=log2(M)*Nsymbols
#' Eb=Es/log2(M)
#' bits <- sample(0:1,Nbits, replace=TRUE)
#' s <- fqpskmod(bits)
#' EbNodB=8
#' No = Eb/(10^(EbNodB/10))
#' n <- fNo(Nsymbols,No,type="complex")
#' r <- s+n
#' bitsr <- fqpskdemod(r)
#' biterrs<-bits[bitsr!=bits]
#' Pberr=length(biterrs)/length(bits)
#' @export
fqpskdemod <- function(r) {
  ri <- Re(r)
  rq <- Im(r)
  r2 <- matrix(rbind(ri,rq),length(ri)+length(rq),byrow=TRUE)
  r <- sapply(r2,function(x) if (x>=0) r=1 else r=0)
  return(r)
}

#' 8-PSK Modulator
#'
#' Receives a vector of bits (1's and 0's). The received vector is mapped
#' to in-phase (real) and quadrature (imaginary) components, according to
#' a Binary Reflective Gray Code (BRGC, see reference). Each received pair of bits
#' are are mapped to 8-PSK symbols,
#' with \eqn{sqrt(E_{s})}{sqrt(Es)} (symbol energy) = 1.The bit to symbol mapping is illustrated in the following constellation diagram.
#' \tabular{cc}{
#' input \tab output \cr
#' 000 \tab \eqn{  0 }  \cr
#' 001 \tab \eqn{  \pi/4 } \cr
#' 011 \tab \eqn{  \pi/2} \cr
#' 010 \tab \eqn{  3 \pi/4} \cr
#' 110 \tab  \eqn{ pi} \cr
#' 111 \tab  \eqn{ -3 \pi/4} \cr
#' 101 \tab  \eqn{ - \pi/2} \cr
#' 100 \tab  \eqn{ - \pi/4} \cr
#' }
#' Reference: E. Agrell, J Lassing, E. Strom, and T. Ottosson, Gray Coding for Multilevel Constellations In Gaussian Noise, IEEE Transactions on Communications, Vol. 53, No. 1, January 2007
#' @param bits - vector of bits (0's and 1's).
#' @param Ns - N samples per symbol (default, Ns = 1)
#' @param p - a vector defining the pulse shape of the transmitted waveform (default, p = 1)
#' @return Returns a complex vector of length = (length(bits) mod 3), 8-PSK symbols. If Ns > 1
#' then the returned signal is shaped with pulse shape, p.
#' @examples
#' M=8
#' Nsymbols=10
#' Nbits=log2(M)*Nsymbols
#' bits <- sample(0:1,Nbits, replace=TRUE)
#' s <- f8pskmod(bits)
#' @family rwirelesscom functions
#' @export
f8pskmod <- function(bits,Ns=1,p=1) {
  # receive symbolbits matrix (Nsym rows x Log2(M) cols )
  # transform to symbolcodes
  # transform to xreal, yimag
  M=8
  Nsymbols = length(bits) %/% log2(M)
  symbolbits<-matrix(bits,Nsymbols,log2(M),byrow=TRUE)
  m8=cbind(c(0,0,1),c(0,2,0),c(4,0,0))
  symbolbitsm8 = symbolbits %*% m8
  symbolcodes <- apply(symbolbitsm8,1,sum)
  s<-sapply(symbolcodes,ft8pskbitmap)
  if (Ns >= 1) {
    Nsymbols=length(s)
    x <-as.vector(rbind(s,matrix(c(0),Ns-1,Nsymbols)))
    s= convolve(x,p,type="open")
  }
  return(s)
}

ft8pskbitmap <- function(x) {                                        #b2b1b0
  if (x == 0) r = complex(real=1,imaginary=0)                         # 000
  else if (x == 1) r = complex(real=1/sqrt(2),imaginary= 1/sqrt(2))   # 001
  else if (x == 3) r = complex(real=0, imaginary=1)                   # 011
  else if (x == 2) r = complex(real=-1/sqrt(2),imaginary=1/sqrt(2))   # 010
  else if (x == 6) r = complex(real=-1,imaginary=0)                   # 110
  else if (x == 7) r = complex(real=-1/sqrt(2), imaginary=-1/sqrt(2)) # 111
  else if (x == 5) r = complex(real=0, imaginary=-1)                  # 101
  else if (x == 4) r = complex(real=1/sqrt(2),imaginary=-1/sqrt(2))   # 100
  else {print(x)
    r=100}
  return(r)
}

#' 8-PSK Demodulator
#'
#' Receives a vector of complex values, r, corresponding to an
#' 8-PSK modulated signal transmitted through a communications channel
#' (e.g., signal plus noise). Three bits are output for each received symbol
#' according to the following decision rules
#' \tabular{cc}{
#' input \tab output \cr
#'  \eqn{ -\pi/8 \ge Arg(r) < \pi/8} \tab 000  \cr
#'  \eqn{  \pi/8 \ge Arg(r) < 3 \pi/8} \tab 001 \cr
#'  \eqn{ 3 \pi/8 \ge Arg(r) < 5 \pi/8} \tab 011 \cr
#'  \eqn{ 5 \pi/8 \ge Arg(r) < 7 \pi/8} \tab 010 \cr
#'  \eqn{ 7 \pi/8 \ge Arg(r) < 9 \pi/8} \tab 110 \cr
#'  \eqn{ -7 \pi/8 \ge Arg(r) < -5 \pi/8} \tab 111 \cr
#'  \eqn{ -5 \pi/8 \ge Arg(r) < -3 \pi/8} \tab 101 \cr
#'  \eqn{ -3 \pi/8 \ge Arg(r) < - \pi/8}\tab  100
#' }
#' @param r - received signal
#' @return returns a vector of 1's and 0's, 3 bits per input element (i.e., 8-PSK symbol)
#' @examples
#' M=8
#' Es=1
#' Eb = Es/log2(M)
#' Nsymbols=10
#' Nbits=log2(M)*Nsymbols
#' bits <- sample(0:1,Nbits, replace=TRUE)
#' s <- f8pskmod(bits)
#' EbNodB=7
#' No = Eb/(10^(EbNodB/10))
#' n <- fNo(Nsymbols,No,type="complex")
#' r <- s+n
#' bitsr <- f8pskdemod(r)
#' biterrs<-bits[bitsr!=bits]
#' b<-factor(bits)
#' Pberr=length(biterrs)/length(bits)
#' @family rwirelesscom functions
#' @export
#'
f8pskdemod <- function(r) {
  r2<-sapply(r,fr8pskbitmap)
  bits2 <- as.vector(matrix(r2,1,length(r2)))
  return(bits2)
}

fr8pskbitmap <- function(r) {
  argr <- Arg(r)
  if ( argr >= -pi/8 && argr < pi/8 ) symbolbits = c(0,0,0)
  else if ( argr >= pi/8 && argr < 3*pi/8 ) symbolbits = c(0,0,1)
  else if ( argr >= 3*pi/8 && argr < 5*pi/8) symbolbits = c(0,1,1)
  else if ( argr >= 5*pi/8 && argr < 7*pi/8) symbolbits = c(0,1,0)
  else if ( (argr >= 7*pi/8 && argr < 9*pi/8 ) || (argr > -9*pi/8 && argr < -7*pi/8 )) symbolbits = c(1,1,0)
  else if ( argr >= -7*pi/8 && argr  < -5*pi/8) symbolbits = c(1,1,1)
  else if ( argr >= -5*pi/8 && argr  < -3*pi/8) symbolbits = c(1,0,1)
  else if ( argr >= -3*pi/8 && argr  < -pi/8) symbolbits = c(1,0,0)
  else symbolbits=c(-1,-1,-1)
  return(symbolbits)
}

#' 16-PSK Modulator
#'
#' Receives a vector of bits (1's and 0's). The received vector is mapped
#' to in-phase (real) and quadrature (imaginary) components, according to
#' a Binary Reflective Gray Code (BRGC, see reference). Each received pair of bits
#' are are mapped to 16-PSK symbols,
#' with Es (symbol energy) = 1.The bit to symbol mapping is illustrated in the following table.
#' \tabular{cc}{
#' input \tab output \cr
#' 0000 \tab \eqn{    0 }  \cr
#' 0001 \tab \eqn{   \pi/8 } \cr
#' 0011 \tab \eqn{   \pi/4} \cr
#' 0010 \tab \eqn{   3\pi/8} \cr
#' 0110 \tab  \eqn{  \pi/2} \cr
#' 0111 \tab  \eqn{ 5 \pi/8} \cr
#' 0101 \tab  \eqn{ 3 \pi/4} \cr
#' 0100 \tab  \eqn{ 7 \pi/8} \cr
#' 1100 \tab \eqn{  -1 }  \cr
#' 1101 \tab \eqn{  -7 \pi/8 } \cr
#' 1111 \tab \eqn{  -3\pi/4} \cr
#' 1110 \tab \eqn{  -5 \pi/8} \cr
#' 1010 \tab  \eqn{ - \pi/2} \cr
#' 1011 \tab  \eqn{ -3 \pi/8} \cr
#' 1001 \tab  \eqn{ - \pi/4} \cr
#' 1000 \tab  \eqn{ - \pi/8}
#' }
#' Reference: E. Agrell, J Lassing, E. Strom, and T. Ottosson, Gray Coding for Multilevel Constellations In Gaussian Noise, IEEE Transactions on Communications, Vol. 53, No. 1, January 2007
#' @param bits - vector of bits (0's and 1's).
#' @param Ns - N samples per symbol (default, Ns = 1)
#' @param p - a vector defining the pulse shape of the transmitted waveform (default, p = 1)
#' @return Returns a complex vector of length of 8-PSK symbols. If Ns > 1
#' then the returned signal is shaped with pulse shape, p.
#' @examples
#' M=16
#' Nsymbols=20
#' Nbits=log2(M)*Nsymbols
#' bits <- sample(0:1,Nbits, replace=TRUE)
#' s <- f16pskmod(bits)
#' @family rwirelesscom functions
#' @export
f16pskmod <- function(bits,Ns=1,p=1) {
  # receive symbolbits matrix (Nsym rows x Log2(M) cols )
  # transform to symbolcodes
  # transform to xreal, yimag
  M=16
  Nsymbols = length(bits) %/% log2(M)
  symbolbits<-matrix(bits,Nsymbols,log2(M),byrow=TRUE)
  m16=cbind(c(0,0,0,1),c(0,0,2,0),c(0,4,0,0),c(8,0,0,0))
  symbolbitsm16 = symbolbits %*% m16
  symbolcodes <- apply(symbolbitsm16,1,sum)
  s<-sapply(symbolcodes,ft16pskbitmap)
  if (Ns >= 1) {
    Nsymbols=length(s)
    x <-as.vector(rbind(s,matrix(c(0),Ns-1,Nsymbols)))
    s= convolve(x,p,type="open")
  }
  return(s)
}

ft16pskbitmap <- function(x) {                                           #b3b2b1b0
  if (x == 0) r = complex(real=1,imaginary=0)                             # 0000
  else if (x == 1) r = complex(real=cos(pi/8),imaginary=sin(pi/8))        # 0001
  else if (x == 3) r = complex(real=cos(pi/4), imaginary=sin(pi/4))       # 0011
  else if (x == 2) r = complex(real=cos(3*pi/8),imaginary=sin(3*pi/8))    # 0010
  else if (x == 6) r = complex(real=0,imaginary=1)                        # 0110
  else if (x == 7) r = complex(real=cos(5*pi/8), imaginary=sin(5*pi/8))   # 0111
  else if (x == 5) r = complex(real=cos(3*pi/4), imaginary=sin(3*pi/4))   # 0101
  else if (x == 4) r = complex(real=cos(7*pi/8),imaginary=sin(7*pi/8))    # 0100
  else if (x == 12) r = complex(real=-1,imaginary=0)                      # 1100
  else if (x == 13) r = complex(real=cos(-7*pi/8),imaginary=sin(-7*pi/8)) # 1101
  else if (x == 15) r = complex(real=cos(-3*pi/4),imaginary=sin(-3*pi/4)) # 1111
  else if (x == 14) r = complex(real=cos(-5*pi/8),imaginary=sin(-5*pi/8)) # 1110
  else if (x == 10) r = complex(real=0,imaginary=-1)                      # 1010
  else if (x == 11) r = complex(real=cos(-3*pi/8),imaginary=sin(-3*pi/8)) # 1011
  else if (x == 9) r = complex(real=cos(- pi/4),imaginary=sin(- pi/4))    # 1001
  else if (x == 8) r = complex(real=cos(-pi/8),imaginary=sin(-pi/8))      # 1000
  else {print(x)
    r=100}
  return(r)
}

#' 16-PSK Demodulator
#'
#' Receives a vector of complex values, r, corresponding to an
#' 8-PSK modulated signal transmitted through a communications channel
#' (e.g., signal plus noise). Three bits are output for each received symbol
#' according to the following decision rules
#' \tabular{cc}{
#' input \tab output \cr
#'  \eqn{ -\pi/16 \ge Arg(r) < \pi/16}        \tab 0000  \cr
#'  \eqn{  \pi/16 \ge Arg(r) < 3 \pi/16}     \tab 0001 \cr
#'  \eqn{ 3 \pi/16 \ge Arg(r) < 5 \pi/16}     \tab 0011 \cr
#'  \eqn{ 5 \pi/16 \ge Arg(r) < 7 \pi/16}     \tab 0010 \cr
#'  \eqn{ 7 \pi/16 \ge Arg(r) < 9 \pi/16}    \tab 0110 \cr
#'  \eqn{ 9 \pi/16 \ge Arg(r) < 11 \pi/16}   \tab 0111 \cr
#'  \eqn{ 11 \pi/16 \ge Arg(r) < 13 \pi/16}   \tab 0101 \cr
#'  \eqn{ 13 \pi/16 \ge Arg(r) < 15 \pi/16}   \tab 0100 \cr
#'  \eqn{ 15 \pi/16 \ge Arg(r) < 17 \pi/16}   \tab 1100  \cr
#'  \eqn{ -15 \pi/16 \ge Arg(r) < -13 \pi/16}  \tab 1101 \cr
#'  \eqn{ -13 \pi/16 \ge Arg(r) < -11 \pi/16} \tab 1111 \cr
#'  \eqn{ -11 \pi/16 \ge Arg(r) < -9 \pi/16}   \tab 1110 \cr
#'  \eqn{ -9 \pi/16 \ge Arg(r) < -7 \pi/16}    \tab 1010 \cr
#'  \eqn{ -7 \pi/16 \ge Arg(r) < -5 \pi/16}    \tab 1011 \cr
#'  \eqn{ -5 \pi/16 \ge Arg(r) < -3 \pi/16}    \tab 1001 \cr
#'  \eqn{ -3 \pi/16 \ge Arg(r) < - \pi/16}     \tab 1000
#' }
#' @param r - received signal
#' @return returns a vector of 1's and 0's, 3 bits per input element (i.e., 8-PSK symbol)
#' @examples
#' M=16
#' Es=1
#' Eb = Es/log2(M)
#' Nsymbols=20
#' Nbits=log2(M)*Nsymbols
#' bits <- sample(0:1,Nbits, replace=TRUE)
#' s <- f16pskmod(bits)
#' EbNodB=7
#' No = Eb/(10^(EbNodB/10))
#' n <- fNo(Nsymbols,No,type="complex")
#' r <- s+n
#' bitsr <- f16pskdemod(r)
#' biterrs<-bits[bitsr!=bits]
#' b<-factor(bits)
#' Pberr=length(biterrs)/length(bits)
#' @family rwirelesscom functions
#' @export
#'
f16pskdemod <- function(r) {
  r2<-sapply(r,fr16pskbitmap)
  bits2 <- as.vector(matrix(r2,1,length(r2)))
  return(bits2)
}

fr16pskbitmap <- function(r) {
  argr <- Arg(r)
  if ( argr >= -pi/16 && argr < pi/16 ) symbolbits = c(0,0,0,0)
  else if ( argr >= pi/16 && argr < 3*pi/16 ) symbolbits = c(0,0,0,1)
  else if ( argr >= 3*pi/16 && argr < 5*pi/16) symbolbits = c(0,0,1,1)
  else if ( argr >= 5*pi/16 && argr < 7*pi/16) symbolbits = c(0,0,1,0)
  else if ( argr >= 7*pi/16 && argr < 9*pi/16)  symbolbits= c(0,1,1,0)
  else if ( argr >= 9*pi/16 && argr  < 11*pi/16) symbolbits = c(0,1,1,1)
  else if ( argr >= 11*pi/16 && argr  < 13*pi/16) symbolbits = c(0,1,0,1)
  else if ( argr >= 13*pi/16 && argr  < 15*pi/16) symbolbits = c(0,1,0,0)

  else if ( (argr >= 15*pi/16 && argr < 17*pi/16 ) || (argr > -17*pi/16 && argr < -15*pi/16 )) symbolbits = c(1,1,0,0)
  else if ( argr >= -15*pi/16 && argr < -13*pi/16 ) symbolbits = c(1,1,0,1)
  else if ( argr >= -13*pi/16 && argr < -11*pi/16) symbolbits = c(1,1,1,1)
  else if ( argr >= -11*pi/16 && argr < -9*pi/16) symbolbits = c(1,1,1,0)
  else if ( argr >= -9*pi/16 && argr < -7*pi/16) symbolbits = c(1,0,1,0)
  else if ( argr >= -7*pi/16 && argr  < -5*pi/16) symbolbits = c(1,0,1,1)
  else if ( argr >= -5*pi/16 && argr  < -3*pi/16) symbolbits = c(1,0,0,1)
  else if ( argr >= -3*pi/16 && argr  < -pi/16) symbolbits = c(1,0,0,0)
  else  {
       symbolbits=c(-1,-1,-1,-1)
       print(r)
  }
  return(symbolbits)
}

#' 16-QAM Modulator
#'
#' Receives a vector of bits (1's and 0's). The received vector is mapped
#' to an in-phase (real) and quadrature (imaginary) 16-QAM (4 bit) symbol according to a
#' a Binary Reflective Gray Code (BRGC, see reference). Each symbol has an average
#' symbol energy Es = 10, where in-phase and quadrature constellation points
#' take on values -3, -1, +1, +3, respectively. The bit to symbol mapping is illustrated in the following constellation diagram.
#' \tabular{cccc}{
#'  -3+3i  \tab -1+3i \tab +1+3i  \tab +3+3i \cr
#'  (1000) \tab (1001) \tab (1011) \tab (1010) \cr
#'  \tab \tab \tab  \cr
#'  -3+1i  \tab -1+1i \tab +1+1i  \tab +3+1i \cr
#'  (1100) \tab (1101) \tab (1111) \tab (1110) \cr
#'  \tab \tab \tab  \cr
#'  -3-1i  \tab -1-1i \tab +1-1i  \tab +3-1i \cr
#'  (0100) \tab (0101) \tab (0111) \tab (0110) \cr
#'  \tab \tab \tab  \cr
#'  -3-3i  \tab -1-3i \tab +1-3i  \tab +3-3i \cr
#'  (0000) \tab (0001) \tab (0011) \tab (0010)
#' }
#' Reference: E. Agrell, J Lassing, E. Strom, and T. Ottosson, Gray Coding for Multilevel Constellations In Gaussian Noise, IEEE Transactions on Communications, Vol. 53, No. 1, January 2007
#' @param bits - received vector of bits (0's and 1's).
#' @param Ns - N samples per symbol (default, Ns = 1)
#' @param p - a vector defining the pulse shape of the transmitted waveform (default, p = 1)
#' @return Returns a complex vector of 16-QAM symbols. If Ns > 1
#' then the returned signal is shaped with pulse shape, p.
#' @examples
#' M=16
#' Nsymbols=100
#' Nbits=log2(M)*Nsymbols
#' bits <- sample(0:1,Nbits, replace=TRUE)
#' s <- f16qammod(bits)
#' @family rwirelesscom functions
#' @export
f16qammod <- function(bits,Ns=1,p=1) {
  # constellation quadrant x = 1+1j, 3+1j, 1+3j, 3+3j ...
  # sqrt(Es) = Sum((x * Conj(x)))/4 = 10
  # receive bits
  # transform to symbolbits (Nsym rows) length(bits)/4 x Log2(M) cols
  # transform to symbolcodes
  # transform to xreal(inphase), yimag(quadrature)
  M=16
  Nsymbols = length(bits) %/% log2(M)
  symbolbits<-matrix(bits,Nsymbols,log2(M),byrow=TRUE)
  m16=cbind(c(0,0,0,1),c(0,0,2,0),c(0,4,0,0),c(8,0,0,0))
  symbolbitsm16 = symbolbits %*% m16
  symbolcodes <- apply(symbolbitsm16,1,sum)
  s<-sapply(symbolcodes,ft16qambitmap)
  if (Ns >= 1) {
    Nsymbols=length(s)
    x <-as.vector(rbind(s,matrix(c(0),Ns-1,Nsymbols)))
    s= convolve(x,p,type="open")
  }
  return(s)
}
ft16qambitmap <- function(x) {
  # b1b0  R := 00 -> -3, 01 -> -1, 11-> +1, 10 -> +3             Im   Re
  # b3b2  I := 00 -> -3, 01 -> -1, 11-> +1, 10 -> +3             b3b2 b1b0
  if (x == 0)      r = complex(real=-3, imaginary=-3)            #00 00
  else if (x == 1) r = complex(real=-1, imaginary=-3)            #00 01
  else if (x == 3) r = complex(real=+1, imaginary=-3)            #00 11
  else if (x == 2) r = complex(real=+3, imaginary=-3)            #00 10
  else if (x == 4) r = complex(real=-3, imaginary=-1)            #01 00
  else if (x == 5) r = complex(real=-1, imaginary=-1)            #01 01
  else if (x == 7) r = complex(real=+1, imaginary=-1)            #01 11
  else if (x == 6) r = complex(real=+3, imaginary=-1)            #01 10
  else if (x == 12) r = complex(real=-3, imaginary=+1)           #11 00
  else if (x == 13) r = complex(real=-1, imaginary=+1)           #11 01
  else if (x == 15) r = complex(real=+1,imaginary=+1)            #11 11
  else if (x == 14) r = complex(real=+3,imaginary=+1)            #11 10
  else if (x == 8) r = complex(real=-3,imaginary=+3)             #10 00
  else if (x == 9) r = complex(real=-1,imaginary=+3)             #10 01
  else if (x == 11) r = complex(real=+1,imaginary=+3)            #10 11
  else if (x == 10) r = complex(real=+3,imaginary=+3)            #10 10
  else {print(x)
    r=100}
  return(r)
}

#' 16 QAM Demodulator
#'
#' Receives a vector of complex values, r, corresponding to a
#' 16-QAM modulated signal transmitted through a communications channel
#' (e.g., signal plus noise). Each received 16-QAM symbol (one symbol per r sample) is decoded into
#' 4 bits. The 16-QAM symbol decision regions are defined with respect to
#' the constellation generated by f16qammod(), where in-phase and quadrature
#' constellation points take on values -3, -1, +1, +3, respectively.
#' @param r - vector of complex vector
#' @return a vector of 1's and 0's, 4 bits per input element (16-QAM symbol)
#' @examples
#' M=16
#' Es=10
#' Eb = Es/log2(M)
#' Nsymbols=100
#' Nbits=log2(M)*Nsymbols
#' bits <- sample(0:1,Nbits, replace=TRUE)
#' s <- f16qammod(bits)
#' EbNodB=8
#' No = Eb/(10^(EbNodB/10))
#' n <- fNo(Nsymbols,No,type="complex")
#' r <- s+n
#' bitsr <- f16qamdemod(r)
#' biterrs<-bits[bitsr!=bits]
#' Pberr=length(biterrs)/length(bits)
#' @family rwirelesscom functions
#' @export
f16qamdemod <- function(r) {
  r2 <- sapply(r,fr16qambitmap)
  bitsr <- as.vector(matrix(r2,1,length(r2)))
  return(bitsr)
}
fr16qambitmap <- function(r) {
  if (Re(r) < -2 ) {
    if (Im(r) < -2) bits=c(0,0,0,0)
    else if (Im(r) < 0) bits=c(0,1,0,0)
    else if (Im(r) < 2) bits=c(1,1,0,0)
    else bits=c(1,0,0,0)
  }
  else if (Re(r) < 0 ) {
    if (Im(r) < -2) bits=c(0,0,0,1)
    else if (Im(r) < 0) bits=c(0,1,0,1)
    else if (Im(r) < 2) bits=c(1,1,0,1)
    else bits=c(1,0,0,1)
  }
  else if (Re(r) < 2) {
    if (Im(r) < -2) bits=c(0,0,1,1)
    else if (Im(r) < 0) bits=c(0,1,1,1)
    else if (Im(r) < 2) bits=c(1,1,1,1)
    else bits=c(1,0,1,1)
  }
  else if (Re(r) >=2 ) {
    if (Im(r) < -2) bits=c(0,0,1,0)
    else if (Im(r) < 0) bits=c(0,1,1,0)
    else if (Im(r) < 2) bits=c(1,1,1,0)
    else bits=c(1,0,1,0)
  }
  return(bits)
}

#' 64-QAM Modulator
#'
#' Receives a vector of bits (1's and 0's). The received vector is mapped to an
#' in-phase (real) and quadrature (imaginary) 64-QAM (6 bit) symbol according to
#' a Binary Reflective Gray Code (BRGC, see reference). In-phase and quadrature constellation points
#' take on values -7, -5, -3, -1, +1, +3, +5, +7, respectively, corresponding to a symbol energy Es = 42.
#' Mapping of bits to 64-QAM constellation points is illustrated in the following constellation.
#' \tabular{cccccccc}{
#' -7+7i  \tab -5+7i \tab +3+7i  \tab -1+7i \tab +1+7i  \tab +3+7i \tab +5+7i  \tab +7+7i \cr
#'  (100000) \tab (100001) \tab (100011) \tab (100010) \tab (100110) \tab (100111) \tab (100101) \tab  (100100) \cr
#'  \tab \tab \tab \tab \tab  \cr
#'   -7+5i  \tab -5+5i \tab +3+5i  \tab -1+5i \tab +1+5i  \tab +3+5i \tab +5+5i  \tab +7+5i \cr
#' (101000) \tab (101001) \tab (101011) \tab (101010) \tab (101110) \tab (101111) \tab (101101) \tab  (101100) \cr
#'   \tab \tab \tab \tab \tab  \cr
#'  -7+3i  \tab -5+3i \tab +3+3i  \tab -1+3i \tab +1+3i  \tab +3+3i \tab +5+3i  \tab +7+3i \cr
#' (111000) \tab (111001) \tab (111011) \tab (111010) \tab (111110) \tab (111111) \tab (111101) \tab (111100) \cr
#'   \tab \tab \tab \tab \tab   \cr
#'  -7+1i  \tab -5+1i \tab +3+1i  \tab -1+1i \tab +1+1i  \tab +3+1i \tab +5+1i  \tab +7+1i \cr
#'  (110000) \tab (110001) \tab (110011) \tab (110010) \tab (110110) \tab (110111) \tab (110101) \tab (110100) \cr
#'   \tab \tab \tab \tab \tab   \cr
#' -7-1i  \tab -5-1i \tab +3-1i  \tab -1-1i \tab +1-1i  \tab +3-1i \tab +5-1i  \tab +7-1i \cr
#'  (010000) \tab (010001) \tab (010011) \tab (010010) \tab (010110) \tab (010111) \tab (010101) \tab  (010100) \cr
#'  \tab \tab \tab \tab \tab  \cr
#'   -7-3i  \tab -5-3i \tab +3-3i  \tab -1-3i \tab +1-3i  \tab +3-3i \tab +5-3i  \tab +7-3i \cr
#' (011000) \tab (011001) \tab (011011) \tab (011010) \tab (011110) \tab (011111) \tab (011101) \tab  (011100) \cr
#'   \tab \tab \tab \tab \tab  \cr
#'  -7-5i  \tab -5-5i \tab +3-5i  \tab -1-5i \tab +1-5i  \tab +3-5i \tab +5-5i  \tab +7-5i \cr
#' (001000) \tab (001001) \tab (001011) \tab (001010) \tab (001110) \tab (001111) \tab (001101) \tab (001100) \cr
#'   \tab \tab \tab \tab \tab   \cr
#'  -7-7i  \tab -5-7i \tab +3-7i  \tab -1-7i \tab +1-7i  \tab +3-7i \tab +5-7i  \tab +7-7i \cr
#'  (000000) \tab (000001) \tab (000011) \tab (000010) \tab (000110) \tab (000111) \tab (000101) \tab (000100)
#' }
#' Reference: E. Agrell, J Lassing, E. Strom, and T. Ottosson, Gray Coding for Multilevel Constellations In Gaussian Noise, IEEE Transactions on Communications, Vol. 53, No. 1, January 2007
#' @param bits - received vector of bits (0's and 1's).
#' @param Ns - N samples per symbol (default, Ns = 1)
#' @param p - a vector defining the pulse shape of the transmitted waveform (default, p = 1)
#' @examples
#' M=64
#' Es=42
#' Eb = Es/log2(M)
#' Nsymbols=1000
#' Nbits=log2(M)*Nsymbols
#' bits <- sample(0:1,Nbits, replace=TRUE)
#' s <- f64qammod(bits)
#' @return Returns a complex vector of  64-QAM symbols. If Ns > 1
#' then the returned signal is shaped with pulse shape, p.
#' @family rwirelesscom functions
#' @export
#'
f64qammod <- function(bits,Ns=1,p=1) {
  # x quad1 = 1+1i 1+3i 1+5i 1+7i 3+1i 3+3i 3+5i 3+7i 5+1i 5+3i 5+5i 5+7i 7+1i 7+3i 7+5i 7+7i
  # Es = sum(Re(x*Conj(x)))/16 = 42
  # receive bits
  # transform to symbolbits (Nsym rows) length(bits)/4 x Log2(M) cols
  # transform to symbolcodes
  # transform to xreal, yimag
  M=64
  Nsymbols = length(bits) %/% log2(M)
  symbolbits<-matrix(bits,Nsymbols,log2(M),byrow=TRUE)
  m64=cbind(c(0,0,0,0,0,1),c(0,0,0,0,2,0),c(0,0,0,4,0,0),c(0,0,8,0,0,0),c(0,16,0,0,0,0),c(32,0,0,0,0,0))
  symbolbits = symbolbits %*% m64
  symbolcodes <- apply(symbolbits,1,sum)
  s <- sapply(symbolcodes,ft64qambitmap)
  if (Ns >= 1) {
    Nsymbols=length(s)
    x <-as.vector(rbind(s,matrix(c(0),Ns-1,Nsymbols)))
    s= convolve(x,p,type="open")
  }
  return(s)
}

ft64qambitmap <- function(x) {
  #   imag      real
                                                      # b5b6b4b3    b2b1b0
  if (x == 0)      r = complex(real=-7, imaginary=-7)   #  000     000
  else if (x == 1) r = complex(real=-5, imaginary=-7)   #  000     001
  else if (x == 3) r = complex(real=-3, imaginary=-7)   #  000     011
  else if (x == 2) r = complex(real=-1, imaginary=-7)   #  000     010
  else if (x == 6) r = complex(real=+1, imaginary=-7)   #  000     110
  else if (x == 7) r = complex(real=+3, imaginary=-7)   #  000     111
  else if (x == 5) r = complex(real=+5, imaginary=-7)   #  000     101
  else if (x == 4) r = complex(real=+7, imaginary=-7)   #  000     100

  else if (x == 8) r = complex(real=-7, imaginary=-5)   #  001     000
  else if (x == 9) r = complex(real=-5, imaginary=-5)   #  001     001
  else if (x == 11) r = complex(real=-3, imaginary=-5)  #  001     011
  else if (x == 10) r = complex(real=-1, imaginary=-5)  #  001     010
  else if (x == 14) r = complex(real=+1, imaginary=-5)  #  001     110
  else if (x == 15) r = complex(real=+3, imaginary=-5)  #  001     111
  else if (x == 13) r = complex(real=+5, imaginary=-5)  #  001     101
  else if (x == 12) r = complex(real=+7, imaginary=-5)  #  001     100

  else if (x == 24) r = complex(real=-7, imaginary=-3)  #  011     000
  else if (x == 25) r = complex(real=-5, imaginary=-3)  #  011     001
  else if (x == 27) r = complex(real=-3, imaginary=-3)  #  011     011
  else if (x == 26) r = complex(real=-1, imaginary=-3)  #  011     010
  else if (x == 30) r = complex(real=+1, imaginary=-3)  #  011     110
  else if (x == 31) r = complex(real=+3, imaginary=-3)  #  011     111
  else if (x == 29) r = complex(real=+5, imaginary=-3)  #  011     101
  else if (x == 28) r = complex(real=+7, imaginary=-3)  #  011     100

  else if (x == 16) r = complex(real=-7, imaginary=-1)  #  010     000
  else if (x == 17) r = complex(real=-5, imaginary=-1)  #  010     001
  else if (x == 19) r = complex(real=-3, imaginary=-1)  #  010     011
  else if (x == 18) r = complex(real=-1, imaginary=-1)  #  010     010
  else if (x == 22) r = complex(real=+1, imaginary=-1)  #  010     110
  else if (x == 23) r = complex(real=+3, imaginary=-1)  #  010     111
  else if (x == 21) r = complex(real=+5, imaginary=-1)  #  010     101
  else if (x == 20) r = complex(real=+7, imaginary=-1)  #  010     100

  else if (x == 48) r = complex(real=-7, imaginary=+1)  #  110     000
  else if (x == 49) r = complex(real=-5, imaginary=+1)  #  110     001
  else if (x == 51) r = complex(real=-3, imaginary=+1)  #  110     011
  else if (x == 50) r = complex(real=-1, imaginary=+1)  #  110     010
  else if (x == 54) r = complex(real=+1, imaginary=+1)  #  110     110
  else if (x == 55) r = complex(real=+3, imaginary=+1)  #  110     111
  else if (x == 53) r = complex(real=+5, imaginary=+1)  #  110     101
  else if (x == 52) r = complex(real=+7, imaginary=+1)  #  110     100

  else if (x == 56) r = complex(real=-7, imaginary=+3)  #  111     000
  else if (x == 57) r = complex(real=-5, imaginary=+3)  #  111     001
  else if (x == 59) r = complex(real=-3, imaginary=+3)  #  111     011
  else if (x == 58) r = complex(real=-1, imaginary=+3)  #  111     010
  else if (x == 62) r = complex(real=+1, imaginary=+3)  #  111     110
  else if (x == 63) r = complex(real=+3, imaginary=+3)  #  111     111
  else if (x == 61) r = complex(real=+5, imaginary=+3)  #  111     101
  else if (x == 60) r = complex(real=+7, imaginary=+3)  #  111     100

  else if (x == 40) r = complex(real=-7, imaginary=+5)  #  101     000
  else if (x == 41) r = complex(real=-5, imaginary=+5)  #  101     001
  else if (x == 43) r = complex(real=-3, imaginary=+5)  #  101     011
  else if (x == 42) r = complex(real=-1, imaginary=+5)  #  101     010
  else if (x == 46) r = complex(real=+1, imaginary=+5)  #  101     110
  else if (x == 47) r = complex(real=+3, imaginary=+5)  #  101     111
  else if (x == 45) r = complex(real=+5, imaginary=+5)  #  101     101
  else if (x == 44) r = complex(real=+7, imaginary=+5)  #  101     100

  else if (x == 32) r = complex(real=-7, imaginary=+7)  #  100     000
  else if (x == 33) r = complex(real=-5, imaginary=+7)  #  100     001
  else if (x == 35) r = complex(real=-3, imaginary=+7)  #  100     011
  else if (x == 34) r = complex(real=-1, imaginary=+7)  #  100     010
  else if (x == 38) r = complex(real=+1, imaginary=+7)  #  100     110
  else if (x == 39) r = complex(real=+3, imaginary=+7)  #  100     111
  else if (x == 37) r = complex(real=+5, imaginary=+7)  #  100     101
  else if (x == 36) r = complex(real=+7, imaginary=+7)  #  100     100

  else {print(x)
    r=100}
  return(r)
}

#' 64-QAM Demodulator
#'
#' Receives a vector of complex values, r, corresponding to a
#' 64-QAM modulated signal transmitted through a communications channel
#' (e.g., signal plus noise). Each received 64-QAM symbol (one symbol per r sample) is decoded into
#' 6 bits. The 64-QAM symbol decision regions are defined with respect to
#' the constellation generated by f64qammod(), where in-phase and quadrature
#' constellation points take on values -7, -5, -3, -1, +1, +3, +5, +7, respectively.
#' @param r - complex valued input vector
#' @family rwirelesscom functions
#' @return a vector of 1's and 0's, 6 bits per input element (64-QAM symbol)
#' @examples
#' M=64
#' Es=42
#' Eb = Es/log2(M)
#' Nsymbols=1000
#' Nbits=log2(M)*Nsymbols
#' bits <- sample(0:1,Nbits, replace=TRUE)
#' s <- f64qammod(bits)
#' EbNodB=12
#' No = Eb/(10^(EbNodB/10))
#' n <- fNo(Nsymbols,No,type="complex")
#' r <- s+n
#' bitsr <- f64qamdemod(r)
#' biterrs<-bits[bitsr!=bits]
#' Pberr=length(biterrs)/length(bits)
#' @export
f64qamdemod <- function(r) {
  r2 <- sapply(r,fr64qambitmap)
  bitsr <- as.vector(matrix(r2,1,length(r2)))
  return(bitsr)
}

fr64qambitmap <- function(r) {
  if (Im(r) < -6 ) {
    if (Re(r) < -6) bits=c(0,0,0,0,0,0)
    else if (Re(r) < -4) bits=c(0,0,0,0,0,1)
    else if (Re(r) < -2) bits=c(0,0,0,0,1,1)
    else if (Re(r) < 0)  bits=c(0,0,0,0,1,0)
    else if (Re(r) < 2)  bits=c(0,0,0,1,1,0)
    else if (Re(r) < 4)  bits=c(0,0,0,1,1,1)
    else if (Re(r) < 6)  bits=c(0,0,0,1,0,1)
    else                 bits=c(0,0,0,1,0,0)
  }
  else if (Im(r) < -4 ) {
    if (Re(r) < -6) bits=c(0,0,1,0,0,0)
    else if (Re(r) < -4) bits=c(0,0,1,0,0,1)
    else if (Re(r) < -2) bits=c(0,0,1,0,1,1)
    else if (Re(r) < 0)  bits=c(0,0,1,0,1,0)
    else if (Re(r) < 2)  bits=c(0,0,1,1,1,0)
    else if (Re(r) < 4)  bits=c(0,0,1,1,1,1)
    else if (Re(r) < 6)  bits=c(0,0,1,1,0,1)
    else                 bits=c(0,0,1,1,0,0)
  }
  else if (Im(r) < -2 ) {
    if (Re(r) < -6) bits=c(0,1,1,0,0,0)
    else if (Re(r) < -4) bits=c(0,1,1,0,0,1)
    else if (Re(r) < -2) bits=c(0,1,1,0,1,1)
    else if (Re(r) < 0)  bits=c(0,1,1,0,1,0)
    else if (Re(r) < 2)  bits=c(0,1,1,1,1,0)
    else if (Re(r) < 4)  bits=c(0,1,1,1,1,1)
    else if (Re(r) < 6)  bits=c(0,1,1,1,0,1)
    else                 bits=c(0,1,1,1,0,0)
  }
  else if (Im(r) < 0 ) {
    if (Re(r) < -6) bits=c(0,1,0,0,0,0)
    else if (Re(r) < -4) bits=c(0,1,0,0,0,1)
    else if (Re(r) < -2) bits=c(0,1,0,0,1,1)
    else if (Re(r) < 0)  bits=c(0,1,0,0,1,0)
    else if (Re(r) < 2)  bits=c(0,1,0,1,1,0)
    else if (Re(r) < 4)  bits=c(0,1,0,1,1,1)
    else if (Re(r) < 6)  bits=c(0,1,0,1,0,1)
    else                 bits=c(0,1,0,1,0,0)
  }
  else if (Im(r) < 2 ) {
    if (Re(r) < -6)bits=c(1,1,0,0,0,0)
    else if (Re(r) < -4) bits=c(1,1,0,0,0,1)
    else if (Re(r) < -2) bits=c(1,1,0,0,1,1)
    else if (Re(r) < 0)  bits=c(1,1,0,0,1,0)
    else if (Re(r) < 2)  bits=c(1,1,0,1,1,0)
    else if (Re(r) < 4)  bits=c(1,1,0,1,1,1)
    else if (Re(r) < 6)  bits=c(1,1,0,1,0,1)
    else                 bits=c(1,1,0,1,0,0)
  }
  else if (Im(r) < 4 ) {
    if (Re(r) < -6) bits=c(1,1,1,0,0,0)
    else if (Re(r) < -4) bits=c(1,1,1,0,0,1)
    else if (Re(r) < -2) bits=c(1,1,1,0,1,1)
    else if (Re(r) < 0)  bits=c(1,1,1,0,1,0)
    else if (Re(r) < 2)  bits=c(1,1,1,1,1,0)
    else if (Re(r) < 4)  bits=c(1,1,1,1,1,1)
    else if (Re(r) < 6)  bits=c(1,1,1,1,0,1)
    else                 bits=c(1,1,1,1,0,0)
  }
  else if (Im(r) < 6 ) {
    if (Re(r) < -6)  bits=c(1,0,1,0,0,0)
    else if (Re(r) < -4) bits=c(1,0,1,0,0,1)
    else if (Re(r) < -2) bits=c(1,0,1,0,1,1)
    else if (Re(r) < 0)  bits=c(1,0,1,0,1,0)
    else if (Re(r) < 2)  bits=c(1,0,1,1,1,0)
    else if (Re(r) < 4)  bits=c(1,0,1,1,1,1)
    else if (Re(r) < 6)  bits=c(1,0,1,1,0,1)
    else                 bits=c(1,0,1,1,0,0)
  }
  else if (Im(r) >=6 ) {
    if (Re(r) < -6) bits=c(1,0,0,0,0,0)
    else if (Re(r) < -4) bits=c(1,0,0,0,0,1)
    else if (Re(r) < -2) bits=c(1,0,0,0,1,1)
    else if (Re(r) < 0)  bits=c(1,0,0,0,1,0)
    else if (Re(r) < 2)  bits=c(1,0,0,1,1,0)
    else if (Re(r) < 4)  bits=c(1,0,0,1,1,1)
    else if (Re(r) < 6)  bits=c(1,0,0,1,0,1)
    else                 bits=c(1,0,0,1,0,0)
  }
  return(bits)
}

#' IQ Scatter Plot
#'
#' A convenience function that plots a scatter diagram of Im(r) vs. Re(r). The function is
#' useful for visualizing constellations such as M-PSK or M-QAM.
#'
#' @param  r - complex or real valued vector
#' @examples
#' M=8
#' Es=1
#' Eb = Es/log2(M)
#' Nsymbols=10000
#' Nbits=log2(M)*Nsymbols
#' bits <- sample(0:1,Nbits, replace=TRUE)
#' s <- f8pskmod(bits)
#' EbNodB=7
#' No = Eb/(10^(EbNodB/10))
#' n <- fNo(Nsymbols,No,type="complex")
#' r <- s+n
#' iqscatterplot(r)
#' @family rwirelesscom functions
#' @export

iqscatterplot <- function(r) {
  X1 <- X2 <- NULL
  ggplot(data.frame(cbind(Re(r),Im(r))), aes(x=X1, y=X2)) + geom_point(size=1)+ xlab("In-Phase") + ylab("Quadrature") + theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face="bold"))
}


#' IQ Density Plot
#'
#' A convenience function to plot a density function of a vector containing the in-phase and
#' quadrature signal (plus noise).
#' @param r - complex or real valued vector
#' @param iq - if iq = "r" (default) then plot density of Re(r) else if iq = "q" then plot density of Im(r)
#' @examples
#' M=4
#' Es=1
#' Eb = Es/log2(M)
#' Nsymbols=1000
#' Nbits=log2(M)*Nsymbols
#' bits <- sample(0:1,Nbits, replace=TRUE)
#' s <- fqpskmod(bits)
#' EbNodB=4
#' No = Eb/(10^(EbNodB/10))
#' n <- fNo(Nsymbols,No,type="complex")
#' r <- s+n
#  iqdensityplot(r,iq="q")
#'  @family rwirelesscom functions
#' @export
iqdensityplot <- function(r,iq="r") {
  ..density.. <- NULL
  if (iq=="r") { # Real Part
    ggplot(data.frame(r), aes(x=Re(r))) +  geom_histogram(binwidth=0.05, colour="black", fill=NA,position="identity",aes(y=..density..)) +
        xlab("") + theme_bw() + theme(axis.text = element_text( size=16), axis.title=element_text(size=16,face="bold"))
  } else {    #Imaginary Part
    ggplot(data.frame(r), aes(x=Im(r))) +  geom_histogram(binwidth=0.05, colour="black", fill=NA,position="identity",aes(y=..density..)) +
       xlab("") + theme_bw() + theme(axis.text = element_text( size=16), axis.title=element_text(size=16,face="bold"))
  }
}

# http://www.r-bloggers.com/matlab-style-stem-plot-with-r/
# An example using stemplot function below

#' Stem Plot
#'
#' Receives a vector of x and  y values and plots a stemplot (line and "points").
#'
#' #' Reference: M, Pastell, http://www.r-bloggers.com/matlab-style-stem-plot-with-r
#' @param x - vector of x axis points
#' @param y - vector of y axis points
#' @param pch - plot character default = 19
#' @param linecol - default line color = 1 (black)
#' @param linew - default line width = 1
#' @param ... - graphical environment parameters are input to stemplot
#' @family rwirelesscom functions
#' @examples
#' x <- seq(-3*pi, 3*pi, by = 0.2)
#' y <- sinc(x)
#' stemplot(x/pi,y,ylim=c(-0.3,1.1), pch=19, cex=0.3, ylab="y", xlab="x/pi")
#' @export
stemplot <- function(x,y,pch=16,linecol=1,linew=1,...){
  if (missing(y)){
    y = x
    x = 1:length(x) }
  plot(x,y,pch=pch,yaxs="i",...)
  for (i in 1:length(x)){
    lines(c(x[i],x[i]), c(0,y[i]),col=linecol,lwd=linew)
  }
  lines(c(x[1]-2,x[length(x)]+2), c(0,0),col=linecol,lwd=linew)
}

#' Eye Diagram
#'
#' Receives a vector x of real or complex points and plots the Re or Im part of x in the form of an "eyediagram."
#' The symbol period is indicated by the input Ns (samples per symbol) and horizontal sweep Np indicates
#' the the number of symbol periods to plot along the horizontal access. The eyediagram
#' is useful for evaluating inter-symbol interference associated with the in-phase or quadrature part
#' of a modulated signal.
#' @param x - vector of real or complex points
#' @param Ns - number of samples per symbol period
#' @param Np - number of symbol periods to plot along the horizontal axis
#' @param No - offset (n points) alignment of the eyediagram along the horizontal asis
#' @param iq - parameter indicates whether to plot the in-phase Re(x) (iq="r" default) or quadrature Im(x) (iq="q")
#' @param pch - Graphical parameter pch (plotting character) set to 19 by default ("point")
#' @param cex - Graphical parameter cex magnificaiton of plotting symbols relative to 1 default set to 0.1.
#' @param ... - graphical environment parameters are input to eyediagram
#' @family rwirelesscom functions
#' @examples
#'
#' # Step 1: generate random set of bits
#' M=4
#' Nsymbols=10000
#' Nbits=log2(M)*Nsymbols
#' bits <- sample(0:1,Nbits, replace=TRUE)
#'
#' # Step 2: Generate a BPSK modulated signal including raised cosine
#' #  pulse shaping sampled at 64 samples per symbol period and roll-off
#' #  factor of 0.5.
#' Ns=64
#' B=0.5
#' hx=seq(-4*Ns,4*Ns,by=1)
#' h=rcosine(hx,B,Ns)
#' s <- fqpskmod(bits,Ns,h)
#'
#' # Step 3: Plot the transmitted waveform with the eyediagram function.
#' #    Remove the initial 1000 and tail points for a cleaner diagram, i.e.
#' #    without startup and tail artifacts
#' Np=3
#' No=27
#'
#' # real part (in-phase)
#' eyediagram(s[1000:60000],Ns,Np,No,iq="r",xlab="Ts",ylab="I")
#' # imaginary part (quadrature)
#' eyediagram(s[1000:60000],Ns,Np,No,iq="q",xlab="Ts",ylab="Q")
#' @export
eyediagram <- function(x,Ns=1,Np=3,No=1,iq="r",pch=19,cex=0.1,...) {
  if (iq=="r")  y=Re(x[No:length(x)])
  else if (iq=="q") y=Im(x[No:length(x)])

  xx=(seq(1:length(y)) %% (Ns*Np))/Ns

  plot(xx,y,pch=19,cex=0.1,...)
}
