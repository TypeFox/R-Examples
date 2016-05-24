TK95 <- function(N=1000, alpha=1.5) {
    # length of time series, parameter of power law (pink noise: alpha=1)
    # Fourier frequencies
    f <- seq(from=0, to=pi, length.out=(N/2+1))[-c(1,(N/2+1))]
    # power law
    f_ <- 1/f^alpha
    # real of Fourier transform
    RW <- sqrt(0.5*f_)*rnorm(N/2-1)
    # imaginary of Fourier transform
    IW <- sqrt(0.5*f_)*rnorm(N/2-1)
    fR <- complex(length.out=N, real=c(rnorm(1), RW, rnorm(1), RW[(N/2-1):1]), imaginary=c(0, IW, 0, -IW[(N/2-1):1]))
    # (
    # complex numbers to be retransformed
    # sequence of frequencies: 0,2pi/N, 2*2pi/N,...,pi,...,2pi-1/N
    # two frequencies f_1 = pi + c , f_2 = pi - c shall be complex conjugates
    # frequencies at 0 or pi shall have imaginary=0
    # )
    # transform to times
    reihe<-fft(fR, inverse=TRUE)
    # series as reals
    return(Re(reihe))
}
