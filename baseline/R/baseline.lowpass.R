baseline.lowpass <- function(spectra, steep = 2, half = 5) {
## Low-pass filter based on fast fourier transform
## as presented by Atakan, Blass and Jennings
## Coded by Kristian Hovde Liland and Bjørn-Helge Mevik
## $Id: baseline.lowpass.R 170 2011-01-03 20:38:25Z bhm $
#
# INPUT:
# spectra - rows of spectra
# steep   - steepness of filter curve
# half    - halfway point of filter curve
#
# OUTPUT:
# corrected - low-pass filtered spectra (beware of scale changes)

    if (!is.matrix(spectra)) spectra <- t(spectra)        # In case of single spectrum
    np <- ncol(spectra)

    ## Low-pass filter
    f <- 0.5 * (1 + tanh(steep * (1:np - half))) *
        (1 + tanh(steep * (rev(1:np - (half - 1)))))

    ## Fast Fourier Transform
    z <- mvfft(t(spectra))

    ## Filtering modulus and transforming back to spectra
    z <- Re(mvfft(f * z, inverse = TRUE))
    corrected <- t(z) * rowMeans(spectra) / colMeans(z)
    list(baseline = spectra - corrected, corrected = corrected)
}
