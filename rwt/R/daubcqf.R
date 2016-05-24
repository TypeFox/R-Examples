###
### $Id: daubcqf.R 22 2014-06-20 20:59:33Z plroebuck $
### Computes Daubechies' scaling and wavelet filters (normalized to sqrt(2)).
###

PHASE.MINIMUM <- 'min'
PHASE.MID     <- 'mid'
PHASE.MAXIMUM <- 'max'


##-----------------------------------------------------------------------------
daubcqf <- function(N, type=PHASE.MINIMUM) {
#%    [h_0,h_1] = daubcqf(N,TYPE);
#%
#%    Function computes the Daubechies' scaling and wavelet filters
#%    (normalized to sqrt(2)).
#%
#%    Input:
#%       N    : Length of filter (must be even)
#%       TYPE : Optional parameter that distinguishes the minimum phase,
#%              maximum phase and mid-phase solutions. If no argument is
#%              specified, the minimum phase solution is used.
#%
#%    Output:
#%       h_0 : Minimal phase Daubechies' scaling filter
#%       h_1 : Minimal phase Daubechies' wavelet filter
#%
#%    Example:
#%       N = 4;
#%       TYPE = 'min';
#%       [h_0,h_1] = daubcqf(N,TYPE)
#%       h_0 = 0.4830 0.8365 0.2241 -0.1294
#%       h_1 = 0.1294 0.2241 -0.8365 0.4830
#%
#%    Reference: "Orthonormal Bases of Compactly Supported Wavelets",
#%                CPAM, Oct.89
#%

    if (!(is.numeric(N) && length(N) == 1)) {
        stop(sprintf("argument %s must be numeric scalar",
                     sQuote("N")))
    }

    if (N <= 0) {
        stop(sprintf("filter length cannot be negative (%s = %s)",
                     sQuote("N"),
                     N))
    }

    if (type == PHASE.MID) {
        stop(sprintf("code not implemented for (%s = PHASE.MID)",
                     sQuote("type")))
    }

    h.0 <- .getDaubechiesCoefficients(N)

    if (type == PHASE.MAXIMUM) {
        h.0 <- matlab::fliplr(t(as.matrix(h.0)))
    }

    if (abs(sum(h.0 ^ 2)) - 1 > 1e-4) {
        stop(sprintf("numerically unstable for this value of %s",
                     sQuote("N")))
    }

    h.1 <- matlab::rot90(as.matrix(h.0), 2)
    h.1[seq(1, N, by = 2)] <- -h.1[seq(1, N, by = 2)]

    list(h.0 = as.vector(h.0),
         h.1 = as.vector(h.1))
}

