
# The distances in TSclust are calculated by means of the following wrapper functions. The functions of the TSclust package are directly used.

# Autocorrelation based similarity

ACFDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.ACF(x, y, ...))}, 
  error=function(e) {print(e); NA})  
}

# Partial autocorrelation based similarity

PACFDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.PACF(x, y, ...))},
  error=function(e) {print(e); NA})  
}

# Dissimilarity Based on LPC Cepstral Coefficients

ARLPCCepsDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.AR.LPC.CEPS(x, y, ...))},
  error=function(e) {print(e); NA})  
}

# Model-based Dissimilarity Proposed by Maharaj (1996, 2000)
ARMahDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  diss.AR.MAH(x, y, ...)},
  error=function(e) {print(e); NA})  
  }

ARMahStatisticDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.AR.MAH(x, y, ...)$statistic)},
  error=function(e) {print(e); NA})   
}

ARMahPvalueDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.AR.MAH(x, y, ...)$p_value)},
  error=function(e) {print(e); NA})   
}

# Model-based Dissimilarity Measure Proposed by Piccolo (1990)

ARPicDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.AR.PIC(x, y, ...))},
  error=function(e) {print(e); NA})   
}

# Compression-based Dissimilarity measure

CDMDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.CDM(x, y, ...))},
  error=function(e) {print(e); NA})   
}

# Complexity-Invariant Distance Measure For Time Series

CIDDistance <- function(x, y) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.CID(x, y))},
  error=function(e) {print(e); NA})   
}

# Correlation-based Dissimilarity

CorDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.COR(x, y, ...))},
  error=function(e) {print(e); NA})   
}

# Dissimilarity Index Combining Temporal Correlation and Raw Values
# Behaviors

CortDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.CORT(x, y, ...))},
  error=function(e) {print(e); NA})   
}

# Dissimilarity for Time Series Based on Wavelet Feature Extraction
WavDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.DWT(rbind(x, y)))},
  error=function(e) {print(e); NA})   
}


# Integrated Periodogram Based Dissimilarity
IntPerDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.INT.PER(x, y, ...))},
  error=function(e) {print(e); NA})   
}

# Periodogram Based Dissimilarity

PerDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.PER(x, y, ...))},
  error=function(e) {print(e); NA})   
}


# Symbolic Aggregate Aproximation related functions

MindistSaxDistance <- function(x, y, w, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.MINDIST.SAX(x, y, w, ...))},
  error=function(e) {print(e); NA})   
}

# Normalized Compression Distance

NCDDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.NCD(x, y, ...))},
  error=function(e) {print(e); NA})   
}

# Dissimilarity Measure Based on Nonparametric Forecast

PredDistance <- function(x, y, h, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.PRED(x, y, h, ...)$L1dist)},
  error=function(e) {print(e); NA})   
}

# Dissimilarity based on the Generalized Likelihood Ratio Test

SpecGLKDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss(rbind(x,y), "SPEC.GLK", ...))},
  error=function(e) {print(e); NA})   
}

# Dissimilarity Based on the Integrated Squared Difference between the
# Log-Spectra

SpecISDDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.SPEC.ISD(x, y, ...))},
  error=function(e) {print(e); NA})   
}


# General Spectral Dissimilarity Measure Using Local-Linear Estima-
# tion of the Log-Spectra

SpecLLRDistance <- function(x, y, ...) {
  # If there is an error, NA is returned and the error message 
  # is printed. This enables executing in batch mode, without stops.
  tryCatch ({
  as.numeric(diss.SPEC.LLR(x, y, ...))},
  error=function(e) {print(e); NA})   
}

