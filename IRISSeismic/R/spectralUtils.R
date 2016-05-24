##
##    Utility functions for the 'IRISSeismic' package.
##
##    Copyright (C) 2013  Mazama Science, Inc.
##    by Jonathan Callahan, jonathan@mazamascience.com
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


################################################################################
# Helpful documents explaining spectral analysis and R's spec.pgram function.
################################################################################

# "Spectral Analysis -- Smoothed Periodogram Method"
#
#   http://www.ltrr.arizona.edu/~dmeko/notes_6.pdf
#
#   Simple explanation of DFT and smoothed periodogram with an explanation of 
#   Daniell smoothers:

################################################################################
# crossSpectrum is a function specifically for calculating cross-spectral values.
# This function marries functionality from R's spec.pgram with MATLAB's pwelch.
################################################################################

# Most of this function comes directly from R-3.0.1/src/library/stats/R/spectrum.R
#
# Other features mimic the code in Octave's pwelch() function:
#
# http://octave-signal.sourcearchive.com/documentation/1.0.7/pwelch_8m-source.html

crossSpectrum <- function (x, spans = NULL, kernel = NULL, taper = 0.1,
                           pad = 0, fast = TRUE,
                           demean = FALSE, detrend = TRUE,
                           na.action = stats::na.fail) {
  
  # BEGIN: Identical to spec.pgram ---------------------------------------------
  
  ## Estimate spectral density from (smoothed) periodogram.
  series <- deparse(substitute(x))
  x <- na.action(stats::as.ts(x))
  xfreq <- stats::frequency(x)
  x <- as.matrix(x)
  N <- N0 <- nrow(x)
  nser <- ncol(x)
  if(!is.null(spans)) { # allow user to mistake order of args
    kernel <- {
      if(stats::is.tskernel(spans)) spans else
        kernel("modified.daniell", spans %/% 2)
    }
  }
  if(!is.null(kernel) && !stats::is.tskernel(kernel)) {
    stop("must specify 'spans' or a valid kernel")
  }
  if (detrend) {
    t <- 1L:N - (N + 1)/2
    sumt2 <- N * (N^2 - 1)/12
    for (i in 1L:ncol(x))
      x[, i] <- x[, i] - mean(x[, i]) - sum(x[, i] * t) * t/sumt2
  } else if (demean) {
    x <- sweep(x, 2, colMeans(x), check.margin=FALSE)
  }
  ## apply taper:
  x <- stats::spec.taper(x, taper)
  ## to correct for tapering: Bloomfield (1976, p. 194)
  ## Total taper is taper*2
  u2 <- (1 - (5/8)*taper*2)
  u4 <- (1 - (93/128)*taper*2)
  if (pad > 0) {
    x <- rbind(x, matrix(0, nrow = N * pad, ncol = ncol(x)))
    N <- nrow(x)
  }
  NewN <- if(fast) stats::nextn(N) else N
  x <- rbind(x, matrix(0, nrow = (NewN - N), ncol = ncol(x)))
  N <- nrow(x)
  Nspec <- floor(N/2)
  freq <- seq.int(from = xfreq/N, by = xfreq/N, length.out = Nspec)

  # END:  Identical to spec.pgram ----------------------------------------------

  # NOTE:  The rest of this functions contains lines from the original spec.pgram that
  # NOTE:  are commented out and replaced by code using notation more akin to that
  # NOTE:  found in the Octave pwelch() function.
  
  if (nser > 2) {
    stop("crossSpectrum: x contains",nser,"timeseries, max 2 allowed.")
  }
  
#   xfft <- mvfft(x) 
#   pgram <- array(NA, dim = c(N, ncol(x), ncol(x)))
#   for (i in 1L:ncol(x)) {
#     for (j in 1L:ncol(x)) { # N0 = #{non-0-padded}
#       pgram[, i, j] <- xfft[, i] * Conj(xfft[, j])/(N0*xfreq)
#       ## value at zero is invalid as mean has been removed, so interpolate:
#       pgram[1, i, j] <- 0.5*(pgram[2, i, j] + pgram[N, i, j])
#     }
#   }
  
  fft_x <- as.complex(stats::fft(x[,1]))  
  # Periodogram
  Pxx <- fft_x * Conj(fft_x)/(N0*xfreq)
  Pxx[1] <- 0.5*(Pxx[2] + Pxx[N])

  if (nser == 1) {
    Pyy <- Pxy <- NULL
    # NOTE:  Pyx is identical to Pxy
  } else {
    fft_y <- as.complex(stats::fft(x[,2]))
    Pyy <- fft_y * Conj(fft_y)/(N0*xfreq)
    Pyy[1] <- 0.5*(Pyy[2] + Pyy[N])
    # Cross-power spectrum
    Pxy <- fft_x * Conj(fft_y)/(N0*xfreq)
    Pxy[1] <- 0.5*(Pxy[2] + Pxy[N]) # TODO:  Is this correct?
   }
  
#   if(!is.null(kernel)) {
#     for (i in 1L:ncol(x)) {
#       for (j in 1L:ncol(x)) {
#         pgram[, i, j] <- stats::kernapply(pgram[, i, j], kernel, circular = TRUE)
#       }
#     }
#     df <- df.kernel(kernel)
#     bandwidth <- bandwidth.kernel(kernel)
#   } else { # raw periodogram
#     df <- 2
#     bandwidth <- sqrt(1/12)
#   }
  
  # NOTE:  We won't be returning 'df' or 'bandwidth'
  if (!is.null(kernel)) {
    Pxx <- stats::kernapply(Pxx, kernel, circular=TRUE)
    if (nser == 2) {
      Pyy <- stats::kernapply(Pyy, kernel, circular=TRUE)
      Pxy <- stats::kernapply(Pxy, kernel, circular=TRUE) # TODO:  Is this correct?
    }
  }
  
#   df <- df/(u4/u2^2)
#   df <- df  * (N0 / N) ## << since R 1.9.0
#   bandwidth <- bandwidth * xfreq/N
#   pgram <- pgram[2:(Nspec+1),,, drop=FALSE]
#   spec <- matrix(NA, nrow = Nspec, ncol = nser)
#   for (i in 1L:nser) spec[, i] <- Re(pgram[1L:Nspec, i, i])
  
  # Force spectra to be real
  Pxx <- Pxx[2:(Nspec+1)]
  spec1 <- pracma::Real(Pxx[1:Nspec])
  if (nser == 2) {
    Pyy <- Pyy[2:(Nspec+1)]
    Pxy <- Pxy[2:(Nspec+1)] # TODO:  Is this correct?
    # NOTE:  Pyx is identical to Pxy
    spec2 <- pracma::Real(Pyy[1:Nspec])  
  }

  # TODO:  Would this be the place to multiply by a factor of 2 to
  # TODO:  get a 'one-sided' spectrum?
  
#   if (nser == 1) {
#     coh <- phase <- NULL
#   } else {
#     coh <- phase <- matrix(NA, nrow = Nspec, ncol = nser * (nser - 1)/2)
#     for (i in 1L:(nser - 1)) {
#       for (j in (i + 1):nser) {
#         coh[, i + (j - 1) * (j - 2)/2] <-
#           Mod(pgram[, i, j])^2/(spec[, i] * spec[, j])
#         phase[, i + (j - 1) * (j - 2)/2] <- Arg(pgram[, i, j])
#       }
#     }
#   }
 
  if (nser == 1) {
    coh <- phase <- NULL
  } else {
    coh <- Mod(Pxy)^2 / (spec1 * spec2)
    phase <- Arg(Pxy)
  }
  
#   ## correct for tapering
#   for (i in 1L:nser) spec[, i] <- spec[, i]/u2
  
  spec1 <- spec1/u2
  spec2 <- spec2/u2
  
#   spec <- drop(spec)
#   spg.out <-
#     list(freq = freq, spec = spec, coh = coh, phase = phase,
#          kernel = kernel, df = df,
#          bandwidth = bandwidth, n.used = N, orig.n = N0,# "n.orig" = "n..."
#          series = series, snames = colnames(x),
#          method = ifelse(!is.null(kernel), "Smoothed Periodogram",
#                          "Raw Periodogram"),
#          taper = taper, pad = pad, detrend = detrend, demean = demean)
#   class(spg.out) <- "spec"
#   
#   return(spg.out)
  
  DF <- data.frame(freq=freq, spec1=spec1, spec2=spec2, coh=coh, phase=phase,
                   Pxx=Pxx, Pyy=Pyy, Pxy=Pxy)

  return(DF)

}


################################################################################
# Calculate average values in bins using the McNamara algorithm:
#   http://pubs.usgs.gov/of/2005/1438/pdf/OFR-1438.pdf
################################################################################

McNamaraBins <- function(df, loFreq=.005, hiFreq=10, alignFreq=0.1) {
   
  # Frequencies for binning will be generated at 1/8 octave intervals aligned to alignFreq
  
  # Frequency smoothe the averaged spectral values over 1-octave intervals at 1/8-octave increments, 
  #     (reducing # frequency samples by factor of 169; include 10 Hz as one of the geometric mean f values to sync f sampling)
  # 
  # Additional notes from Mary Templeton on 2013-05-31:
  #    The frequency smoothing (second to last step) involves averaging power values over an octave using a geometric mean 
  #    and repeating at 1/8 octave intervals in frequency (to refine my wording above).
  
  # Define "center" frequencies at 1/8 octave intervals from 10 Hz down to loFreq
  if (alignFreq >= hiFreq) {
    octaves <- seq(log2(alignFreq),log2(loFreq),-0.125)
    octaves <- octaves[octaves <= log2(hiFreq)]
  } else {
    loOctaves <- seq(log2(alignFreq),log2(loFreq),-0.125)
    hiOctaves <- seq(log2(alignFreq),log2(hiFreq),0.125)
    octaves <- sort(unique(c(loOctaves,hiOctaves)))
  }
  binFreq <- sort(2^octaves)

  # Set up empty dataframe for return values
  DF <- as.data.frame(matrix(nrow=length(binFreq),ncol=ncol(df)))
  names(DF) <- names(df)
  
  # Define bins starting at zero so that freq #1 has bin 0 to the left and bin 1 to the right
  halfOctaveBelow <- 2^(log2(loFreq)-0.5)
  halfOctaveAbove <- 2^(log2(hiFreq)+0.5)
  
  # Use .bincode to assign every frequency to a bin
  # NOTE:  First bin begins at zero to ensure that the bin numbering aligns with the frequency numbering
  bins <- .bincode(df$freq, c(0,binFreq,halfOctaveAbove)) - 1
  
  # NOTE:  If the user says they want hiFreq=10 when max(freq) is 20 or higher, that should
  # NOTE:  be allowed.  Values above maxFreq will be assigned NA rather than a bin number and
  # NOTE:  will simply be ignored.  Frequency values below loFreq will be similarly ignored.
  
  # Loop through the bin numbers and calculate the geometric mean of all values in
  # an octave (8 bins) centered on that freq. For example, freq #1 should average all
  # frequencies in bins 0,1,2,3,4; freq #6 should have bins 2,3,4,5,6,7,8,9
  for (i in seq(length(binFreq))){
    
    # NOTE:  The lowest bin (bin 0 ~= DC) has values that corrupt the output and is ignored.
    # NOTE:  Instead, we start with bin 1.
    # NOTE:  See the extended email exchange between Rob C, Mary T and Jonathan C. in May 2014.
    
    loBin <- max(1,i-4)
    hiBin <- min(i+3,max(bins,na.rm=TRUE))
    
    # Determine which values are in this octave
    indices <- which(bins %in% seq(loBin,hiBin))
    
    # For eaach column of the dataframe, calculate the mean value for this bin.
    for (j in seq(ncol(df))) {
      DF[i,j] <- sum(df[indices,j])/length(indices)
    }
    
  } # End of binFreq loop

  # Replace averaged freq values with pre-determined ones
  DF$freq <- binFreq
  
  return(DF)
}


################################################################################
# Calculate the PSD using the McNamara algorithm:
#   http://pubs.usgs.gov/of/2005/1438/pdf/OFR-1438.pdf
################################################################################

McNamaraPSD <- function(tr, loFreq=.005, hiFreq=10, alignFreq=0.1, binned=TRUE) {
  
  # NOTE:  Each incoming trace contains datapoints. After truncation to N datapoints, we
  # NOTE:  divide the N datapoints into smaller segments of length Nseg = 4/16 N. The first 
  # NOTE:  segment begins at 0/16 and ends at 4/16.  The 13'th segment begins at 12/16
  # NOTE:  and ends at 16/16.  The small segments overlap like this:
  # NOTE:  
  # NOTE:  1---5---9---3---
  # NOTE:   2---6---0---
  # NOTE:    3---7---1---
  # NOTE:     4---8---2---
  # NOTE:  
  # NOTE:  We implement this by calculating the truncation point and the size of each segment and looping
  # NOTE:  over the 13 segments thus avoiding creation of additional variables that take up space in memory.
  
  # Truncate segment to nearest power of 2 samples
  pow2 <- floor(log(length(tr),2))
  truncatedLength <- 2^pow2
  
  # specSum will contain the sum of the spectra returned by spec.pgram()
  specSum <- 0
  
  # Divide each truncated Z-hour segment into 13 segments with 75% overlap
  for (i in 0:12) {
    
    first <- i * truncatedLength/16 + 1
    last <- (i+4) * truncatedLength/16
    
    # NOTE:  spans=NULL means that we are not smoothing at this stage
    # The spec.pgram() function takes care of all of the following:
    #   Demean
    #   Detrend
    #   Apply 10% sine taper
    #   FFT 
    #   PSD
    #   Multiply PSD by sine taper scale factor (= 1.142857 for 10% taper)
    #   Normalize the power at each PSD frequency, multiplying it by (2*dt/Nseg) where Nseg is the number of samples in the segment
    sp <- stats::spec.pgram(stats::ts(tr@data[first:last],frequency=tr@stats@sampling_rate),
                            spans=NULL,
                            taper=0.1,
                            pad=0,
                            fast=TRUE,
                            demean=TRUE,
                            detrend=TRUE,
                            plot=FALSE)
    
    # NOTE:  R's spec.pgram() returns 'two-sided' spectra and needs to be multiplied
    # NOTE:  by 2.0 when the time series being evaluated is Real rather than Complex.
    # NOTE:  Please see these two excellent sources, especially the 'psd' vignette:
    # NOTE:  
    # NOTE:    http://www.stanford.edu/class/ee262/software/signal_tb.pdf
    # NOTE:    http://cran.r-project.org/web/packages/psd/vignettes/normalization.pdf
    specSum <- specSum + 2 * sp$spec
    
  }
  
  # Average
  sp$spec <- specSum / 13
  
  # Bin specturm if requested ----------
  
  if (binned) {
    df <- data.frame(freq=sp$freq,spec=sp$spec)
    DF <- McNamaraBins(df,loFreq,hiFreq,alignFreq)
    freq <- DF$freq
    spec <- DF$spec  
  } else {
    # Return the raw spectrum
    freq <- sp$freq
    spec <- sp$spec
  }

  # Conversion to dB *AFTER* binning, not before.
  spec <- 10*log10(spec)
  
  # Create the PSD object we will return
  psd <- list(freq=freq,
              spec=spec,
              snclq=paste(tr@stats@network,tr@stats@station,tr@stats@location,tr@stats@channel,tr@stats@quality,sep="."),
              starttime=tr@stats@starttime,
              endtime=tr@stats@endtime)
  
  return(psd) 
  
}

################################################################################
# Functions for working with PSDs and PDFS
################################################################################

# NOTE:  These functions work with PSDs generated by McNamaraPSD
# NOTE:  and provide functionality for generating plots as seen in:
# NOTE:  
# NOTE:    http://pubs.usgs.gov/of/2005/1438/pdf/OFR-1438.pdf



# psdList ----------------------------------------------------------------------

# Create a psdList from an incoming seismic signal

psdList <- function(st) {
  
  # Fill any gaps in the stream with zeroes (quick because it returns immediately if no gaps are found)
  st_merged <- mergeTraces(st,fillMethod="fillZero")
  tr_merged <- st_merged@traces[[1]]
  
  # Choose chunk size based on the chanel 'band code'(=sampling rate)
  # See:  http://www.iris.edu/manuals/SEED_appA.htm
  # NOTE:  This choice was recommended by Mary Templeton
  channel <- st_merged@traces[[1]]@stats@channel
  if (stringr::str_detect(channel,"^L")) {
    Z <- 3 * 3600
    loFreq <- 0.001
    hiFreq <- 0.5 * tr_merged@stats@sampling_rate
  } else if (stringr::str_detect(channel,"^M")) {
    Z <- 2 * 3600
    loFreq <- 0.0025
    hiFreq <- 0.5 * tr_merged@stats@sampling_rate
  } else {
    Z <- 3600
    loFreq <- 0.005
    hiFreq <- 0.5 * tr_merged@stats@sampling_rate
  }
  
  # Set an alignment frequency from which octaves will be generated
  alignFreq <- 0.1
  
  # Initialize
  psdList <- list()
  index <- 1
  psdCount <- 0
  start <- tr_merged@stats@starttime
  end <- start + Z
  secsLeft <- difftime(tr_merged@stats@endtime,start,units="secs")
  
  # NOTE:  At 20 Hz we can be left with secsLeft=3599.95 which should pass the test.
  # NOTE:  Continue if at least 99% of a chunk still remains. 
  while (secsLeft >= 0.99*Z) {
    
    # Slice a chunk out of the trace
    tr <- slice(tr_merged,start,end)
    
    if (!isDC(tr)) {
      # NOTE:  Each psd has elements: freq, spec, snclq, starttime, endtime
      psds <- McNamaraPSD(tr, loFreq, hiFreq, alignFreq)
      if (! "-Inf" %in% psds$spec) {
          psdCount <- psdCount + 1
          psdList[[psdCount]] <- McNamaraPSD(tr, loFreq, hiFreq, alignFreq)
      }
    }
    
    # Increment the window with 50% overlap
    index <- index + 1
    start <- start + Z/2
    end <- end + Z/2
    secsLeft <- difftime(tr_merged@stats@endtime,start,units="secs")
    
  }
  
  return(psdList)
  
}
  

# psdList2NoiseMatrix ------------------------------------------------------------

# Create instrument corrected noiseMatrix from psdList

# NOTE:  The incoming psdList contains raw, uncorrected PSDs.
# NOTE:  The returned matrix contains corrected PSDs, one per row.
psdList2NoiseMatrix <- function(psdList) {
  
  # Need to ensure that the IrisClient object exists as R has a default
  # "iris" dataset.  See help(iris, package="datasets")
  if (class(iris) == "data.frame") {
    iris <- new("IrisClient")
  }
  
  # TODO:  psdList2NoiseMatrix should check to make sure that there is no
  # TODO:  end-of-epoch between the first PSD in the list and the last.
  # TODO:  This will only be important if there are significant changes in
  # TODO:  the instrument response.
  
  psdCount <- length(psdList)
  
  # Get information needed for getEvalresp
  # NOTE:  Each element of psdList is a list
  snclq <- as.character(psdList[[1]]$snclq)
  parts <- unlist(stringr::str_split(snclq,"[.]"))
  network <- parts[1]
  station <- parts[2]
  location <- parts[3]
  channel <- parts[4]
  quality <- parts[5]
  starttime <- psdList[[1]]$starttime
  endtime <- psdList[[psdCount]]$endtime
  minfreq <- min(psdList[[1]]$freq)
  maxfreq <- max(psdList[[1]]$freq)
  nfreq <- length(psdList[[1]]$freq)
  units <- "acc"
  
  # Get instrument response
  evalresp <- getEvalresp(iris,network,station,location,channel,time=starttime,
                          minfreq=minfreq,maxfreq=maxfreq,nfreq=nfreq,units=units)

  # NOTE:  No try block needed as getEvalresp will generate appropriate errors
  
  # NOTE:  Because we're operating in dB space we need to think in terms of logarithms.
  
  # Instrument response converted to dB and then squared is:
  correction <- 10*log10(evalresp$amp) * 2
  
  # TODO:  This is the place we could split the list up based on snclq if we are
  # TODO:  doing a comparison between different snclq's
  
  # Need to tranpose here to have rows=psds and columns=freqs
  rawNoiseMatrix <- t(sapply(psdList, getElement, "spec"))
  
  # Sanity checks
  if ( nrow(rawNoiseMatrix) != psdCount ) {
    stop(paste("psdList2NoiseMatrix: psdCount =", psdCount,
               "and nrow(rawNoiseMatrix) =", nrow(rawNoiseMatrix),
               "are not equal."))
  }
  if ( ncol(rawNoiseMatrix) != length(evalresp$freq) ) {
    stop(paste("psdList2NoiseMatrix: length(evalresp$freq) =", length(evalresp$freq),
               "and ncol(rawNoiseMatrix) =", ncol(rawNoiseMatrix),
               "are not equal."))
  }
  
  
  # Create a correciton matrix that matches the dimensions of noiseMatrix
  correctionMatrix <- matrix(rep(correction,times=psdCount), nrow=psdCount, byrow=TRUE)
  
  # Dividing by the correction, in dB space, is just subtracting
  noiseMatrix <- rawNoiseMatrix - correctionMatrix
  
  return( noiseMatrix )
}


# psdDF2NoiseMatrix ------------------------------------------------------------

# Create instrument corrected noiseMatrix from psdDF dataframe

# psdDF is obtained from getPSDMeasurements (currently) in the seismicMetrics package
#
# Example rows from psdDF
#
#             target           starttime             endtime  frequency amplitude phase
# 1 IU.ANMO.00.BHZ.M 2013-05-01 00:30:00 2013-05-01 01:30:00 0.00532474   26.6653     0
# 2 IU.ANMO.00.BHZ.M 2013-05-01 00:30:00 2013-05-01 01:30:00 0.00580668   26.6653     0
# 3 IU.ANMO.00.BHZ.M 2013-05-01 00:30:00 2013-05-01 01:30:00 0.00633222   25.9748     0
# 4 IU.ANMO.00.BHZ.M 2013-05-01 00:30:00 2013-05-01 01:30:00 0.00690534   25.9748     0
# 5 IU.ANMO.00.BHZ.M 2013-05-01 00:30:00 2013-05-01 01:30:00 0.00753033   23.0990     0
# ... ... through all frequencies
# ... repeated for each PSD, typically associated with an hour long chunk of signal

psdDF2NoiseMatrix <- function(DF) {
  
  # Need to ensure that the IrisClient object exists as R has a default
  # "iris" dataset.  See help(iris, package="datasets")
  if (class(iris) == "data.frame") {
    iris <- new("IrisClient")
  }
  
  # TODO:  psdDF2NoiseMatrix should check to make sure that there is no
  # TODO:  end-of-epoch between the first PSD in the list and the last.
  # TODO:  This will only be important if there are significant changes in
  # TODO:  the instrument response.
  
  targetCount <- length(unique(DF$target))
  if (targetCount > 1) {
    stop(paste("psdLDF2NoiseMatrix: more than one target in PSD dataframe"))
  }
  
  psdCount <- length(unique(DF$starttime))
  
  # Get information needed for getEvalresp
  snclq <- DF$target[1]
  parts <- unlist(stringr::str_split(snclq,"[.]"))
  network <- parts[1]
  station <- parts[2]
  location <- parts[3]
  channel <- parts[4]
  quality <- parts[5]
  starttime <- min(DF$starttime)
  endtime <- max(DF$endtime)
  
  minfreq <- min(DF$frequency)
  maxfreq <- max(DF$frequency)
  nfreq <- length(unique(DF$frequency))
  units <- "acc"
  
  evalresp <- getEvalresp(iris,network,station,location,channel,time=starttime,
                          minfreq=minfreq,maxfreq=maxfreq,nfreq=nfreq,units=units)
  
  # NOTE:  No try block needed as getEvalresp will generate appropriate errors
  
  # NOTE:  Because we're operating in dB space we need to think in terms of logarithms.
  
  # Instrument response converted to dB and then squared is:
  correction <- 10*log10(evalresp$amp) * 2
  
  # Reshape DF amplitudes into a matrix with correct ordering of frequencies and correction applied
  
  noiseMatrix <- matrix(nrow=psdCount, ncol=nfreq)
  
  starttimes <- sort(unique(DF$starttime))
  
  for (i in seq(psdCount)) {    
    mask <- DF$starttime == starttimes[i]
    psd <- DF[mask,]
    psd <- psd[order(psd$frequency),]    
    noiseMatrix[i,] <- psd$amplitude - correction
  }
  
  return(noiseMatrix)
  
}


# noiseMatrix2PdfMatrix ---------------------------------------------------

# Create PDF matrix from instrument corrected noiseMatrix

#  INPUT:  matrix where columns are frequencies and rows are individual, corrected PSDs
# OUTPUT:  matrix where columns are frequencies and rows are counts of how many input PSDs have that power level

noiseMatrix2PdfMatrix <- function(noiseMatrix, lo=-200, hi=-50, binSize=1) {
   
  # NOTE:  Define a function to convert one column of noiseMatrix into a 
  # NOTE:  column of histogram counts.  However many rows exist in noiseMatrix,
  # NOTE:  pdfMatrix rows will always represt seq(lo,hi,binSize).
  
  histCounts <- function(x, lo, hi, binSize) {
    breaks <- seq(lo,hi,binSize)
    nbins <- length(breaks) - 1
    discretizedValue <- .bincode(x, breaks)
    return(tabulate(discretizedValue, nbins))  
  }
  
  pdfMatrix <- apply(noiseMatrix, 2, histCounts, lo, hi, binSize)
  
  return(pdfMatrix)
  
}


# noiseModels -----------------------------------------------------------

# Generate NHNM and NHLM

noiseModels <- function(freq) {
  
  if (missing(freq)) {
    stop("noiseModels: argument 'freq' is missing")  
  }
  
  period <- 1/freq
  
  # NOTE:  Original New High/Low Noise Models in Peterson paper:
  # NOTE:    http://www.mttmllr.com/ADS/DATA/peterson_usgs_seismic_noise_ofr93-322.pdf
  # NOTE:
  # NOTE:  IRIS DMC 2005 paper
  # NOTE:    http://www.earth.northwestern.edu/people/seth/327/HV/McNamaraetal_AmbientNoise.3.0.pdf
  # NOTE:
  # NOTE:  "Ambient Noise Levels in the Continental US", 2003 
  # NOTE:     http://geohazards.usgs.gov/staffweb/mcnamara/PDFweb/McNamaraBuland.pdf
  # NOTE:
  # NOTE:  Source code to compare:
  # NOTE:    https://github.com/g2e/seizmo/blob/master/noise/nlnm.m
  
  # New Low Noise Model in acceleration: NLNMacc = A + B*log10(T) referred to 1 (m/s^2)^2/Hz
  # T ( minimum period)
  
  # Create the NLNM ----------------------
  
  NLNM_table <- utils::read.table(header=TRUE, text="
minPeriod       A       B
     0.10 -162.36    5.64
     0.17 -166.70    0.00
     0.40 -170.00   -8.30
     0.80 -166.40   28.90
     1.24 -168.60   52.48
     2.40 -159.98   29.81
     4.30 -141.10    0.00
     5.00  -71.36  -99.77
     6.00  -97.26  -66.49
    10.00 -132.18  -31.57
    12.00 -205.27   36.16
    15.60  -37.65 -104.33
    21.90 -114.37  -47.10
    31.60 -160.58  -16.28
    45.00 -187.50    0.00
    70.00 -216.47   15.70
   101.00 -185.00    0.00
   154.00 -168.34   -7.61
   328.00 -217.43   11.90
   600.00 -258.28   26.60
   100000 -346.88   48.75
")

  # We create breaks based on the first column of each table to figure out
  # the appropriate A and B for each period.
  NLNM_breaks <- c(NLNM_table$minPeriod)
  NLNM_rows <- .bincode(period,NLNM_breaks,right=FALSE,include.lowest=TRUE)
  
  # Here is the vectorized version of "A + B*log10(period)":
  NLNM <- NLNM_table$A[NLNM_rows] + NLNM_table$B[NLNM_rows]*log10(period)
  
  # Create the NHNM ----------------------
  
  NHNM_table <- utils::read.table(header=TRUE, text="
 minPeriod       A       B
      0.10 -108.73  -17.23
      0.22 -150.34  -80.50
      0.32 -122.31  -23.87
      0.80 -116.85   32.51
      3.80 -108.48   18.08
      4.60  -74.66  -32.95
      6.30    0.66 -127.18 
      7.90  -93.37  -22.42
     15.40   73.54 -162.98
     20.00 -151.52   10.01
    100000 -206.66   31.63       
")

  # We create breaks based on the first column of each table to figure out
  # the appropriate A and B for each period.
  NHNM_breaks <- c(NHNM_table$minPeriod)
  NHNM_rows <- .bincode(period,NHNM_breaks,right=FALSE,include.lowest=TRUE)
  
  # Here is the vectorized version of "A + B*log10(period)":
  NHNM <- NHNM_table$A[NHNM_rows] + NHNM_table$B[NHNM_rows]*log10(period)
  
  return(data.frame(nlnm=NLNM,nhnm=NHNM))
  
}


# psdStatistics -------------------------------------------------------

# Calculate basic statistics on all PSDs in list or dataframe.

psdStatistics <- function(PSDs) {
  
  # Get list of frequencies and a noiseMatrix
  
  if (class(PSDs) == "list") {   
    
    freq <- PSDs[[1]]$freq
    noiseMatrix <- psdList2NoiseMatrix(PSDs) 
    
  } else if (class(PSDs) == "data.frame") {  
    
    freq <- sort(unique(PSDs$freq))
    noiseMatrix <- psdDF2NoiseMatrix(PSDs)  
    
  }
    
  nrow <- nrow(noiseMatrix)
  ncol <- ncol(noiseMatrix)
  
  # NOTE:  The noiseMatrix has one column per frequency and one row per PSD.
  # NOTE:  We will create vectors to store per frequency values and then fill 
  # NOTE:  them in one frequency at a time.
  
  colMax <- rep(NA,ncol)
  colMin <- rep(NA,ncol)
  colMean <- rep(NA,ncol)
  colMedian <- rep(NA,ncol)
  
  for (j in seq(ncol)) {
    colMax[j] <- max(noiseMatrix[,j])
    colMin[j] <- min(noiseMatrix[,j])
    colMean[j] <- mean(noiseMatrix[,j])
    colMedian[j] <- median(noiseMatrix[,j])
  }
  
  # Calculate pctAboveNHNM and pctBelowNLNM ------------
  
  noiseModels <- noiseModels(freq)
  
  aboveCount <- rep(0,ncol)
  belowCount <- rep(0,ncol)
  
  for (j in seq(ncol)) {    
    aboveCount[j] <- length(which(noiseMatrix[,j] > noiseModels$nhnm[j]))
    belowCount[j] <- length(which(noiseMatrix[,j] < noiseModels$nlnm[j]))
  }
  
  pct_above <- 100 * aboveCount/nrow
  pct_below <- 100 * belowCount/nrow
  
  
  # Now calculate mode -------------------------------------
  
  # For mode, we need to convert to the discretized pdfMatrix
  # that contains in each bin, the count of PSDs that have 
  # that value.  The mode is just the bin with the highest count.
  
  # Default dbBins = seq(-200,-50,1)
  pdfBins <- seq(-200,-50,1)
  pdfMatrix <- noiseMatrix2PdfMatrix(noiseMatrix)
  
  colMode <- rep(NA,ncol)
  for (j in seq(ncol)) {
    modeIndex <- which(pdfMatrix[,j] == max(pdfMatrix[,j]))[1]
    colMode[j] <- pdfBins[modeIndex]
  }    
  
  return(list(noiseMatrix=noiseMatrix,
              pdfMatrix=pdfMatrix,
              freq=freq,
              pdfBins=pdfBins,
              max=colMax,
              min=colMin,
              mean=colMean,
              median=colMedian,
              mode=colMode,
              nlnm=noiseModels$nlnm,
              nhnm=noiseModels$nhnm,
              pct_above=pct_above,
              pct_below=pct_below))  
}


# psdPlot --------------------------------------------------------------

# Plot instrument corrected noise values of PSDs in list or dataframe

psdPlot <- function(PSDs, style='psd', showNoiseModel=TRUE, showMaxMin=TRUE, showMode=TRUE, showMean=FALSE, ...) {
  
  if (class(PSDs) == "list") {
    
    psdCount <- length(PSDs)
    
    # Get information from the PSDs
    # NOTE:  Each element of psdList is a list
    snclq <- as.character(PSDs[[1]]$snclq)
    parts <- unlist(stringr::str_split(snclq,"[.]"))
    network <- parts[1]
    station <- parts[2]
    location <- parts[3]
    channel <- parts[4]
    quality <- parts[5]
    starttime <- PSDs[[1]]$starttime
    endtime <- PSDs[[psdCount]]$endtime
    freq <- PSDs[[1]]$freq
        
  } else if (class(PSDs) == "data.frame") {
    
    psdCount <- length(unique(PSDs$starttime))

    # Get information needed from the psdDF
    snclq <- PSDs$target[1]
    parts <- unlist(stringr::str_split(snclq,"[.]"))
    network <- parts[1]
    station <- parts[2]
    location <- parts[3]
    channel <- parts[4]
    quality <- parts[5]
    starttime <- min(PSDs$starttime)
    endtime <- max(PSDs$endtime)    
    freq <- sort(unique(PSDs$freq))
    
  } else {
    
    stop(paste("psdPlot: PSDs is of class '",class(PSDs),"' instead of 'list' or 'data.frame'.",sep=""))
    
  }
  
  # Generate basic statics as well as noiseMatrix and pdfMatrix
  
  stats <- psdStatistics(PSDs)
  
  noiseMatrix <- stats$noiseMatrix
  pdfMatrix <- stats$pdfMatrix
    
  # Choose appropriate limits for period
  period <- 1/freq
  if (stringr::str_detect(channel,"^B")) {
    xlim <- c(0.1,100)
    verticalLines <- c(seq(.1,1,length.out=10),
                       seq(1,10,length.out=10),
                       seq(10,100,length.out=10))
  } else {
    xlim <- c(min(period),max(period))
    verticalLines <- c(seq(.1,1,length.out=10),
                       seq(1,10,length.out=10),
                       seq(10,100,length.out=10))
  }
  
  # Choose appropriate limits for dB
  ylo <- -200
  yhi <- -50
  ylim <- c(ylo,yhi)
  horizontalLines <- seq(ylo,yhi,10)
  
  if (style == 'psd') {
    
    # Plot the first PSD
    plot(period, noiseMatrix[1,], log="x",
         xlim=xlim, ylim=ylim,
         xlab="Period (Sec)", ylab="Power (dB)",
         las=1,
         main=paste(psdCount,"corrected, hourly PSDs for",snclq),
         ...)
    
    # Add all additional PSDs
    for (i in seq(nrow(noiseMatrix))) {
      graphics::points(period, noiseMatrix[i,], ...)
    }
    
  } else if (style == 'pdf') {
    
    # Set up colors and breaks
    cols <- c('grey90', rev(grDevices::rainbow(30))[4:30])
    breaks <- c(0,seq(0.001,max(pdfMatrix),length.out=length(cols)))
    
    # NOTE:  To get image() to plot with the same axes as the data displayed as a table you
    # NOTE:  have to transpose and reverse the order of the columns.  This means we also have
    # NOTE:  to reverse the order of the X-axis locations associated with the columns.
    
    # Initial plot
    graphics::image(t(pdfMatrix)[ncol(pdfMatrix):1,],
                    x=rev(period),
                    y=seq(ylo,yhi),
                    breaks=breaks,
                    col=cols,
                    las=1, log="x",
                    xlab="Period (Sec)", ylab="Power dB",
                    main=paste("PDF plot of",psdCount,"corrected, hourly PSDs for",snclq),
                    ...)      
    
  } else {
    
    stop(paste("psdPlot: style '",style,"' is not recognized.",sep=""))
    
  }

  
  # Add various lines from the statistics
  legend <- c()
  lwd <- c()
  col <- c()
  if (showNoiseModel) {
    graphics::points(period, stats$nlnm, type='l', col='gray50', lwd=2)
    graphics::points(period, stats$nhnm, type='l', col='gray50', lwd=2)
    legend <- c(legend,"NLNM, NHNM")
    lwd <- c(lwd,2)
    col <- c(col,'gray50')
  }
  if (showMaxMin) {
    graphics::points(period, stats$max, type='l', col='blue', lwd=2)
    graphics::points(period, stats$min, type='l', col='red', lwd=2)
    legend <- c(legend,"max","min")
    lwd <- c(lwd,2,2)
    col <- c(col,'blue','red')
  }
  if (showMode) {
    graphics::points(period, stats$mode, type='l', col='yellow', lwd=3)
    legend <- c(legend,"mode")
    lwd <- c(lwd,3)
    col <- c(col,'yellow')
  }
  if (showMean) {
    graphics::points(period, stats$mean, type='l', col='orange', lwd=4)
    legend <- c(legend,"mean")
    lwd <- c(lwd,4)
    col <- c(col,'orange')    
  }
  
  # Add grid lines
  graphics::abline(h=horizontalLines, v=verticalLines, col='gray50', lty='dotted')
  
  # Add a legend
  legend("topleft", bg='white',
         legend=legend,
         lwd=lwd,
         col=col)
  
  # Add the starttime and endtime
  legend("topright", bg='white',
         legend=c(paste(starttime,"  start "),
                  paste(endtime,"  end ")))
  
}


################################################################################
# Hilbert method from the seewave package
#   cran.r-project.org/web/packages/seewave
#
# This function is called to calculate hilbert and envelope transforms.
################################################################################

hilbertFFT <- function(x) {
  # Data should already be demeaned, detrended and tapered
  n <- length(x)
  ff <- stats::fft(x)
  h <- rep(0,n)
  
  if (n>0 & 2*floor(n/2)==n) {
    h[c(1, n/2+1)] <- 1
    h[2:n/2] <- 2
  } else {
    if (n>0) {
      h[1] <- 1
      h[2:(n+1)/2] <- 2
    }
  }
  
  hfft <- stats::fft(ff*h, inverse=TRUE)/length(ff)
  
  return(hfft)
}

################################################################################
# END
################################################################################
