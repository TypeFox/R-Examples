baseline.medianWindow <- function(spectra, hwm, hws, end=FALSE){
## Implementation and extention of Marks S. Friedrichs model-free algorithm
## Coded by Kristian Hovde Liland and Bjørn-Helge Mevik
## $Id: baseline.medianWindow.R 170 2011-01-03 20:38:25Z bhm $
#
# INPUT:
# spectra - rows of spectra
# hwm     - half width of median window
# hws     - half width of smoothing window
# end     - original endrule (boolean)
#
# OUTPUT:
# baseline  - proposed baseline
# corrected - baseline corrected spectra

  if (missing(hws)) hws <- hwm
  if(end)
	endrule <- "constant"
  else
    endrule <- "median"
  np <- dim(spectra)
  spect    <- matrix(0,np[1],np[2])
  baseline <- matrix(0,np[1],np[2])

  # Keeps first and last values of original spectrum
  spect[,1]     <- spectra[,1]
  spect[,np[2]] <- spectra[,np[2]]

  ## FIXME: Idea: subsampling followed by interpolation?

  # Median of windows for each variable
  ## FIXME: This gives different values at the ends.  That can be fixed by
  ## using the old implementation at the ends.
  for(j in 1:np[1]){
      spect[j,] <- runmed(spectra[j,], k = 2*hwm + 1, endrule = endrule)
  }

  # Smoothing by gaussian weighting
  gau <- dnorm(-hws:hws,0,hws/2)
  for(i in 1:np[2]){
    cutL <- i-(hws+1)
    if(cutL < 0) # Low cut-off
      g <- gau[(-cutL+1):length(gau)]
    else if(np[2] < (i+hws)) # High cut-off
        g <- gau[1:(length(gau)-((i+hws)-np[2]))]
      else
        g <- gau
    g <- g/sum(g)
    baseline[,i] <- spect[,max(1,i-hws):min(i+hws,np[2]),drop=FALSE]%*%g
  }

  corrected <- spectra-baseline
  list(baseline=baseline, corrected=corrected)
}
