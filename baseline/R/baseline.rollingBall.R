## $Id: baseline.rollingBall.R 170 2011-01-03 20:38:25Z bhm $
baseline.rollingBall <- function(spectra, wm, ws){
## Ideas from Rolling Ball algorithm for X-ray spectra by M.A.Kneen
## and H.J. Annegarn. Variable window width has been left out
#
# INPUT:
# spectra - rows of spectra
# wm      - width of local window for minimization/maximization
# ws      - width of local window for smoothing
#
# OUTPUT:
# baseline  - proposed baseline
# corrected - baseline corrected spectra

  np <- dim(spectra)

  # Initialization
  y        <- numeric(np[2]) # Spectrum
  T1       <- numeric(np[2]) # Minimizers
  T2       <- numeric(np[2]) # Maximizers
  basel    <- numeric(np[2]) # Smoothers
  baseline <- matrix(0,np[1],np[2]) # Final baseline
  u1       <- 0 # Place holder #1
  u2       <- 0 # Place holder #2
  v        <- 0 # Sum holder

  # Loop through spectra
  for(a in 1:np[1]){
    y <- spectra[a,]
    # Minimize
	u1 <- ceiling((wm+1)/2)+1
	T1[1] <- min(y[1:u1])
    for(i in 2:wm){                  # -- Start of spectrum --
	  u2 <- u1+1+(i%%2)
	  T1[i] <- min(y[(u1+1):(u2)],T1[i-1]) # Check if new is smaller
	  u1 <- u2
    }
    for(i in (wm+1):(np[2]-wm)){     # -- Main part of spectrum --
	  if((y[u1+1]<=T1[i-1]) && (y[u1-wm]!=T1[i-1]))
		T1[i] <- y[u1+1]   # Next is smaller
      else
        T1[i] <- min(y[(i-wm):(i+wm)])
	  u1 <- u1+1
    }
	u1 <- np[2]-2*wm-1
    for(i in (np[2]-wm+1):np[2]){    # -- End of spectrum --
      u2<- u1+1+((i+1)%%2)
	  if(min(y[u1:(u2-1)])>T1[i-1])
	    T1[i] <- T1[i-1]   # Removed is larger
	  else
        T1[i] <- min(y[(u2):np[2]])
	  u1 <- u2
    }

    # Maximize
	u1 <- ceiling((wm+1)/2)+1
	T2[1] <- max(T1[1:u1])
    for(i in 2:wm){                  # -- Start of spectrum --
	  u2 <- u1+1+(i%%2)
	  T2[i] <- max(T1[(u1+1):(u2)],T2[i-1]) # Check if new is larger
	  u1 <- u2
    }
    for(i in (wm+1):(np[2]-wm)){     # -- Main part of spectrum --
	  if((T1[u1+1] >= T2[i-1]) && (T1[u1-wm]!=T2[i-1]))
	    T2[i] <- T1[u1+1] # Next is larger
      else
	    T2[i] <- max(T1[(i-wm):(i+wm)])
	  u1 <- u1+1
    }
	u1 <- np[2]-2*wm-1
    for(i in (np[2]-wm+1):np[2]){    # -- End of spectrum --
      u2<- u1+1+((i+1)%%2)
	  if(max(T1[u1:(u2-1)])<T2[i-1])
	    T2[i] <- T2[i-1]   # Removed is smaller
	  else
        T2[i] <- max(T1[(u2):np[2]])
	  u1 <- u2
    }

    # Smooth
	u1 <- ceiling(ws/2)
	v <- sum(T2[1:u1])
    for(i in 1:ws){                  # -- Start of spectrum --
	  u2<-u1+1+(i%%2)
	  v <- v+sum(T2[(u1+1):u2])      # Sum so far
	  basel[i] <- v/u2               # Mean so far
	  u1 <- u2
    }
	v <- sum(T2[1:(2*ws+1)])
	basel[ws+1] <- v/(2*ws+1)
    for(i in (ws+2):(np[2]-ws)){     # -- Main part of spectrum --
	  v <- v - T2[i-ws-1] + T2[i+ws] # Sum so far
	  basel[i] <- v/(2*ws+1)         # Mean so far
    }
	u1 <- np[2]-2*ws+1
	v <- v-T2[u1]                    # Sum so far
	basel[np[2]-ws+1] <- v/(2*ws)    # Mean so far
    for(i in (np[2]-ws+2):np[2]){    # -- End of spectrum --
	  u2 <- u1+1+(i+1)%%2
	  v <- v-sum(T2[(u1):(u2-1)])    # Sum so far
	  basel[i] <- v/(np[2]-u2+1)     # Mean so far
	  u1 <- u2
    }
	baseline[a,] <- basel
  }

  corrected <- spectra-baseline
  list(baseline=baseline, corrected=corrected)
}
