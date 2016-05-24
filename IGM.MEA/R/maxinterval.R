## maxinterval.R --- maxinterval burst detection (from Neuroexplorer).
## Author: Stephen Eglen
## Copyright: GPL
## Fri 23 Feb 2007


mi.find.bursts <- function(spikes,mi.par) {

  ## For one spike train, find the burst using max interval method.
  ## e.g.
  ## .find.bursts(s$spikes[[5]])
  ## init.
  ## params currently in MI.PAR
  ##

  ##no.bursts = NA;                       #value to return if no bursts found.
  no.bursts = matrix(nrow=0,ncol=1)     #emtpy value nrow()=length() = 0.

  par = mi.par
  beg.isi =    par$beg.isi
  end.isi =    par$end.isi
  min.ibi =    par$min.ibi
  min.durn =   par$min.durn
  min.spikes = par$min.spikes
  
  nspikes = length(spikes)

  ## Create a temp array for the storage of the bursts.  Assume that
  ## it will not be longer than Nspikes/2 since we need at least two
  ## spikes to be in a burst.
  
  max.bursts <- floor(nspikes/2)
  bursts <- matrix(NA, nrow=max.bursts, ncol=3)
  colnames(bursts) = c("beg", "end", "IBI")
  burst <- 0                            #current burst number

  ## Phase 1 -- burst detection.  Here a burst is defined as starting
  ## when two consecutive spikes have an ISI less than BEG.ISI apart.
  ## The end of the burst is given when two spikes have an ISI greater
  ## than END.ISI.
  
  ## Find ISIs closer than beg.isi, and end with end.isi.


  ## LAST.END is the time of the last spike in the previous burst.
  ## This is used to calculate the IBI.
  ## For the first burst, this is no previous IBI
  last.end = NA;                        #for first burst, there is no IBI.

  n = 2
  in.burst = FALSE
  
  while ( n <= nspikes) {
    
    next.isi = spikes[n] - spikes[n-1]
    if (in.burst) {
      if (next.isi > end.isi) {
        ## end of burst
        end = n-1; in.burst = FALSE

        
        ibi =  spikes[beg] - last.end; last.end = spikes[end]
        res = c(beg, end, ibi)
        burst = burst + 1
        if (burst > max.bursts) {
          print("too many bursts!!!")
          browser()
        }
        bursts[burst,] <- res
      }
    } else {
      ## not yet in burst.
      if (next.isi < beg.isi) {
        ## Found the start of a new burst.
        beg = n-1; in.burst = TRUE
      }
    }
    n = n+1
  }

  ## At the end of the burst, check if we were in a burst when the
  ## train finished.
  if (in.burst) {
    end = nspikes
    ibi =  spikes[beg] - last.end
    res = c(beg, end, ibi)
    burst = burst + 1
    if (burst > max.bursts) {
      print("too many bursts!!!")
      browser()
    }
    bursts[burst,] <- res
  }

  ## Check if any bursts were found.
  if (burst > 0 ) {
    ## truncate to right length, as bursts will typically be very long.
    bursts = bursts[1:burst,,drop=FALSE]
  } else {
    ## no bursts were found, so return an empty structure.
    return(no.bursts)
  }
  
  
  ## Phase 2 -- merging of bursts.  Here we see if any pair of bursts
  ## have an IBI less than MIN.IBI; if so, we then merge the bursts.
  ## We specifically need to check when say three bursts are merged
  ## into one.
  
  
  ibis = bursts[,"IBI"]
  merge.bursts = which(ibis < min.ibi)
  
  if (any(merge.bursts)) {
    ## Merge bursts efficiently.  Work backwards through the list, and
    ## then delete the merged lines afterwards.  This works when we
    ## have say 3+ consecutive bursts that merge into one.

    for (burst in rev(merge.bursts)) {
      bursts[burst-1, "end"] = bursts[burst, "end"]
      bursts[burst, "end"] = NA         #not needed, but helpful.
    }
    bursts = bursts[-merge.bursts,,drop=FALSE] #delete the unwanted info.
  }

  ## Phase 3 -- remove small bursts: less than min duration (MIN.DURN), or
  ## having too few spikes (less than MIN.SPIKES).
  ## In this phase we have the possibility of deleting all spikes.

  ## LEN = number of spikes in a burst.
  ## DURN = duration of burst.
  len = bursts[,"end"] - bursts[,"beg"] + 1
  durn = spikes[bursts[,"end"]] - spikes[bursts[,"beg"]]
  bursts = cbind(bursts, len, durn)

  rejects = which ( (durn < min.durn) | ( len < min.spikes) )
  
  if (any(rejects)) {
    bursts = bursts[-rejects,,drop=FALSE]
  }

  if (nrow(bursts) == 0) {
    ## All the bursts were removed during phase 3.
    bursts = no.bursts
  } else {
    ## Compute mean ISIS
    len = bursts[,"end"] - bursts[,"beg"] + 1
    durn = spikes[bursts[,"end"]] - spikes[bursts[,"beg"]]
    mean.isis = durn/(len-1)

    ## Recompute IBI (only needed if phase 3 deleted some cells).
    if (nrow(bursts)>1) {
      ibi2 = c(NA, .calc.ibi(spikes, bursts))
    } else {
      ibi2 = NA
    }
    bursts[,"IBI"] = ibi2
    
    SI = rep(1, length(mean.isis ))
    bursts = cbind(bursts, mean.isis, SI)
  }
  
  ## End -- return burst structure.
  bursts
  
}

