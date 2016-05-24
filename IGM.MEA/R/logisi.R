## Implement logisi method for burst detection.
## Author: Zhengzheng Zhang
## Copyright: GPL.

.logisi.find.burst <- function(spikes, debug=FALSE) {

  ## For one spike train, find the burst using log isi method.
  ## e.g.
  ## .find.bursts(s$spikes[[5]])
  ## init.
  ## params currently in LOGISI.PAR
  ##

  no.bursts = NA;                       #value to return if no bursts found.

  logisi.par <- list(min.ibi=0.800,   min.durn=0.05, min.spikes=6,
                     isi.low=0.02)
  par = logisi.par
  ##beg.isi =    par$beg.isi
  ##end.isi =    par$end.isi
  min.ibi =      par$min.ibi
  min.durn =     par$min.durn
  min.spikes =   par$min.spikes
  isi.low =      par$isi.low
  
  nspikes = length(spikes)

  ## Create a temp array for the storage of the bursts.  Assume that
  ## it will not be longer than Nspikes/2 since we need at least two
  ## spikes to be in a burst.
  
  max.bursts <- floor(nspikes/2)
  bursts <- matrix(NA, nrow=max.bursts, ncol=3)
  colnames(bursts) = c("beg", "end", "IBI")
  burst <- 0                            #current burst number

  ## Phase 1 -- burst detection. Each interspike interval of the data 
  ## is compared with the threshold THRE. If the interval is greater 
  ## than the threshold value, it can not be part of a burst; if the 
  ## interval is smaller or equal to the threhold, the interval may be 
  ## part of a burst.
 


  ## LAST.END is the time of the last spike in the previous burst.
  ## This is used to calculate the IBI.
  ## For the first burst, this is no previous IBI
  last.end = NA;                        #for first burst, there is no IBI.

  n = 2
  in.burst = FALSE
  
  while ( n < nspikes) {
    
    next.isi = spikes[n] - spikes[n-1]
    if (in.burst) {
      if (next.isi > isi.low) {
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
      if (next.isi <= isi.low) {
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
  
  if (debug) {
    print("End of phase1\n")
    print(bursts)
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

  if (debug) {
    print("End of phase 2\n")
    print(bursts)
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


## Peak finding algorithm; taken from R-help mailing list.
## Author: Brian Ripley.
.locpeaks <- function (series, span = 3)
{
    z <- embed(series, span)
    s <- span%/%2
    v <- max.col(z) == 1 + s
    result <- c(rep(FALSE, s), v)
    result <- result[1:(length(result) - s)]
    which(result)
}

.logisi.compute <- function(s, min.nspikes = 10,
                           breaks.max = 100,
                           channel = ncells+1,
                           span = 1, span.max = 50, Rat = 0.08, plot = FALSE)
{

  ## N --> MIN.NSPIKES == minimum number of spikes in a channel.
  ## br.max --> breaks.max
  ## channel.
  ## 
  ## Given the spike data structure S,
  ## compute the log ISI transform and return useful information.

  ## This function should be expanded to find the peaks in the
  ## histogram either of each channel or of the grand average.
  
  h = list()                                                 # hist objects  
  total.isi = NULL
  ncells <- s$NCells
  if (plot){
      par(mfrow=c(8,8), mar=c(3,3,1,1), ask=FALSE, oma=c(0,1.5,0,0))
  }
  
  for (i in 1:ncells){
    if (s$nspikes[i] >= min.nspikes ) {
      ## Channel has "enough" spikes for method.
      isi = diff( s$spikes[[i]] )
      total.isi =c(total.isi, isi) 
      B = sqrt(s$nspikes[i])
      if (B > breaks.max){
          B = breaks.max
      }
      if (plot){
          title = sprintf("%d # %d", i, s$nspikes[i])
          h[[i]] = hist(log(isi), br = B, main = title, xlab = "logisi")
          abline(h = mean(h[[i]]$counts), col = 'blue')
      } else{
             h[[i]] = hist(log(isi), br = B, plot = FALSE)
      }
    } else {
      ## insufficient spikes to compute ISI histogram.
      if (plot){
          title = sprintf("%d # %d", i, s$nspikes[i])
          plot(1, type='n', main=title)
      }
    }
  }

  ## The log histogram for the grand average.
  B = sqrt(length(total.isi))
  if (B > breaks.max){
    B = breaks.max
  }
  file = s$file
  if (plot){
      last.h = hist(log(total.isi), br = B, main = "All", xlab = "logisi")                             
      abline(h = mean(last.h$counts), col = 'red')
      mtext(file, side=2, outer=TRUE)
  }else{
        last.h = hist(log(total.isi), br = B, plot = FALSE) 
  }
  h[[ncells+1]] = last.h
  counts = h[[channel]]$counts
  breaks = h[[channel]]$breaks

  if (plot){
      ## Find the peaks in the histogram 
      ## either of each channel or of the grand average.
      par(mfrow=c(1,1))
      plot(counts, type='l', xlab="Logisi (ms)", ylab="Frequency", xaxt="n")
      axis(1, 0:length(counts), format(exp(breaks), sci=TRUE))   
  }   
  peaks = .locpeaks(counts, span)
 
  ## Return some dummy values; these might be times of thresholds.
  ## e.g. max1 might be the time of the first peak interval.
  ## res <- list(max1=.1, max2=.4)
  MAX = max(counts)
  if (length(peaks)==0){
      peak.max = -Inf
  }else{
        peak.max = max(counts[peaks])
  }     

  if (span.max >= length(counts)){
      span.max = length(counts) -1 
  }
    
  ## Find the no. of peaks no more than 6, and
  ## the golobal max is one of peaks.
  while(length(peaks) >6 || MAX!=peak.max){
        span= span + 1
        peaks = .locpeaks(counts, span)
        if (length(peaks)==0){
            peak.max = -Inf
        }else{
              peak.max = max(counts[peaks])
        }    
        if (span > span.max){
            peaks = 0
            break
        }
  }


  if (length(peaks)!=1 || (length(peaks)==1&& peaks!=0)){
      if (plot){  
          points( peaks, counts[peaks], pch=19, col='blue')
      }
  }else{
        browser()
  }     

  ## Find the local minimums between two successive peaks, and report the lowest.
  ## If the peak finding algorithm gives some unlikely peaks between them,
  ## then the peaks will be filtered out.
  ## Rat = 0.08        # a threhold for filtering unreasonable peaks

  pos = -1             # flag
  len = length(peaks)  # initial length
  j = 1
  mini = NULL
  R= NULL

  while (pos==-1 || j < len){
         len = length(peaks)
         if (len >= 2){
             loc.min = min(counts[peaks[j]:peaks[j+1]])
             temp = c(peaks[j]:peaks[j+1])
             pos = temp[counts[peaks[j]:peaks[j+1]]==loc.min]
             pos = pos[length(pos)]           # last local min
             pair = c(counts[peaks[j]], counts[peaks[j+1]])
             smallest = c(j,j+1)[which.min(pair)]  
             Diff = counts[peaks[smallest]] - counts[pos]
             ## If the second peaks occurs after the first in the next 3 
             ## breaks, then remove the smallest peak.
             if (diff(peaks[j:(j+1)])<=3){
                 peaks = peaks[-smallest]
                 pos = -1
                 j=1
             }else{ 
                   if (Diff==0){
                       peaks = peaks[-smallest]
                       pos = -1
                       j=1
                   }else{
                         ## define a ratio
                         ratio = Diff/(max(counts[peaks]) - counts[pos])
                         ## If the ratio is less than Rat, remove the smallest peak.
                         if (ratio < Rat){
                             peaks = peaks[-smallest]
                             pos = -1
                             j=1
                         }else{
                               if (ratio < 0 || ratio > 1){
                               browser()
                               }
                               mini = c(mini, pos)
                               R = c(R, ratio)
                               j = j+1
                         }
                   }  
             } 

         }else{
               lowest = -2
               break
         }  
  }
  if (length(mini)!=0){
      M = min(counts[mini])
      lowest = mini[counts[mini] == M][1]           # choose the first
  }else{
        lowest = -2
  }

  if (lowest != -2){
      if (plot){
                points( lowest, counts[lowest], pch=19, col='red')
      }
      a1 = h[[channel]]$breaks[lowest]
      a2 = h[[channel]]$breaks[lowest+1]
      av.a = (a1+a2)/2
      loc.min = exp(av.a)
  }else{
      loc.min = NA                     
  }  



  b1 = h[[channel]]$breaks[peaks]
  b2 = h[[channel]]$breaks[peaks+1]
  av.b = (b1+b2)/2
  isi.peaks = exp(av.b)                             # time in seconds   
  res = list(max1=isi.peaks[1], max2=isi.peaks[2], max3=isi.peaks[3])
       
  return(list(Max=res,Locmin=loc.min))

}
