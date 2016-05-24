## Jennifer Simonotto
## 2011-06-07

.episode.detection <- function(s, max.burst.length=40) {
  ## max.burst.length is the longest allowed burst.
  ## Any burst longer than this is ignored. Set this value to "Inf" if
  ## you do not want any filtering of burst times.

  ## Takes in a "s" data structure

  ## Needed in global environment?
  ## mi.par <- list(beg.isi =    0.1,
  ##                 end.isi =    0.25,
  ##                 min.ibi =    0.8,
  ##                 min.durn =   0.05,
  ##                 min.spikes = 6)

  mi_array<-sapply(s$spikes,mi.find.bursts,s$parameters$mi.par)

  ## do episode detection
  ## First, determine the number of active channels that have spike bursts 
  
  channelslist = 1:length(s$channels)
  burst_channels = channelslist[which(unlist(lapply(mi_array,is.matrix))*1>0)]
  
  ## Initialize counters and variables to 0 and initialize a vector
  ## for the number of bursts within active channels
	
  jj = 0
  jjjjj = 0
  spikebursts_start = 0
  spikebursts_end = 0
  duration_checker = 0
  numbursts=rep(0,length(burst_channels))

  ## Extract the number of bursts per active channel as well as start,
  ## stop and duration of each burst to its own vector

  for (j in burst_channels){
    jj = jj+1;
    numbursts[jj] = length(mi_array[[j]][,1])
    spikebursts_start = append(spikebursts_start,s$spikes[[j]][mi_array[[j]][,1]])
    spikebursts_end = append(spikebursts_end,s$spikes[[j]][mi_array[[j]][,2]])
    duration_checker = append(duration_checker,mi_array[[j]][,5])
    
  }

  ## As we have initialized all these vectors to 0, remove this first
  ## zero from each vector.  Also check duration, and remove any
  ## bursts which are longer than 40 seconds

  spikebursts_start = spikebursts_start[2:length(spikebursts_start)]
  spikebursts_end = spikebursts_end[2:length(spikebursts_end)]
  duration_checker = duration_checker[2:length(duration_checker)]
  duration_checked = which(duration_checker<max.burst.length)
  
  spikebursts_start_longremoved = spikebursts_start[duration_checked]
  spikebursts_end_longremoved = spikebursts_end[duration_checked]
  
  ## Sort the start times of the cleaned up vectors, extracting also
  ## the sorted start index and applying it to the end times vector.
	
  llama = sort(spikebursts_start_longremoved,index.return=TRUE)
  spikes_endsorted = spikebursts_end_longremoved[llama$ix]
  
  ## Initialize index vectors, with size of the length of the cleaned up vectors
	
  jj=0
  proper_indexes1 = rep(0,length(spikebursts_start_longremoved))
  proper_indexes = rep(0,length(duration_checked))
  
  ## Create an index that indicates the channel number for the number
  ## of bursts for that channel
	
  for (j in 1:length(burst_channels)){
    for (k in 1:numbursts[j]){
      jj = jj+1
      proper_indexes1[jj] = burst_channels[j];
    }
  }
  
  ## Sort channel indices according to the start time

  N = length(duration_checked )
  proper_indexes = proper_indexes1[duration_checked]
  proper_indexes_sorted = proper_indexes[llama$ix]
  
  ## algorithm for checking if one spike burst overlaps the next
  ## Initialize things and create lists
  
  j = 1;
  tempjjj=0
  savedpath<-list()
  savedstart<-list()
  savedend<-list()
  
  ## while within the number of detected bursts for whole file, set
  ## path, start and end points to current burst
  
  while ( j < N ){
    path = 0
    pathIdx = 1;
    
    path[pathIdx] = proper_indexes_sorted[j];
    StrPnt          = llama$x[j];
    EndPnt          = spikes_endsorted[j];
    
    ## Default is that this is not an episode, and so ActivePath is
    ## false Look forward in sorted index files to see if the start
    ## time of the next burst happens before the present end point If
    ## it does, then we have an episode, and ActivePath is set to
    ## true, the path index is updated, as is the endpoint otherwise
    ## we look at the next burst, and so on.
		
    j2 = j + 1;
    
    ActivePath = FALSE;
    while( j2 < N ){
      if ( llama$x[j2] < EndPnt ){
        ActivePath = TRUE;
        pathIdx = pathIdx + 1; 
        path[pathIdx ] = proper_indexes_sorted[j2];  
        j = j2; which
        
        if (spikes_endsorted[j2] > EndPnt){
          EndPnt = spikes_endsorted[j2];
        }
        
      } else(break)
      j2 = j2 + 1;
    }
    
    ## If the total path length of burst that occurred concurrently is
    ## greater than one, save the path, start and end times of each
    ## episode

    if (length( path ) > 1) {
      jjjjj=jjjjj+1
      tempjjj=tempjjj+1
      savedpath[[jjjjj]] = path
      savedstart[[jjjjj]] = StrPnt 
      savedend[[jjjjj]]= EndPnt
    }
    j = j + 1;
  }
  
  
  ## return results, including episode start and end time, as well as
  ## list all channels who participated

  ## SJE: make sure results are vectors wherever possible.
  beg <- unlist(savedstart)
  if (is.null(beg)) {
    ## if no episodes found.
    res <- NULL
  } else {
    res <- list(beg=beg,
                end=unlist(savedend),
                savedpath=savedpath)
  }
  
  res
  
}
