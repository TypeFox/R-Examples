## SJE version of finding episodes by overlapping bursts.
## 2011-06-09

.findEpisodes <- function(s, max.burst.length=40) {
  ## Find the episodes of overlapping bursting activity.
  ##
  ## MAX.BURST.LENGTH is the longest allowed burst.
  ## Any burst longer than this is ignored. Set this value to "Inf" if
  ## you do not want any filtering of burst times.


  ## First, perform burst analysis (could check if not already done.)
  if (is.null(s$allb)) {
    s$allb <- lapply(s$spikes, mi.find.bursts,s$parameters$mi.par)
  }    
  allb <- s$allb


  ## Find the number of bursts in each channel; if a channel has no bursts, just
  ## let it be of length 1.
  nbursts <- sapply(s$allb, function(m) {n <- nrow(m); ifelse(is.null(n), 1, n)})
  rep(names(nbursts), times=nbursts)
  burst.ids <- as.numeric(rep(names(nbursts), times=nbursts))

  ## Flatten the "allb" burst information into one big matrix, to
  ## create the BURSTS matrix.
  ## Each row of the burst matrix contains three columns (beg,end, id):
  ## BEG, END: Spike number of the first and last spike in the burst.
  ## ID: channel number from where this burst came from

  allb.flat <- do.call(rbind, s$allb)
  allb.flat2 <- cbind(allb.flat[,1:2], burst.ids)

  ## Convert the spike number into a spike time, by looking up times
  ## in the spike trains for each channel.

  all.ids <- mapply(.extract.burst.times.with.id,
                    s$spikes, allb, as.numeric(names(allb)))
  flat <- do.call(rbind, all.ids)

  ## remove any electrodes which had no bursts.
  zero.bursts <- which(is.na(flat[,1]))
  if (any(zero.bursts))
    flat <- flat[-zero.bursts,]
  colnames(flat) <- c("beg", "end", "id")
  
  ## Remove any bursts which are too long.
  long.bursts <- apply(flat, 1, function(x) { (x["end"] - x["beg"]) > max.burst.length})
  if (any(long.bursts)) {
    .printf("%d long bursts removed\n", sum(long.bursts))
    flat <- flat[!long.bursts,]
  }

  
  ## Now that we have the BURSTS information in the right format, we
  ## can now find the episodes of overlapping bursts.

  if (nrow(flat)==0 || is.null(flat)) {
    ## We have no bursts to make episodes from
    ## (e.g. if max.burst.length is too small).
    episodes <- NULL
  } else {
    episodes <- .burstmatrix.to.episode(flat, channels=s$channels)
  }
      
  episodes
}
                       

.burstmatrix.to.episode <- function(bursts, channels) {
  ## Given a BURST matrix, sort it according to the burst start time
  ## and divide it up into smaller episodes, such that each block
  ## represents the bursts within one episode.


  ## An episode is defined as all the bursts that overlap in time;
  ## e.g. if burst A and B overlap and burst B and C overlap, even
  ## though A and C may not overlap, A,B,C form an episode.
  ##
  ## By definition, if one electrode fires and does not overlap with
  ## any other burst, that counts as an episode, with just one burst.
  ##
  ##

  ## For each episode, we then compute a few statistics, such as the
  ## start and end time of the episode, and the number of electrodes
  ## recruited.


  get.stats <- function(be, channels) {
    burstinfo <- bursts[be[1]:be[2],,drop=FALSE]
    beg <- burstinfo[1,"beg"]
    end <- max(burstinfo[,"end"])
    ## Finding which electrodes were active requires checking against
    ## list of all channels, as the channel names are used, which are characters.
    ##f <- factor(c("b", "b", "d"), levels=c("d", "c", "b", "a"))
    ## table(f)
    active.channels <- table(factor(burstinfo[,"id"], levels=channels))
    n.active <- length(unique( burstinfo[,"id"]))
    if (n.active != sum(active.channels > 0))
      browser()
    c(active.channels, n.active=n.active, beg, end=end)
  }

  ## Sort the burst information according to the burst onset time.
  ## This makes episode detection straightforward: bursts i and i+1 do
  ## not overlap (and hence do not belong to the same episode) if the
  ## start time for burst i+1 is greater than the end time for burst i.
  
  bursts <- bursts[order(bursts[,"beg"]), ]



  B <- nrow(bursts)
  


  new.episode <- TRUE
  looking <- TRUE
  i <- 1

  ## worst case: each burst is on its own, so max number
  ## of episodes is equal to number of bursts
  breaks <- rep(NA, B)                    #store start of new episode
  ##breaks[1] <- 1
  j <- 1
  end <- bursts[1, "end"]

  while (looking) {

    if ( bursts[i+1, "beg"] <= end) {
      ## keep the chain going: possible extend end of episode.
      end = max(end, bursts[i+1, "end"])
    } else {
      ## come to the end of an episode.
      breaks[j] <- i
      j <- j + 1
      end <- bursts[i+1, "end"]
    }    

    i <- i + 1
    if (i == B) {
      ## finished all the bursts.
      breaks[j] <- i
      looking <- FALSE
    }
  }

  ends <- breaks[1:j]
  starts <- c(1, 1+ends[-length(ends)])
  episode.begend <- cbind(starts, ends)

  episodes <- apply(episode.begend, 1, get.stats, channels)

  ## Was previously a list, but better as a matrix for thresholding I think.
  ## res <- cbind(beg=episodes[1,],
  ##             end=episodes[2,],
  ##             num.channels=episodes[3,])
  ##res
  t(episodes)
}


.extract.burst.times.with.id <- function(spikes, burstinfo, id) {
  ## Given the beg+end indexes of a burst, lookup the times of the burst
  ## and return this along with the channel name.
  nbursts <- nrow(burstinfo)
  if(!is.null(nbursts)) {
    begt = spikes[burstinfo[,"beg"]]
    endt = spikes[burstinfo[,"end"]]
    id2 = rep(id, nbursts)
  }
  else {
    ## no bursts
    begt = NA
    endt = NA
    id2 = id
  }
  cbind(begt, endt, id2)
}


