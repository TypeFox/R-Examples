##' Compute the tiling coefficient
##'
##' The tiling coefficient is computed for two spike trains A and B.  Spike
##' trains must be sorted (smallest first.)
##
##' @param a first spike train (vector of spike times)
##' @param b second spike train
##' @param dt Time-window (in seconds) for coincident activity.
##' @param rec.time 2-vector stating beg and end of recording.  If omitted
##' we take the min(max) spike time as the beg(end).
##' @return Tiling coefficient, which should be between -1 and 1.
##' @author Catherine Cutts and Stephen Eglen
##' @useDynLib IGM.MEA, .registration=TRUE, .fixes = "C_"
.tiling.corr <- function(a, b, dt=0.05, rec.time=NULL) {
  ## Return the tiling correlation between two spike trains A and B.
  ## DT is the key parameter.
  ##
  ## rec.time should be a length-2 vector stating the beginning and end
  ## time of the recording.  If NULL, beg (end) is set to the time of
  ## the first (last) spike in either train.
  if(is.null(rec.time)) {
    rec.time <- range( c(a,b))
  }

  z <- .C("run_TM",
          as.integer(length(a)),
          as.integer(length(b)),
          as.double(dt),
          rec.time,
          coeff=double(1),              #return value
          as.double(a),
          as.double(b))
  z$coeff
}

  
  
.tiling.allpairwise.old <- function(s, dt=0.05) {
  ## Return matrix of all pairwise tiling correlations from object s.
  ## TODO: only compute upper triangular elements completed.
  ## This is slow as lots of looping happens in R.
  N <- s$NCells
  indices=array(0,dim=c(N,N))

  rec.time <- s$rec.time
  for(i in 1:N) {
    ni <- s$nspikes[i]
    for(j in 1:N) {
      nj <- s$nspikes[j]
      z <- .C("run_TM",
              as.integer(ni),
              as.integer(nj),
              as.double(dt),
              as.double(rec.time),
              coeff=double(1),
              as.double(s$spikes[[i]]),
              as.double(s$spikes[[j]]))
      indices[i,j] <- z$coeff
    }
  }
  indices
}


##' Compute tiling coefficient for an MEA recording.
##'
##' Given an s object, we return all pairwise correlations.
##' @param s  The spike object.
##' @param dt Time-window (in seconds) for coincident activity.
##' @return Upper triangular matrix of tiling coefficients.
##' @author Stephen Eglen
.tiling.allpairwise <- function(s, dt=0.05) {
  n <- length(s$spikes)

  all.spikes <- unlist(s$spikes)
  nspikes <- sapply(s$spikes, length)
  first.spike <- c(0, cumsum(nspikes)[-n])
  z <- .C("tiling_arr",
          as.double(all.spikes),
          as.integer(n),
          as.integer(nspikes),
          as.integer(first.spike),
          as.double(s$rec.time),
          as.double(dt),
          res = double(n*n))

  ## return the result.
  m <- array(z$res, dim=c(n,n))
  m[lower.tri(m)] <- NA                 #we didn't do lower triangle, so ignore those.
  m
}
