# For aligning peaks across multiple templates so templates can be compared
# Modified: 2013 DEC 05

compareTemplates <-
function(
  detection.obj,  # A detectionList object
  cutoff.return,  # Events with no values above will not be returned
  cutoff.ignore,  # Peaks with values below will be dropped at the start
  tol,            # Tolerance (s). If a peak is within tol of a peak from another template, they are in the same event
  n.drop=0        # Rows with this many templates or fewer will be dropped. n.drop=0 drops none.
) {

  # To avoid CRAN check notes
  score <- template <- NULL

  # Check arguments
  if(missing(detection.obj)) stop('Required argument detection.obj is missing')
  if(missing(cutoff.return)) {
    i <- which.min(lapply(detection.obj@templates, slot, "score.cutoff"))
    warning('cutoff.return argument missing, so using score.cutoff from template \"', names(detection.obj@templates)[i], '\", which is ', cutoff.return <- detection.obj@templates[[i]]@score.cutoff, '.')
  }
  if(missing(cutoff.ignore)) {
    i <- which.min(lapply(detection.obj@templates, slot, "score.cutoff"))
    warning('cutoff.ignore argument missing, so using 1/10th score.cutoff from template \"', names(detection.obj@templates)[i], '\", which is ', cutoff.ignore <- detection.obj@templates[[i]]@score.cutoff/10, '.')
  }
  if(missing(tol)) {
    i <- which.min(lapply(detection.obj@templates, slot, "duration"))
    warning('tol argument missing, so using the template duration from \"', names(detection.obj@templates)[i], '\", which is ', tol <- signif(detection.obj@templates[[i]]@duration, 3), ' s.')
  }
  if(tol<0.001) warning('Did you really mean for tol argument to be so low (', tol, ')? If you see errors try a higher value.')

  n.templates <- length(detection.obj@templates)
  if(n.templates<2) stop('The detection.obj argument, ', deparse(substitute(detection.obj)), ', has results for only one template--there is nothing to compare it to.')
  pks <- getPeaks(detection.obj)
  pks <- subset(pks, score>cutoff.ignore)
  n.pks <- nrow(pks)
  pks.L <- list()
  for(i in names(detection.obj@templates)) {
    pks.L[[i]] <- subset(pks, template == i)
  }

  # Combine all peaks together initially with no merging--one template per row
  pk.idx.mat <- times.mat <- matrix(NA, nrow=n.pks, ncol=n.templates, dimnames=list(time.mean=1:n.pks, template=names(detection.obj@templates)))
  n.in <- 0
  for(i in names(detection.obj@templates)) {
    pks <- pks.L[[i]]
    n.pks <- nrow(pks)
    times.mat[1:n.pks + n.in, i] <- pks$time 
    pk.idx.mat[1:n.pks + n.in, i] <- 1:n.pks
    n.in <- n.in + n.pks
  }

  t.mean <- rowMeans(times.mat, na.rm=TRUE)

  # Sort matrices by mean time
  times.mat <- times.mat[order(t.mean), ]
  pk.idx.mat <- pk.idx.mat[order(t.mean), ]
  t.mean <- rowMeans(times.mat, na.rm=TRUE)
  n.events <- nrow(times.mat)

  # Now go through the times.mat matrix row by row, and then cell by cell, moving times to the row below if they are close to the mean
  dropped <- TRUE
  while(dropped) {
    dropped <- FALSE 
    i <- 1
    while(i<nrow(times.mat)) {
      for(j in 1:n.templates) {
        if(!is.na(times.mat[i, j]) && is.na(times.mat[i+1, j]) && abs(times.mat[i, j] - t.mean[i+1])<tol) {
          times.mat[i+1, j] <- times.mat[i, j]
          times.mat[i, j] <- NA
          pk.idx.mat[i+1, j] <- pk.idx.mat[i, j]
          pk.idx.mat[i, j] <- NA
        }
      }
      # Drop the row i if only NAs are left 
      if(sum(is.na(times.mat[i, ])) == n.templates) {
        times.mat <- times.mat[-i, ]
        pk.idx.mat <- pk.idx.mat[-i, ]
        dropped <- TRUE
      }
      t.mean <- rowMeans(times.mat, na.rm=TRUE)
      i <- i + 1
    }
  }

  n.events <- nrow(times.mat)

  # Put together matrices of scores and absolute time
  scores.mat <- date.time.mat <- pk.idx.mat
  #colnames(scores.mat) <- colnames(date.time.mat) <- colnames(times.mat) <- names(pks.L)
  dimnames(scores.mat) <- dimnames(date.time.mat) <- dimnames(times.mat) <- list(time.mean=t.mean, template=names(detection.obj@templates))
  for(i in 1:n.templates) {
    t <- pks.L[[i]]$date.time
    s <- pks.L[[i]]$score
    date.time.mat[, i] <- as.character(t[pk.idx.mat[, i]])
    scores.mat[, i] <- s[pk.idx.mat[, i]]
  }

  # Find those rows with at least one score above the cutoff and enough templates to not be dropped
  which.rows <- apply(scores.mat, 1,max, na.rm=TRUE)>cutoff.return & apply(times.mat, 1,function(x) sum(!is.na(x)))>n.drop

  # Trim down all matrices
  times.mat <- times.mat[which.rows, ]
  scores.mat <- scores.mat[which.rows, ]
  date.time.mat <- date.time.mat[which.rows, ]
  t.mean <- rowMeans(times.mat, na.rm=TRUE)
  rownames(scores.mat) <- rownames(date.time.mat) <- rownames(times.mat) <- round(t.mean, 3)

  # RETURNS A LIST
  return(list(times.mean=t.mean, times=times.mat, date.time=date.time.mat, scores=scores.mat))
}
