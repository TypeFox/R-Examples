# For finding peaks in output from corMatch or binMatch and identifying detections
# Modified: 2015 APR 2

findPeaks <- 
function(
  score.obj,       # A scoreList object, complete output from corMatch or binMatch
  fd.rat=1,        # Factor to multiply template duration by for determining frame width
  frame,           # Or directly specify the frame width
  parallel=FALSE
) {

  if(missing(score.obj)) stop('Required argument score.obj is missing')

  # Required packages
  if(parallel) {
    lapplyfun <- function(X, FUN, mc.cores) parallel::mclapply(X, FUN, mc.cores=parallel::detectCores())
  } else lapplyfun <- lapply

  # Finally start working on peaks. 
  score.names <- as.list(names(score.obj@templates))
  names(score.names) <- c(score.names)
  # Create an object to hold results.
  results <- list()

  # Work through all templates
  if(missing(frame)) frame <- NA
  results <- lapplyfun(
    X=score.names, 
    FUN=function(i) {
      dat <- score.obj@scores[[i]]
      time <- dat$time
      score <- dat$score
      if(is.na(frame)) frame <- fd.rat*score.obj@templates[[i]]@duration
      score.cutoff <- score.obj@templates[[i]]@score.cutoff

      # Time step between points in time, score data
      time.step <- time[2]-time[1]

      # Convert frame from seconds to number of time bins, and call span
      span <- min(length(score), floor(frame/time.step)) # embed function won't work for span > length(x)
      # Make span odd--ensures that peaks (and not peak neighbors) are returned  
      if(span%%2 == 0) span <- span+1
      halfspan <- span%/%2 # Should always be (span-1)/2, so the numbers of bins on either side of point to check.
      # Extend score data so points close or at ends can be identified as peaks
      score.extended <- c(rep(0, halfspan), score, rep(0, halfspan))
      lagmat <- embed(score.extended, span) # Makes a matrix where score vector is lagged by one position between columns
      result <- max.col(lagmat) == 1 + halfspan # Returns TRUE when center column of matrix has the maximum value
      pks <- dat[result, ] 
      rownames(pks) <- 1:nrow(pks)

      pks$detection <- pks$score>=score.cutoff
      hits <- pks[pks$score>=score.cutoff, ]
      hits$detection <- NULL
      if(nrow(hits)>0) rownames(hits) <- 1:nrow(hits)

      # Summarize peak results
      summary <- c(n.peaks=round(nrow(pks)), n.hits=round(nrow(hits)), max.score=signif(max(score), 3), min.score=signif(min(score), 3))

      cat('\nDone with ', i)
      return(list(peaks=pks, detections=hits))
    }
  )

  # Separate peaks and detections
  peaks <- lapply(results, `[[`, 'peaks')
  detections <- lapply(results, `[[`, 'detections')

  object <- new('detectionList', survey.name=score.obj@survey.name, survey=score.obj@survey, survey.data=score.obj@survey.data, templates=score.obj@templates, scores=score.obj@scores, peaks=peaks, detections=detections)
  cat('\nDone\n')
  return(object)
}


