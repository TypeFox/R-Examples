FillInGaps <-
function(dataset, impt.times) {
  Gaps <- function(observed.times, impt.times) {
    max.time <- max(observed.times)
    min.time <- min(observed.times)
    fill.times <- impt.times[impt.times > min.time & impt.times < max.time]    
    if(any(fill.times %in% observed.times))
      fill.times <- fill.times[ - which(fill.times %in% observed.times)]
    return(fill.times)
  }  
  tp.gaps <- tapply(dataset[, 2], dataset[, 1], function(x) Gaps(x, impt.times))
  no.missing.obs <- unlist(lapply(tp.gaps, length))
  tp.gaps <- unlist(tp.gaps)
  unique.id <- unique(dataset[, 1])
  id.gaps.long <- rep(unique.id, no.missing.obs) 
  dataset.gaps <- as.data.frame(matrix(NA, nrow=length(tp.gaps), 
                                       ncol=dim(dataset)[2]))
  dataset.gaps[, 1] <- id.gaps.long
  dataset.gaps[, 2] <- tp.gaps
  names(dataset.gaps) <- names(dataset)
  dataset.gaps$obs.type <- rep(2, dim(dataset.gaps)[1])
  new.dataset <- rbind(dataset, dataset.gaps)
  new.dataset <- new.dataset[order(new.dataset[, 1], new.dataset[, 2]), ]
  return(new.dataset)
}
