Mould_dur <-
function(possible.fix, Hz, plot = FALSE, default = 100){
  fixa <- possible.fix
  dat <- log(fixa / (Hz/1000))
  nsd <- seq(1, max(dat), length.out = round(max(fixa / (Hz/1000)) / 4))
  freq <- sapply(1:round(max(fixa / (Hz/1000))/4), function(x) length(which(dat >= nsd[x] & dat < nsd[x+1])))
  nsd <- nsd[freq != 0]
  freq <- freq[freq != 0]
  
  h <- .2
  pred <- predict(loess(freq ~ nsd ,span = h))
  while (length(lomax(pred)) > 2) {
    h <- h + .01
    pred <- predict(loess(freq ~ nsd ,span = h ))
  }
  
  if(length(lomax(pred)) == 2){
    thres_dur <- exp(nsd[which.min(pred[lomax(pred)[1] : lomax(pred)[2]])])
  } else {
    thres_dur <- default
    warning('The duration threshold could not be estimated, due to smoothing problems, therefore the threshold specified in function is used (default = 100ms)')
  }
  
  if(plot == T){
    hist(fixa / (Hz/1000), breaks = 200, col = 'grey30', border = 'grey30', main = '', xlab = 'Non-saccadic duration, ms')
    segments(thres_dur, 30, thres_dur, 0, col = 1, lwd = 5)
    lines(pred ~ exp(nsd), lwd = 6, col = 2)
  }
  ## Set a maximal value to overcome smoothing problems
  if(thres_dur > 250) thres_dur <- 250
  return(thres_dur)
}
