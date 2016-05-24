get.window.estimate <-
function(dataset, seed, iterations=5){
  win <- vector(length=iterations)
  for(i in 1:iterations){
    set1 <- get.cycle.breaks(dataset, window=seed)
    cyclemat <- get.cycle.matrix(dataset, set1$CycleBreaks)
    cycledurs <- get.cycle.durations(cyclemat)
    win[i] <- get.window(cycledurs)
    print(win[i])
    seed <- win[i]
    if(i > 1){
      if(win[i] == win[i-1]){
        WindowEst <- win[i]
        break
      }
    }
  }
  return(WindowEst)
}

