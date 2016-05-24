getreadcounts <-
function(hits, windowlength, chrlength){
  windows <- seq(1,chrlength,windowlength)
  readcount <- numeric(length=length(windows))
  count <- 0
  i <- 1
  hits <- sort(hits)
  for (hit in hits){
    if (!isbetween(hit,windows[i],(windows[i]+windowlength-1))){
      readcount[i] <- count
      count <- 0
      while(!isbetween(hit,windows[i],(windows[i]+windowlength-1))){
        i <- i + 1
      }
    }
    count <- count+1
  }
  readcount[i] <- count
  return(readcount)
}
