octav <- function(x, oct, preston=FALSE){
  if(is(x, "fitsad"))
    y <- x@data$x
  else if(is(x,"fitrad"))
    y <- x@rad.tab$abund
  else if(is(x,"numeric"))
    y <- x[x>0]
  if(missing(oct)){
    oct <- 0:(ceiling(max(log2(y)))+1)
    if(any(y < 1)){
      octlower <- floor(min(log2((y)))):-1
      oct <- c(octlower, oct)
    }
  }
  else{
      if(min(oct)>min(log2(y))||max(oct)<max(log2(y))) stop(" 'oct' should include all abundance values in 'x' ")
      oct <- min(oct):max(oct)
  }
  oct <- unique(oct)
  N <- 2^(oct)
  oct.hist <- hist(y, breaks=c(0,N), plot=FALSE)
  res <- data.frame(octave = oct, upper = oct.hist$breaks[-1], Freq = oct.hist$counts)
  if(preston) res <- prestonfy(res, y)
  new("octav", res)
}

# Helper function to create Preston octaves. Is also used by octavpred
prestonfy <- function(res, y) {
  N <- 2^(res$octave)
  j <- N[-length(N)]
  w <- y[y%in%j]
  ties <- table(factor(w, levels=j))
  res[-1, 3] <- res[-1, 3]+ties/2
  res[-length(N), 3] <- res[-length(N), 3]-ties/2
  return(res)
}

