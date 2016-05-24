"nonexceeds" <-
function(f01=FALSE, less=FALSE, sig6=FALSE) {
  if(sig6) {
     f <- pnorm(seq(0,6,by=0.02))
  } else {
     if(! less) {
        f <- c(0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.96,0.98,0.99,0.996,0.998)
     } else {
        f <- c(0.6,0.7,0.8,0.9,0.95,0.98,0.99,0.995,0.998, 0.999)
     }
  }
  cf <- sort(1-f)
  f <- c(cf,0.5,f)
  f <- unique(f) # just incase duplicates get in (sig6)
  ifelse(f01, return(c(0,f,1)), return(f))
}
