# matrix of all patterns
patternmat2<-function(nobj)      #pattStr<-names(countpattern(dat[1,]))
{
      nvar<-nobj*(nobj-1)/2
      Y <- matrix(0, 2^nvar, nvar)
      for (i in 1:nvar) Y[, nvar + 1 - i] <- rep(rep(c(0, 1), c(2^(i -
          1), 2^(i - 1))), 2^(nvar - i))
      Y
}
