# matrix of all patterns
patternmat3<-function(nobj)      #pattStr<-names(countpattern(dat[1,]))
{
      nvar<-nobj*(nobj-1)/2
      Y <- matrix(0, 3^nvar, nvar)
      for (i in 1:nvar) Y[, nvar + 1 - i] <- rep(rep(c(0, 1, 2), c(3^(i -
          1), 3^(i - 1),3^(i-1))), 3^(nvar - i))
      Y
}
