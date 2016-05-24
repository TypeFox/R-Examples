trnd<-
  function(l,u)
  {# uses acceptance rejection to simulate from truncated normal
    x=rnorm(length(l)) # sample normal
    # keep list of rejected
    I=which(x<l|x>u);d=length(I);
    while (d>0) {# while there are rejections
      ly=l[I] # find the thresholds of rejected
      uy=u[I]
      y=rnorm(length(ly))
      idx=(y>ly)&(y<uy) # accepted
      x[I[idx]]=y[idx] # store the accepted
      I=I[!idx] # remove accepted from list
      d=length(I) # number of rejected
    } 
    return(x)
  }
