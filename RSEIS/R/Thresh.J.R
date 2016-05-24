`Thresh.J` <-
function(y, thresh)
  {
      x = 1:length(y)
      ##  determine cut off for ratio curve
      ##  this may not be a great idea:  Athresh = min(thresh, 0.8*max(y))
       Athresh = thresh
      k = y>Athresh
      G = rep(0,length(y))
      G[k] =  y[k]
      if(length(y[k])<1) { return(NULL) }
      h = x[k]
      dd = diff(h)
      wd = which(dd>1)
      nw = length(wd)

      ## if there is only one sequence (all ones), return the first value of that sequence
      if(nw<1) { return(list(J=h[1], L=h[1])) }
      J = c(h[1], h[wd[1:(length(wd))]+1])
      
      L = c(h[wd[1:(length(wd))]],  h[length(h)])

      

      return(list(J=J, L=L) )
  }

