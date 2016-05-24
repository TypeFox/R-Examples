insertNAs <-
function(v, w)
  {
    nw = length(w)
    nv = length(v)
    vprime = NULL
    
    ###  vprime[1:(w[1]-1)] = v[1:(w[1]-1)]
    istart = 1
    for(i in 1:nw)
      {
        iend = w[i]
        temp = v[istart:iend]
        vprime = c(vprime, temp, NA)
        istart = w[i]+1
      }
    istart = w[nw]+1
    iend = nv
    temp = v[istart:iend]
    vprime = c(vprime, temp)

    return(vprime)

  }

