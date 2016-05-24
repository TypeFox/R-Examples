selstas <-
function(sta, ind)
  {
    for(i in 1:length(sta))
      {
        g = sta[[i]]
        sta[[i]] = g[ind]

      }
    return(sta)
  }

