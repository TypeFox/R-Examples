`hilow` <-
function(y)
{
 #####  find the peaks and valleys of the signals return the indecies of the locations 

  ind = seq(1,length(y))
  w = y

  aa = peaks(w, span=3)
 
  ab = peaks((-1)*w, span=3)
  
  ka = rep(FALSE, length(aa)+2)
  ka[2:(length(aa)+1)] = aa
  pks = ind[ka]
 ###  print(pks)
  
  kb = rep(FALSE, length(ab)+2)
  kb[2:(length(ab)+1)] = ab 
  vals = ind[kb]
###   print(vals)

  return(list(hi=pks, lo=vals))

}

