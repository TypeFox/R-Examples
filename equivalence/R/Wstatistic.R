"Wstatistic" <-
function(W, eps1, eps2, VarW)
   return(abs(W - 1/2 - (eps2-eps1)/2)/sqrt(VarW))

