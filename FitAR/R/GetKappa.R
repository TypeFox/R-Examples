`GetKappa` <-
function(phi)
 rev(cumsum(rev(TacvfMA(phi, length(phi))[-1])))

