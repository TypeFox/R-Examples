"varAV" <-
function(x)
 {
 rho <- rhoCA(x)
 AA <- CAFAA(x)
 cvf <- CVF(x)
 vav <- rho^2*AA^2*cvf
 return(vav)
 }

