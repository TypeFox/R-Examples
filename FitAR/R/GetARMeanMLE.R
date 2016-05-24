`GetARMeanMLE` <-
function(z, phi){
 stopifnot (length(z)>=2*length(phi))
 g1<-Get1G(phi, length(z))
 sum(g1*z)/sum(g1)
 }

