"CredIntRho" <-
function(x,  aprox=FALSE, level=.95)
 {
 r <- rhoCA(x)
 if (!aprox) {
 vr <- varRCA(x)}
 else { vr <- varRCA(x,T)} 
zs <- - qnorm((1-level)/2)
sdr <- sqrt(vr)
lcl <- r - zs * sdr
ucl <- r + zs * sdr
return(list(lcl,ucl))
}

