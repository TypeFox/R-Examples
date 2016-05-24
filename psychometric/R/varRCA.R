"varRCA" <-
function(x, aprox=FALSE)
 {
 if (!aprox) {vrt <- varResT(x)}
 else {vrt <- varResT(x, T)}
 aa <- CAFAA(x)
 vr <- vrt/aa^2
 return(vr)
 }

