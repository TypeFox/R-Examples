"pvaaa" <-
function(x, aprox=FALSE)
 {
 if (!aprox) {ve <- vare(x)}
 else {ve <- aprox.vare(x)}
vr <- varr(x)
vav <- varAV(x)
 pv <- (ve+vav)/vr*100
 return(pv)
 }

