"varResT" <-
function(x, aprox=FALSE)
 {
 if (!aprox) {ve <- vare(x)}
 else {ve <- aprox.vare(x)}
 vr <- varr(x)
 vav <- varAV(x)
 vrest <- vr - ve - vav
 return(vrest)
 }

