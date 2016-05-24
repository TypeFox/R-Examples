"prob2lrv" <-
function(f) {
   if(! check.fs(f)) return()
   return(-log((1-f)/f))
}

"lrv2prob" <-
function(lrv) {
   return(1/(exp(-1*lrv) + 1))
}
