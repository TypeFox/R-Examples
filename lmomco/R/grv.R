"prob2grv" <-
function(f) {
   if(! check.fs(f)) return()
   return(-log(-log(f)))
}

"grv2prob" <-
function(grv) {
   return(exp(-exp(-1*grv)))
}

