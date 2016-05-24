"to.double" <-
function(x) {
# warning: does not convert single to double in R
   storage.mode(x) <- "double";
   x}

