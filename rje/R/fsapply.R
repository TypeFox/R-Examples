fsapply <-
function (x, FUN) 
unlist(lapply(x, FUN))
