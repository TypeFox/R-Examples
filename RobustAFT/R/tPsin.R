"tPsin" <-
function(r,tl,tu){
# Weight function
nr   <- length(r)
tmp  <- rep(0,nr)
ind  <- (1:nr)[r>tl & r<tu]
tmp[ind] <- r[ind]; tmp}

