"tPsiw" <-
function(r,tl,tu){
nr   <- length(r); tmp  <- rep(0,nr)
ind  <- (1:nr)[r > tl & r < tu] 
tmp[ind] <- exp(r[ind])-1; tmp}

