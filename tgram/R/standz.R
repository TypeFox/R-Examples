standz <-
function(tgl1, G=30){
   cacho <- length(tgl1)
   cosa <- rep(tgl1, each=G)
   r <- NULL
   secu <- seq(1, length(cosa), by=cacho)
    for (i in secu){
      r <- c(r, mean(cosa[i:(i+cacho-1)]))
    }
  return(r)
}

