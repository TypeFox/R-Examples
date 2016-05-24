#returns vector from lo to hi with multiplication steps of incr. Used for making ticks to a log-scaled axis 
seqLog <- function(from, to, base=10){base^(log(from,base):log(to,base))}
