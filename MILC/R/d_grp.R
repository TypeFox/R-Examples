d_grp <-
function(d)
 {
  d_char <- c("NA", "1-10", "11-20", "21+")
  d_num  <- c(1, 11, 21) 
  if (is.na(d)) {dgrp=NA} else if (d>=max(d_num)) {dgrp <- d_char[length(d_char)]} else 
  {ind_d  <- min(which(d_num>d)) ; dgrp <- d_char[ind_d]}
  return(dgrp)
 }
