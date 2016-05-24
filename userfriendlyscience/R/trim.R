### trim simply trims spaces from the start and end of a string
trim <- function(str) {
  ### Based on 'trim' in package Gdata by
  ### Gregory R. Warnes <greg at warnes.net> and others
  str <- sub(pattern="^ +", replacement="", x=str)
  str <- sub(pattern=" +$", replacement="", x=str)
  str <- sub(pattern="^\t+", replacement="", x=str)
  str <- sub(pattern="\t+$", replacement="", x=str)
  return(str);
}