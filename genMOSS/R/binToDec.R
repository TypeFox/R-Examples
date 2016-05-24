binToDec <-
function (x) {

  dec <- sum(x * 2^(rev(seq_along(x)) - 1))
  return(dec)  

}
