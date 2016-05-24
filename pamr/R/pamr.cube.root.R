## Cube root transformation for Affy chips
pamr.cube.root  <- function(x) {
  return(sign(x) * abs(x)^{1/3})
}

