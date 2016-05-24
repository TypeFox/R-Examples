
# is.pc and is.returns are TRUE
#    values are like "1" for 1%
# is.returns is true but not is.pc
#    values are like "0.01" for 1%
# is.returns is false                         in this case is.pc is ignored!
#    values are like 1.01 for 1%
remap.cumprod <- function(z, is.pc=TRUE, is.returns=TRUE, base=100) {
  for (i in 1:ncol(z)) {
    tmp <- z[,i]
    if (is.returns) {
      if (is.pc) {
        tmp <- tmp/100
      }
      tmp <- 1+tmp
    }
    tmp[1] <- base
    z[,i] <- cumprod(tmp)
  }
  z
}
