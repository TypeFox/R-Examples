
# If is.pc then a value like "1" means 0.01
remap.cumsum <- function(z, is.pc=TRUE, base=0) {
  for (i in 1:ncol(z)) {
    tmp <- z[,i]
    if (is.pc) {
      tmp <- tmp/100
    }
    z[,i] <- base+cumsum(tmp)
  }
  z
}
