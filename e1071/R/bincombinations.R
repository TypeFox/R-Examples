## Kopie in mlbench
bincombinations <- function(p) {

  retval <- matrix(0, nrow=2^p, ncol=p)

  for(n in 1:p){
    retval[,n] <- rep(c(rep(0, (2^p/2^n)), rep(1, (2^p/2^n))),
                      length=2^p)
  }
  retval
}



