generateSudoku <- function(Nblank=50, print.it=FALSE) {
  z <- c(1:9,4:9,1:3,7:9,1:6,2:9,1,5:9,1:4,8:9,1:7,3:9,1:2,6:9,1:5,9,1:8)
  z <- matrix(sample(9)[z], 9,9)
  for (i in 1:5) z <- z[replicate(3, sample(3)) + 3*rep(sample(0:2), each=3),
                        replicate(3, sample(3)) + 3*rep(sample(0:2), each=3)]
  for (bi in seq(0,6,3)) for (bj in seq(0,6,3)) {
    idx <- data.matrix(expand.grid(bi + 1:3, bj + 1:3))
    z[idx[sample(1:9, Nblank%/%9), ]] <- 0
  }
  ## Depopulate (if we had a test for uniqueness, we'd put it here):
  while (sum(!z) < Nblank) z[matrix(sample(9,2), 1)] <- 0
  if (print.it) printSudoku(z)
  z
}
