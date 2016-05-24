context("Testing pca()")

test_that("pca functions works", {
  "mysign" <- function(a,b) {
    if(all(sign(a)==sign(b)))
      return(1)
    else
      return(-1)
  }

  attach(transducin)
  inds <- unlist(lapply(c("1TND_A", "1TAG", "1AS0", "1AS2"), grep, pdbs$id))
  pdbs <- trim.pdbs(pdbs, row.inds=inds)
  gaps <- gap.inspect(pdbs$xyz)

  ## Calc modes
  invisible(capture.output(pc <- pca(pdbs)))

  ## check dimensions
  expect_that(dim(pc$U), equals(c(939, 939)))
  expect_that(length(pc$L), equals(939))

  ## check eigenvalues
  Lexpected <- c(1.964689e+02, 1.715903e+02, 7.091482e+01)

  expect_that(head(pc$L, n=3), equals(Lexpected, tolerance=1e-6))

  ## check atom-wise loadings
  AUexpected <- c(0.013422274, 0.023879443, 0.022779307, 0.023490288, 0.006308588, 0.010057694)
  expect_that(head(pc$au[,1], n=6), equals(AUexpected, tolerance=1e-6))

  Z1expected <- c(-3.555218, -12.081135, -4.602888, 20.239240)
  Z1now <- as.numeric(head(pc$z[,1], n=4))
  expect_that(Z1now * mysign(Z1expected, Z1now), equals(Z1expected, tolerance=1e-6))

  Mexpected <- c(30.12193, 67.76449, 43.36594, 27.01919, 69.66411, 44.47434)
  expect_that(head(pc$mean, n=6), equals(Mexpected, tolerance=1e-6))

  detach(transducin)
})
