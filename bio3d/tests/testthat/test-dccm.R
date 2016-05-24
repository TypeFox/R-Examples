context("Testing dccm functions")


test_that("Correlation matrix from NMA", {

  ## Calculate correl mat on a small protein
  file <- system.file("examples/1hel.pdb",package="bio3d")
  invisible(capture.output(pdb <- read.pdb(file)))

  invisible(capture.output(modes <- nma(pdb)))
  invisible(capture.output(cm <- dccm.nma(modes, ncore=1)))
  
  expect_that(cm[1,1], equals(1,           tolerance=1e-6))
  expect_that(cm[1,2], equals(0.4380029,  tolerance=1e-6))
  expect_that(cm[1,3], equals(0.1407395, tolerance=1e-6))
  expect_that(cm[1,3], equals(cm[3,1]))
  expect_that(sum(cm), equals(57.71768, tolerance=1e-6))

  ## Check multicore DCCM  
  invisible(capture.output(cm.mc <- dccm.nma(modes, ncore=NULL)))
  expect_that(cm, equals(cm.mc, tolerance=1e-6))
  
}
          )

test_that("Correlation matrix from XYZ (dccm.xyz)", {
  skip_on_cran()
  ## Calculate correl mat on a short HIV protease simulation
  trjfile <- system.file("examples/hivp.dcd", package="bio3d")
  invisible(capture.output(trj <- read.dcd(trjfile)))
  invisible(capture.output(cm <- dccm(trj, ncore=1)))
  
  expect_that(cm[1,1], equals(1,           tolerance=1e-6))
  expect_that(cm[1,2], equals(0.9965964,  tolerance=1e-6))
  expect_that(cm[1,3], equals(0.992211, tolerance=1e-6))
  expect_that(cm[1,3], equals(cm[3,1]))

  ## Check multicore DCCM
  invisible(capture.output(cm.mc <- dccm(trj, ncore=NULL)))
  expect_that(cm, equals(cm.mc, tolerance=1e-6))
  
})

test_that("Correlation matrix from PCA (dccm.pca)", {
  skip_on_cran()
  ## Calculate correl mat on a short HIV protease simulation
  trjfile <- system.file("examples/hivp.dcd", package="bio3d")
  invisible(capture.output(trj <- read.dcd(trjfile)))
  invisible(capture.output(xyz <- fit.xyz(trj[1, ], trj[1:20, ], 1:ncol(trj), 1:ncol(trj)) )) 
  invisible(capture.output(pca <- pca.xyz(xyz)))
  invisible(capture.output(cm <- dccm(pca, ncore = 1)))
  pca$z <- NULL
  invisible(capture.output(cm2 <- dccm(pca, ncore = 1)))
   
  expect_that(cm[1,1], equals(1,           tolerance=1e-6))
  expect_that(cm[1,2], equals(0.7120510,  tolerance=1e-6))
  expect_that(cm[1,3], equals(0.5455956, tolerance=1e-6))
  expect_that(cm[1,3], equals(cm[3,1]))
  expect_that(cm, equals(cm2, tolerance=1e-6))

  ## Check multicore DCCM
  invisible(capture.output(cm.mc <- dccm(pca, ncore=NULL)))
  expect_that(cm, equals(cm.mc, tolerance=1e-6))
  
})
