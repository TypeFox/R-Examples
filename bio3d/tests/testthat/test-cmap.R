context("Testing cmap function")


test_that("cmap() works properly", {
  ## Simple test with PDB ID 1HEL
  file <- system.file("examples/1hel.pdb",package="bio3d")
  invisible(capture.output(pdb <- read.pdb(file)))
 
  ## Calculate contact map on a small protein
  invisible(capture.output(inds <- atom.select(pdb, "protein")))
  invisible(capture.output(cm <- cmap(pdb$xyz[, inds$xyz],
                                      grpby=pdb$atom[inds$atom, "resno"], ncore=1)))
  
  expect_equal(length(which(cm==1)), 285)
  expect_true(all(is.na(cm[1,1:3])))
  expect_equal(cm[1,4], 0)
  expect_equal(cm[13,129], 1)

  ## Check multicore cmap 
  skip_on_cran()
  trjfile <- system.file("examples/hivp.dcd", package="bio3d")
  invisible(capture.output(trj <- read.dcd(trjfile)))
  invisible(capture.output(cm <- cmap(trj, dcut=6, ncore=1)))
  invisible(capture.output(cm.mc <- cmap(trj, dcut=6, ncore=NULL)))
  expect_that(cm, equals(cm.mc, tolerance=1e-6))
})
