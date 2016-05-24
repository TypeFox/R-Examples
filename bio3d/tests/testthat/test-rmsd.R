context("Testing RMSD function")

test_that("rmsd() gets the same results as PyMOL", {
  file <- system.file(c("examples/1hel.pdb", "examples/1dpx.pdb"), package="bio3d")

   invisible(capture.output(pdb.a <- read.pdb(file[1])))
   invisible(capture.output(pdb.b <- read.pdb(file[2])))

   invisible(capture.output(inds.a <- atom.select(pdb.a, "calpha")))
   invisible(capture.output(inds.b <- atom.select(pdb.b, "calpha")))

  rd1 <- rmsd(a=pdb.a$xyz, b=pdb.b$xyz, a.inds=inds.a$xyz, b.inds=inds.b$xyz, fit=FALSE)
  rd2 <- rmsd(a=pdb.a$xyz, b=pdb.b$xyz, a.inds=inds.a$xyz, b.inds=inds.b$xyz, fit=TRUE)
  
  ## with pymol "pair_fit 1hel and name CA, 1dpx and name CA" = 0.293
  rmsd1 <- 1.386
  rmsd2 <- 0.293

  expect_equal(round(rd1, 3), rmsd1)
  expect_equal(round(rd2, 3), rmsd2)
})

test_that("rmsd() with ncore>1 works properly", {

  file <- system.file(c("examples/1hel.pdb", "examples/1dpx.pdb"), package="bio3d")
  invisible(capture.output(pdb.a <- read.pdb(file[1])))
  invisible(capture.output(pdb.b <- read.pdb(file[2])))
  
  invisible(capture.output(inds.a <- atom.select(pdb.a, "calpha")))
  invisible(capture.output(inds.b <- atom.select(pdb.b, "calpha")))
     
  ## check if ncore > 1 is really faster 
  time1 <- system.time(rmsd1 <- rmsd(a=pdb.a$xyz, b=pdb.b$xyz, a.inds=inds.a$xyz, b.inds=inds.b$xyz, fit=TRUE, ncore=1))
  time2 <- system.time(rmsd2 <- rmsd(a=pdb.a$xyz, b=pdb.b$xyz, a.inds=inds.a$xyz, b.inds=inds.b$xyz, fit=TRUE, ncore=NULL))
  
  ##time1 <- time1["elapsed"]
  ##time2 <- time2["elapsed"]
  
#  expect_equivalent(rmsd1, rmsd2)
  expect_equal(rmsd1, rmsd2, tolerance=1e-6)
  
#   cat("Speed up by", round((time1-time2)/time2, 1)*100, "%", sep="")
#   if(getOption("cores") > 1)
#      expect_true(time1 > time2)
})

