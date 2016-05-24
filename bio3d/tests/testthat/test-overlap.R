
context("Testing overlap functions")


test_that("Overlap functions", {

  ## Simple test with PDB ID 1HEL
  file <- system.file("examples/1hel.pdb",package="bio3d")
  invisible(capture.output(pdb.a <- read.pdb(file)))

  file <- system.file("examples/1dpx.pdb",package="bio3d")
  invisible(capture.output(pdb.b <- read.pdb(file)))
  
  ## Calculate modes with default arguments
  invisible(capture.output(modes <- nma(pdb.a, inds=NULL, ff='calpha',
                                        mass=FALSE, temp=300.0)))

  ca.inds.a <- atom.select(pdb.a, "calpha", verbose=FALSE)
  ca.inds.b <- atom.select(pdb.b, "calpha", verbose=FALSE)

  ## Set new coordinates of pdb.b
  xyz.b <- fit.xyz(pdb.a$xyz, pdb.b$xyz,
                   fixed.inds=ca.inds.a$xyz,
                   mobile.inds=ca.inds.b$xyz)
  pdb.b$xyz <- xyz.b
  
  ## difference vector
  dv <- difference.vector(rbind(pdb.a$xyz[ca.inds.a$xyz],
                                pdb.b$xyz[ca.inds.b$xyz]))
  
  o1 <- overlap(modes, dv, nmodes=(modes$natoms*3)-6)
  expect_that(o1$overlap.cum[(modes$natoms*3)-6],
              equals(1, tolerance=1e-6))
  
  expect_that(o1$overlap.cum[1], equals(0.2786508, tolerance=1e-6))

  o2 <- overlap(modes$U[,7:26], dv)
  expect_that(all((round(o1$overlap[1:20] - o2$overlap, 10)==0)),
              equals(TRUE))

  ## Calculate modes with default arguments
  invisible(capture.output(modes.b <- nma(pdb.b, inds=NULL, ff='calpha',
                                          mass=FALSE, temp=300.0)))

  r <- rmsip(modes, modes.b)
  expect_that(r$overlap[1,1], equals(0.704, tolerance=1e-6))
  expect_that(r$overlap[1,2], equals(0.286, tolerance=1e-6))
  expect_that(r$overlap[2,1], equals(0.289, tolerance=1e-6))
 
}
          )
