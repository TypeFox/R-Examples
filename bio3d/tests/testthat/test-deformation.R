context("Testing deformation analysis")

test_that("still works", {
  
  ## Simple test with PDB ID 1HEL
  file <- system.file("examples/1hel.pdb",package="bio3d")
  invisible(capture.output(pdb <- read.pdb(file)))
  invisible(capture.output(modes <- nma(pdb)))
  
  #sums0 <- c(59.89283, 141.39431, 109.09525, 122.52931, 172.63766, 317.01506)
  sums0 <- c(180.9078, 198.6242, 318.4639, 379.9139, 479.9795, 473.1810)
  
  defe <- deformation.nma(modes)
  expect_that(defe$sums[1:6], equals(sums0, tolerance=1e-6))
  expect_that(defe$sums[1:6], equals(colSums(defe$ei[,1:6]), tolerance=1e-6))
   
})

test_that("fits with MMTK", {
  
  "calpha.mmtk" <- function(r, ...) {
    ## MMTK Units: kJ / mol / nm^2
    a <- 128; b <- 8.6 * 10^5; c <- 2.39 * 10^5;
    ifelse( r<4.0,
           b*(r/10) - c,
           a*(r/10)^(-6) )
  }

  ## Calc modes
  file <- system.file("examples/1hel.pdb",package="bio3d")
  invisible(capture.output(pdb <- read.pdb(file)))

  invisible(capture.output(modes <- nma(pdb, pfc.fun=calpha.mmtk,
                                        addter=FALSE, mmtk=TRUE)))

  ## deformation energies of mode 7 (using MMTK - with PDB id 1etl)
  #def.mmtk <- c(1306.17014108, 524.571239022, 66.6665951865, 820.62710645,
  #              154.703500149, 754.482784094, 382.993752804, 173.118373857,
  #              287.880418213, 205.968139938, 466.277540766, 814.845931887)
  
  def.mmtk <- c(38.416002, 9.468705, 36.652248, 23.372066, 28.379588,  22.746524,
                35.267401,  58.006941,  48.556190,  46.155725,  92.189766,  75.059341)
  
  ## calc deformation energies
  defe <- deformation.nma(modes, mode.inds=seq(7,26), pfc.fun=calpha.mmtk)
  expect_that(defe$ei[1:12,1], equals(def.mmtk, tolerance=1e-6))

  # mode 8
  def.mmtk <- c(92.87263, 142.04833, 208.63627,  77.01778)
  expect_that(head(defe$ei[,2], n=4), equals(def.mmtk, tolerance=1e-6))

  #mode 9
  def.mmtk <- c(250.2483, 183.0401, 362.0342, 255.6288)
  expect_that(head(defe$ei[,3], n=4), equals(def.mmtk, tolerance=1e-6))

})

