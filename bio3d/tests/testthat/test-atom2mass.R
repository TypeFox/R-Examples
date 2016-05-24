context("Testing atom mass functions")


test_that("atom to mass tests", {

  ## Simple test
  atom.names <- c("CA", "O", "N", "OXT")
  ##masses <- c(12.01, 16.00, 14.01, 16.00)
  masses <- c(12.0107, 15.9994, 14.0067, 15.9994)
  expect_that(atom2mass(atom.names), equals(masses, tolerance=1e-6))

  ##masses <- c(42.02, 16.00)
  masses <- c(42.0168, 15.9994)
  expect_that(as.numeric(atom2mass(atom.names, grpby=c(1,1,1,2))), 
              equals(masses, tolerance=1e-6))
  
  ## Should end with error  
  atom.names <- c("CA", "O", "N", "OXT", "CL2", "PT1")
  expect_that(atom2mass(atom.names, rescue=FALSE), throws_error())
  expect_that(atom2mass(atom.names, rescue=TRUE), gives_warning())


  ## Simple test with PDB ID 1HEL
  file <- system.file("examples/1hel.pdb",package="bio3d")
  invisible(capture.output(pdb <- read.pdb(file)))
  
  invisible(capture.output(prot.inds <- atom.select(pdb, "protein")))
  invisible(capture.output(pdb.prot <- trim.pdb(pdb, prot.inds)))

  eletys <- pdb$atom$elety[ pdb$atom$type=="ATOM" ]
  expect_that(sum(atom2mass( eletys )), equals(13346.39, tolerance=1e-6))
  expect_that(sum(atom2mass( eletys )), equals(sum(atom2mass(pdb.prot), tolerance=1e-6)))
  expect_that(sum(atom2mass( pdb.prot )), equals(13346.39, tolerance=1e-6))
  
  ## Try center of mass at the same go
  coma <- c(-0.4991111, 20.5858389, 19.2604674)
  expect_that(c(com(pdb.prot)), equals(coma, tolerance=1e-6))
  
#  coma <- c(-0.5829897, 20.5306061, 19.1081465)
#  expect_that(com(pdb), equals(coma, tolerance=1e-6))
  
  ## Add custom masses
  atom.names <- c("CA", "O", "N", "OXT", "CL2", "PT1")
  masses <- c(12.0107, 15.9994, 14.0067, 15.9994, 35.4530, 195.0780)

  elety.cust <- data.frame(name = c("CL2","PT1"), symb = c("Cl","Pt"))
  ##mass.cust  <- data.frame(symb = c("Cl","Pt"), mass = c(35.45, 195.08))
  
  expect_that(atom2mass(atom.names, elety.custom=elety.cust),  equals(masses, tolerance=1e-6))
  
  ## mass from formula
  form <- "C5 H6 N O3"
  masses <- c(60.050,  6.048, 14.010, 48.000)
  masses <- c(60.05350, 6.04764, 14.00670, 47.99820)
  expect_that(formula2mass(form, sum.mass=FALSE),
              equals(masses, tolerance=1e-6))

  form <- "C5H6"
  expect_that(formula2mass(form), throws_error())

}
          )
