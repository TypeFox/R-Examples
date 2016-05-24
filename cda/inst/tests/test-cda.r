
require(cda)
require(testthat)

context("Checking low-level functions against R implementation")

test_that("Euler rotation matrix is correct", {  
  cpp_result <- cda$euler(pi/3,pi/3,pi/3)
  R_result <- cda:::.euler(c(pi/3,pi/3,pi/3))
  expect_equal(cpp_result, R_result)
})


## ------------------
## define constants
## ------------------

wavelength <- 500
medium <- 1.33
kn <- 2*pi/wavelength*medium

E0L=1/sqrt(2)*c(0,1,1i)
E0R=1/sqrt(2)*c(0,1i,1)
k0=c(1,0,0)

cluster <- list(r = rbind(c(0, 0, 0),
                          c(0, 0, 400)),
                angles = rbind(c(0, 0, 0),
                               c(pi/6, pi/2, 0)),
                sizes = rbind(c(40, 20, 20),
                              c(40, 20, 20)))
N <- nrow(cluster$r)

Beta <- inverse_polarizability(cluster, material=epsAu(wavelength), medium=medium)

test_that("inverse polarizabilities haven't changed", {  
  .Beta <- structure(c(-2.03834438665386e-05-3.67008788359938e-05i, 2.69269434775055e-05-3.67008788359939e-05i, 
                       2.69269434775055e-05-3.67008788359939e-05i, -2.03834438665386e-05-3.67008788359938e-05i, 
                       2.69269434775055e-05-3.67008788359939e-05i, 2.69269434775055e-05-3.67008788359939e-05i
  ), .Dim = c(6L, 1L), .Dimnames = list(c("aa", "ab", "ac", "aa", 
                                          "ab", "ac"), NULL))
  
  expect_equal(Beta, .Beta)
})

Angles <- cbind(c(0, pi/2, 0), # +x is phi=0, psi=0
                c(pi/2, pi/2, 0), # +y is phi=pi/2, psi=0
                c(pi/2, pi/2, pi/2)) # +z is phi=pi/2, psi=pi/2

## ------------------

test_that("Interaction matrix is correct", {  
  R_A <- cda:::.interaction_matrix(cluster$r, kn, Beta, cluster$angles)
  cpp_A <- cda$interaction_matrix(cluster$r, kn, c(Beta), cluster$angles, TRUE)
  
  expect_equal(cpp_A, R_A)
})

A <- cda$interaction_matrix(cluster$r, kn, c(Beta), cluster$angles, TRUE)

test_that("Diagonal blocks of the interaction matrix are correct", {  
  R_Adiag <- cda:::.block_diagonal(Beta, cluster$angles)
  cpp_Adiag <- cda$block_diagonal(c(Beta), cluster$angles)
  
  expect_equal(cpp_Adiag, R_Adiag)
})

Adiag <- cda$block_diagonal(c(Beta), cluster$angles)

test_that("Incident field at multiple angles is correct", {  
  R_Ei <- cda:::.incident_field(E0L, k=kn*k0, r=cluster$r, Angles)
  cpp_Ei <- cda$incident_field(E0L, k=kn*k0, r=cluster$r, Angles)
  
  expect_equal(cpp_Ei, R_Ei)
})


Ei <- cda$incident_field(E0L, k=kn*k0, r=cluster$r, Angles)
P <- solve(A, Ei)

test_that("Cross-sections are correct", {  
  
  cpp_extinction <- as.vector(cda$extinction(kn, P, Ei))
  R_extinction <- cda:::.extinction(kn, P, Ei)  
  cpp_absorption <- as.vector(cda$absorption(kn, P, Adiag))
  R_absorption <- cda:::.absorption(kn, P, Adiag)
  
  expect_equal(cpp_extinction, R_extinction)
  expect_equal(cpp_absorption, R_absorption)
})

