library(testthat)
library(cda)
library(plyr)
context("Checking the intermediate quantities in CD calculations")

wavelength <- c(500, 550, 600)
material <- epsAu(wavelength)
medium <- 1.3  

k0 <- 2*pi/wavelength
kn <- k0*medium

cluster <- list(r = rbind(c(0, 0, 0),
                     c(0, 0, 200)),
           angles = rbind(c(0, 0, 0),
                          c(pi/4, 0, 0)),
           sizes = rbind(c(40, 20, 20),
                         c(40, 20, 20)))

Beta <- inverse_polarizability(cluster, material, 
                                  polarizability_fun=polarizability_ellipsoid, 
                                  medium=medium, kuwata=TRUE)

.Beta <-structure(c(-1.82237397053259e-05-3.56779387878792e-05i, 2.89922059965679e-05-3.56779387878792e-05i, 
                        2.89922059965679e-05-3.56779387878792e-05i, -1.82237397053259e-05-3.56779387878792e-05i, 
                        2.89922059965679e-05-3.56779387878792e-05i, 2.89922059965679e-05-3.56779387878792e-05i, 
                        -1.40064964665508e-05-1.3699744599798e-05i, 3.2848949238791e-05-1.3699744599798e-05i, 
                        3.2848949238791e-05-1.3699744599798e-05i, -1.40064964665508e-05-1.3699744599798e-05i, 
                        3.2848949238791e-05-1.3699744599798e-05i, 3.2848949238791e-05-1.3699744599798e-05i, 
                        -2.39972560681825e-06-5.46975933994133e-06i, 4.41715911001179e-05-5.4697593399413e-06i, 
                        4.41715911001179e-05-5.4697593399413e-06i, -2.39972560681825e-06-5.46975933994133e-06i, 
                        4.41715911001179e-05-5.4697593399413e-06i, 4.41715911001179e-05-5.4697593399413e-06i
), .Dim = c(6L, 3L), .Dimnames = list(c("aa", "ab", "ac", "aa", 
                                        "ab", "ac"), NULL))

# [,1]                        [,2]                        [,3]
# aa -1.822374e-05-3.567794e-05i -1.400650e-05-1.369974e-05i -2.399726e-06-5.469759e-06i
# ab  2.899221e-05-3.567794e-05i  3.284895e-05-1.369974e-05i  4.417159e-05-5.469759e-06i
# ac  2.899221e-05-3.567794e-05i  3.284895e-05-1.369974e-05i  4.417159e-05-5.469759e-06i
# aa -1.822374e-05-3.567794e-05i -1.400650e-05-1.369974e-05i -2.399726e-06-5.469759e-06i
# ab  2.899221e-05-3.567794e-05i  3.284895e-05-1.369974e-05i  4.417159e-05-5.469759e-06i
# ac  2.899221e-05-3.567794e-05i  3.284895e-05-1.369974e-05i  4.417159e-05-5.469759e-06i

test_that("the inverse polarisability is the same as earlier versions", {
  expect_equal(Beta, .Beta)
})



A <- cda$interaction_matrix(cluster$r, kn[1], Beta, cluster$angles, TRUE)
.A <- structure(c(-1.82237397053259e-05-3.56779387878792e-05i, 0+0i, 
                  0+0i, 1.14864729173871e-06+5.5676095065754e-07i, 0+0i, 0+0i, 
                  0+0i, 2.89922059965679e-05-3.56779387878792e-05i, 0+0i, 0+0i, 
                  1.14864729173871e-06+5.5676095065754e-07i, 0+0i, 0+0i, 0+0i, 
                  2.89922059965679e-05-3.56779387878792e-05i, 0+0i, 0+0i, 3.5040262644085e-07-7.79039958472603e-07i, 
                  1.14864729173871e-06+5.5676095065754e-07i, 0+0i, 0+0i, 5.384233145621e-06-3.56779387878792e-05i, 
                  -0.0000236079728509469+0i, 0+0i, 0+0i, 1.14864729173871e-06+5.5676095065754e-07i, 
                  0+0i, -0.0000236079728509469+0i, 5.384233145621e-06-3.56779387878792e-05i, 
                  0+0i, 0+0i, 0+0i, 3.5040262644085e-07-7.79039958472603e-07i, 
                  0+0i, 0+0i, 2.89922059965679e-05-3.56779387878792e-05i), .Dim = c(6L, 
                                                                                    6L))

test_that("the interaction matrix is the same as earlier versions", {
  expect_equal(A, .A)
})

.Adiag <- structure(c(-1.82237397053259e-05-3.56779387878792e-05i, 0+0i, 
                      0+0i, 0+0i, 0+0i, 0+0i, 0+0i, 2.89922059965679e-05-3.56779387878792e-05i, 
                      0+0i, 0+0i, 0+0i, 0+0i, 0+0i, 0+0i, 2.89922059965679e-05-3.56779387878792e-05i, 
                      0+0i, 0+0i, 0+0i, 0+0i, 0+0i, 0+0i, 5.384233145621e-06-3.56779387878792e-05i, 
                      -0.0000236079728509469+0i, 0+0i, 0+0i, 0+0i, 0+0i, -0.0000236079728509469+0i, 
                      5.384233145621e-06-3.56779387878792e-05i, 0+0i, 0+0i, 0+0i, 0+0i, 
                      0+0i, 0+0i, 2.89922059965679e-05-3.56779387878792e-05i), .Dim = c(6L, 
                                                                                        6L))
Adiag <- cda$block_diagonal(c(Beta), cluster$angles)

test_that("diagonal of the interaction matrix is the same as earlier versions", {
  expect_equal(Adiag, .Adiag)
  expect_equal(diag(Adiag), diag(A))
})

## cheap averaging

angles <- rbind(c(0, pi/2, 0), # +x is phi=0, psi=0
                c(pi/2, pi/2, 0), # +y is phi=pi/2, psi=0
                c(pi/2, pi/2, pi/2)) # +z is phi=pi/2, psi=pi/2
weights <- rep(1/nrow(angles), nrow(angles)) 


## single wavelength
res <- cd$averaging(cluster$r, A, Adiag, kn[1],
                           as.matrix(angles), weights)

.res <- structure(c(7629.60523941392, 7629.82209301281, 6966.45784786352, 
                    6963.78305385718), .Dim = c(4L, 1L))


## spectrum
resall <- cd$average_spectrum(kn, Beta, cluster$r, cluster$angles, 
                                       as.matrix(angles), weights,
                                       TRUE, FALSE)

.resall <- structure(c(7629.71366621337, 7180.9666774582, 18093.6538420437, 
                    6965.12045086035, 5975.47136711237, 12499.019321288, 664.593215353016, 
                    1205.49531034583, 5594.63452075571, -0.216853598899434, -26.9159455604804, 
                    -533.349596245193, 2.67479400633692, -30.6246464849028, -361.940147617594, 
                    -2.89164760523636, 3.70870092442237, -171.409448627599), .Dim = c(3L, 6L))


test_that("the result of cheap orientation averaging is the same as earlier versions", {
  expect_equal(res,.res)  
  expect_equal(resall,.resall)
})


