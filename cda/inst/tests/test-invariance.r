context("Checking the invariance upon spatial operations")

gold <- epsAu(seq(400, 600, by=100))
original <- cluster_helix(2, R0=12, pitch=15, 
                      delta=pi/2, delta0=0, right=TRUE,
                      a=5, b=5/2, c=5/2,
                      angles="helix")
## copies
swapped <- shifted <- original

.CD <- circular_dichroism_spectrum(original, gold)

rotated <- cluster_helix(2, R0=12, pitch=15, 
                         delta=pi/2, delta0=pi/3, right=TRUE,
                         a=5, b=5/2, c=5/2,
                         angles="helix")

swapped$r <- original$r[2:1, ]
swapped$angles <- original$angles[2:1, ]

shifted$r <- original$r + 1000

test_that("test invariance by swapping elements of the cluster", {
  CD <- circular_dichroism_spectrum(swapped, gold)
  expect_equal(CD, .CD)
})
test_that("test invariance by translation of the cluster", {
  CD <- circular_dichroism_spectrum(shifted, gold)
  expect_equal(CD, .CD)
})
test_that("test invariance by rotation of the cluster", {
  CD <- circular_dichroism_spectrum(rotated, gold)
  expect_equal(CD, .CD)
})


