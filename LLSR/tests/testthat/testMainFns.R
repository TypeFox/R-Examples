library(LLSR)


test_that("Mechuk's Coefficients", {
  expect_equal(summary(mrchk(peg4kslt))$coefficients[1], 0.94926887081, tolerance=1e-7)
  expect_equal(summary(mrchk(peg4kslt))$coefficients[2], -5.08275261773, tolerance=1e-7)
  expect_equal(summary(mrchk(peg4kslt))$coefficients[3], 787.55369595231, tolerance=1e-7)
})

test_that("Murugesan's Coefficients", {
  expect_equal(summary(mrgsn(peg4kslt))$coefficients[1], 0.90388691803, tolerance=1e-7)
  expect_equal(summary(mrgsn(peg4kslt))$coefficients[2], -3.48971940000, tolerance=1e-7)
  expect_equal(summary(mrgsn(peg4kslt))$coefficients[3], 2.92380736863, tolerance=1e-7)
})

# AQSys.gsnchk(peg4kslt[1:2],peg4kslt[2,3],peg4kslt[2,4],peg4kslt[2,5],
# peg4kslt[2,6],peg4kslt[2,7],peg4kslt[2,8])
