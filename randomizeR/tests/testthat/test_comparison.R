###################################################################
# ----------------------------------------------------------------#
# Tests for the compare function                                  #
# ----------------------------------------------------------------#
###################################################################

context("Comparison")

test_that("compare returns valid object", {
  RP1 <- getAllSeq(pbrPar(bc = c(2, 4, 2)))
  RP2 <- getAllSeq(bsdPar(10, 2))
  
  type <- sample(c("CS", "DS"), 1)
  i1 <- corGuess(type)
  method <- sample(c("sim", "exact"), 1)
  i2 <- selBias(type, 1, method)
  endp <- normEndp(c(0, 0), c(1, 1))
  expect_is(randomizeR::compare(i1, RP1), "comparison")
  expect_is(randomizeR::compare(i2, RP1, RP2, endp = endp), "comparison")
  
  type <- sample(c("imb", "absImb", "loss", "maxImb"), 1)
  i3 <- imbal(type)
  expect_is(randomizeR::compare(i3, RP2), "comparison")
  
  expect_error(randomizeR::compare(RP2))
  expect_error(randomizeR::compare(RP1, "blubb"))
  expect_error(randomizeR::compare(RP3, "issue"))
  expect_error(randomizeR::compare(RP1, 42))
  expect_error(randomizeR::compare(i2, RP2))
})

test_that("issue is computed right", {
  RP1 <- getAllSeq(pbrPar(bc = c(2, 2, 2)))
  RP2 <- getAllSeq(bsdPar(10, 1))
  RP3 <- genSeq(rpbrPar(N = 12, rb = 2), r = 100)
  RP4 <- getAllSeq(mpPar(16,1))
  i1 <- corGuess("CS")
  i2 <- corGuess("DS")
  i3 <- imbal("imb")
  i4<- imbal("absImb")
  i5 <- imbal("loss")
  i6 <- imbal("maxImb")
  
  # test for average value and standard deviation of the convergence strategy
  expect_equal(as.numeric(randomizeR::compare(i1, RP1, RP2, RP3, RP4)@S[1,]),
               rep(0.75, 4))
  expect_equal(as.numeric(randomizeR::compare(i1, RP1, RP2, RP3, RP4)@S[2,]),
               rep(0, 4))
  # test for average value and standard deviation of the divergence strategy
  expect_equal(as.numeric(randomizeR::compare(i2, RP1, RP2, RP3, RP4)@S[1,]),
               rep(0.25, 4))
  expect_equal(as.numeric(randomizeR::compare(i2, RP1, RP2, RP3, RP4)@S[2,]),
               rep(0, 4))
  # test for imbalance at the end
  expect_equal(as.numeric(randomizeR::compare(i3, RP1, RP2, RP3, RP4)@S[1,]),
               rep(0, 4))
  expect_equal(as.numeric(randomizeR::compare(i3, RP1, RP2, RP3, RP4)@S[2,]),
               rep(0, 4))
  # test for the absolute imbalance at the end
  expect_equal(as.numeric(randomizeR::compare(i4, RP1, RP2, RP3, RP4)@S[1,]),
               rep(0, 4))
  expect_equal(as.numeric(randomizeR::compare(i4, RP1, RP2, RP3, RP4)@S[2,]),
               rep(0, 4))
  # test for the loss
  expect_equal(as.numeric(randomizeR::compare(i5, RP1, RP2, RP3, RP4)@S[1,]),
               rep(0, 4))
  expect_equal(as.numeric(randomizeR::compare(i5, RP1, RP2, RP3, RP4)@S[2,]),
               rep(0, 4))
  # test for the maxImb
  expect_equal(as.numeric(randomizeR::compare(i6, RP1, RP2, RP3, RP4)@S[1,]),
               rep(1, 4))
  expect_equal(as.numeric(randomizeR::compare(i6, RP1, RP2, RP3, RP4)@S[2,]),
               rep(0, 4))
})