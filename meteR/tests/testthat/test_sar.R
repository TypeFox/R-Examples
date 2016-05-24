context('SAR works')
## values of S0, N0, A0 and predicted S from Newman et al. 


test_that('predicted SAR values are correct', {
  bell.sar <- meteSAR(S0=32, N0=920, Amin=1, A0=16)$pred$S
  newmanS <- c(12.2560715, 15.8991283, 20.2158684, 25.7856639, 32)
  expect_true(all(round(bell.sar, 3) == round(newmanS, 3)))
})
