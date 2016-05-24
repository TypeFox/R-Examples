
test_that("fxcat() returns the correct values", {
  vals <- c(NA, -3, 0, 1, 2.5, 5, 10)
  allNA <- c(fxcat(expand.grid(nfx=vals, ntot=vals[1:3])),
    fxcat(expand.grid(nfx=vals[1:2], ntot=vals)),
    fxcat(data.frame(nfx=vals[5:7], ntot=vals[5:7]-1)))
  all0 <- fxcat(data.frame(nfx=0, ntot=vals[4:7]))
  all50 <- fxcat(data.frame(nfx=vals[4:7], ntot=2*vals[4:7]))
  all100 <- fxcat(data.frame(nfx=vals[4:7], ntot=vals[4:7]))
  expect_that(all(is.na(allNA)), is_true())
  expect_that(all(all0==0), is_true())
  expect_that(all(all50==50), is_true())
  expect_that(all(all100==100), is_true())
})
