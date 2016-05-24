library(attribrisk)
library(testthat)

temp <- with(chapter.dat, table(cases, hbp))
whisnant <- data.frame(hbp=c(0,0,1,1), 
                       infarct= c(0,1,0,1),
                       y = c(temp))

fit1 <- attribrisk(infarct ~ expos(factor(hbp)), weight=y, data=whisnant )

test_that("attribrisk called using factor indicating dichotomous exposure.", {
  expect_that(is.null(fit1), equals(FALSE))
  expect_that(fit1$attribrisk, equals(0.313059, tolerance=1e-6))
  expect_that(fit1$var, equals(0.001361468, tolerance=1e-6))
})

fit2 <- attribrisk(infarct ~ expos(hbp),  weight=y, data=whisnant)
test_that("attribrisk called using numeric dichotomous exposure.", {
  expect_that(is.null(fit2), equals(FALSE))
  expect_that(fit2$attribrisk, equals(0.313059, tolerance=1e-6))
  expect_that(fit1$var, equals(0.001361468, tolerance=1e-6))
})


fit3 <- attribrisk(infarct ~ expos(1-hbp), weight=y, data=whisnant,baseline=list(hbp=0))
test_that("attribrisk called using numeric dichotomous exposure with specified baseline.", {
  expect_that(is.null(fit3), equals(FALSE))
  expect_that(fit3$attribrisk, equals(0.313059, tolerance=1e-6))
  expect_that(fit1$var, equals(0.001361468, tolerance=1e-6))
})

fit4 <- attribrisk(infarct ~ expos(hbp), weight=y, data=whisnant, var="none")
test_that("attribrisk called with var='none'.", {
  expect_that(is.null(fit4), equals(FALSE))
  expect_that(fit4$attribrisk, equals(0.313059, tolerance=1e-6))
  expect_that(is.null(fit4$var), equals(TRUE))
})

whisnant2 <- list()
for (i in 1:nrow(whisnant)){
  whisnant2[[i]] <- whisnant[ rep(i,whisnant[i,'y']),]
}
whisnant2 <- do.call(rbind, whisnant2)

fit7 <- attribrisk(infarct ~ expos(hbp), weight=y, data=whisnant, model=TRUE, x=TRUE, y=TRUE)
test_that("attribrisk called with model=TRUE, x=TRUE, y=TRUE.", {
  expect_that(is.null(fit7$model), equals(FALSE))
  expect_that(is.null(fit7$x), equals(FALSE))
  expect_that(is.null(fit7$y), equals(FALSE))
})
