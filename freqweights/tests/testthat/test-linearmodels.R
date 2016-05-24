## Emilio Torres Manzanera
## University of Oviedo
## Time-stamp: <2014-04-19 Sat 16:24 emilio on emilio-despacho>
## ============================================================

context("Linear model estimations")

test_that("xxxfreq obtain identical results than linear models",{
ml <- lm(Sepal.Length ~ Sepal.Width, iris)
mfq <- lmfreq(Sepal.Length ~ Sepal.Width, iris)
mbfq <- biglmfreq(Sepal.Length ~ Sepal.Width, iris)

expect_that(coef(mfq), equals(coef(ml)))
expect_that(coef(mbfq), equals(coef(ml)))
expect_that(AIC(mfq), equals(AIC(ml)))

chunk1 <- iris[1:30,]
chunk2 <- iris[-c(1:30),]
mf1 <- biglmfreq(Sepal.Length ~ Sepal.Width, chunk1)
mf2 <- update(mf1, chunk2)
expect_that(coef(mf2), equals(coef(ml)))

})
