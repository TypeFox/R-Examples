require(testthat)
require(lfactors)

mon <- lfactor(1:12,
               levels=1:12,
               labels=c("Jan", "Feb", "Mar", "Apr", "May","Jun",
                        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

context("equal")
test_that("equal", {
  expect_equal(mon == "Feb", c(FALSE,TRUE,rep(FALSE,10)))
  expect_equal(mon == "Feb", mon==2)
  expect_equal(mon[3] == c("Jan", "Feb", "Mar"), mon[3] == 1:3)
  expect_equal(mon[1:2] == c("Feb", "Tuesday"), mon[1:2] == c(2,-4) )
})

context("not equal")
test_that("not equal", {
  expect_equal(mon != "Feb", mon != 2)						
  expect_equal(mon[3] == c("Jan", "Feb", "Mar"), mon[3] == 1:3)
})

context("in")
test_that("in", {
  expect_equal(mon %in% c(2, 3), mon %in% c("Feb", "Mar"))	
  expect_equal(c(-4, 14,3,10) %in% mon, c("not a month", "Third December","Mar","Oct") %in% mon)
})

context("GT and GTE")
test_that("GT and GTE", {
  expect_equal(mon >  3, c(rep(FALSE,3), rep(TRUE,9)))
  expect_equal(mon >= 3, c(rep(FALSE,2), rep(TRUE,10)))
  expect_equal(3 < mon, c(rep(FALSE,3), rep(TRUE,9)))
  expect_equal(3 <= mon, c(rep(FALSE,2), rep(TRUE,10)))
  expect_equal(mon >  "Mar", c(rep(FALSE,3), rep(TRUE,9)))
  expect_equal(mon >= "Mar", c(rep(FALSE,2), rep(TRUE,10)))
  expect_equal("Mar" < mon, c(rep(FALSE,3), rep(TRUE,9)))
  expect_equal("Mar" <= mon, c(rep(FALSE,2), rep(TRUE,10)))
})

context("droplevels")
test_that("droplevels", {
  dl <- droplevels(x <- lfactor(c(1,3),levels=1:3, labels=LETTERS[1:3]))
  expect_is(dl, "lfactor") 
  expect_equal(attributes(dl)$levels, LETTERS[c(1,3)])
  expect_equal(attributes(dl)$llevels, c(1,3))
  dlprime <- lfactor(c(1,3),levels=c(1,3), labels=LETTERS[c(1,3)])
  expect_equal(dlprime, dl)
})

context("relevel")
test_that("relevel", {
  monp <- relevel(mon, "Jun")
  expect_equal(mon=="Jun", monp=="Jun")
  expect_equal(mon=="Jan", monp=="Jan")
  expect_equal(mon=="Feb", monp=="Feb")
  expect_equal(mon==6, monp==6)
  expect_equal(mon==1, monp==1)
  expect_equal(mon==2, monp==2)
})

context("set text")
test_that("set text", {
  monp <- mon
  monp[4] <- "Jun"
  expect_equal(monp=="Jun", monp==6)
  expect_equal(monp=="Jun", c("Jan", "Feb", "Mar", "Jun", "May","Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")=="Jun")
})

context("set num")
test_that("set num", {
  monp <- mon
  monp[4] <- 6
  expect_equal(monp=="Jun", monp==6)
  expect_equal(monp=="Jun", c("Jan", "Feb", "Mar", "Jun", "May","Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")=="Jun")
  monp <- mon
  monp[4:5] <- c(6,"Feb")
  expect_equal(monp=="Jun", monp==6)
  expect_equal(monp=="Feb", monp==2)
  expect_equal(monp=="Jun", c("Jan", "Feb", "Mar", "Jun", "Feb","Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")=="Jun")
  expect_equal(monp=="Feb", c("Jan", "Feb", "Mar", "Jun", "Feb","Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")=="Feb")

})

context("get num")
test_that("get num", {
  let <- lfactor(4:12,
                 levels=4:12,
                 labels=letters[4:12])
  expect_equal(as.numeric(let), 4:12)  
  expect_equal(as.double(let), 4:12)  
  skip_on_cran()
  expect_equal(as.integer(let), 1:9)
})

context("subset drop argument")
test_that("subset drop argument", {
  expect_equal(levels(mon[1:4,drop=TRUE]), c("Jan", "Feb", "Mar", "Apr"))
  expect_equal(levels(mon[1:4,drop=FALSE]), levels(mon))
})

context("subset with no i works")
test_that("subset with no i works", {
  monb <- mon[1:11]
  expect_equal(mon[1:11,drop=TRUE], monb[,drop=TRUE])
})

context("lm with a dropped level")
test_that("lm with a dropped level", {
  data1 <- data.frame(y=1:100, x=seq(1:10), z=lfactor(rep(1:5,each=20),1:5,letters[1:5]))
  data2 <- data.frame(y=1:100, x=seq(1:10), z=factor(rep(1:5,each=20),1:5,letters[1:5]))
  lm1 <- lm(y ~ x+z, data1)
  lm2 <- lm(y ~ x+z, data2)
  expect_equal(coef(summary(lm1)), coef(summary(lm2)))
  lmp1 <- lm(y ~ x+z, subset(data1, z %in% letters[1:4]))
  lmp2 <- lm(y ~ x+z, subset(data2, z %in% letters[1:4]))
  expect_equal(coef(summary(lmp1)), coef(summary(lmp2)))
})


context("sparse.model.matrix with dropped levels")
test_that("lm with a dropped level", {
  skip_on_cran()
  if(!exists("sparse.model.matrix")) {
    skip("Matrix package not loaded")
  }
  df0 <- data.frame(a=1:20, b=lfactor(sample(2:4, 20, replace=TRUE), 1:5, letters[1:5]))
  df0$d <- droplevels(df0$b)
  X1 <- model.matrix(a ~ b, df0)
  X1p <- sparse.model.matrix(a ~ b, df0)
  attributes(X1) <- NULL
  X1p <- as.matrix(X1p)
  attributes(X1p) <- NULL
  expect_equal(X1, X1p)

  X2 <- model.matrix(a ~ d, df0)
  X2p <- sparse.model.matrix(a ~ d, df0)
  attributes(X2) <- NULL
  X2p <- as.matrix(X2p)
  attributes(X2p) <- NULL
  expect_equal(X2, X2p)

  df0$b <- relevel(df0$b, "c")
  X3 <- model.matrix(a ~ b, df0)
  X3p <- sparse.model.matrix(a ~ b, df0)
  attributes(X3) <- NULL
  X3p <- as.matrix(X3p)
  attributes(X3p) <- NULL
  expect_equal(X3, X3p)
})

