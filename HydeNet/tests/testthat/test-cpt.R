context("cpt")

test_that("cpt.list",
{
  n <- 50000
  df <- data.frame(
    di1 = as.factor(1:6 %*% rmultinom(n,1,prob=c(.4,.3,.15,.10,.03,.02))),
    di2 = as.factor(1:6 %*% rmultinom(n,1,prob=rev(c(.4,.3,.15,.10,.03,.02)))),
    di3 = as.factor(1:6 %*% rmultinom(n,1,prob=c(.15,.10,.02,.3,.4,.03)))
  )
  
  expect_that(cpt(list(y = "di3", x = c("di1", "di2")), data= df),
              not(throws_error()))
})

test_that("cpt with weights",
{
  echodata <- cbind(expand.grid(list(echo = c("Negative", "Positive"),
                                     cad = c("No","Yes"))),
                    data.frame(pr=c(0.83,0.17,0.12,0.88)))
  expect_that(cpt(echo ~ cad, data=echodata, wt=echodata$pr),
              not(throws_error()))
})

test_that("cpt with character string naming weights",
{
  echodata <- cbind(expand.grid(list(echo = c("Negative", "Positive"),
                                     cad = c("No","Yes"))),
                    data.frame(pr=c(0.83,0.17,0.12,0.88)))
  expect_that(cpt(echo ~ cad, echodata, wt="pr"),
              not(throws_error()))
})

test_that("cpt with multiple character string naming weights",
{
  echodata <- cbind(expand.grid(list(echo = c("Negative", "Positive"),
                                     cad = c("No","Yes"))),
                    data.frame(pr=c(0.83,0.17,0.12,0.88)))
  expect_warning(cpt(echo ~ cad, echodata, wt=c("pr", "wt")))
})

test_that("cpt with multiple character string naming weights, but choosing non-existent variable",
{
  echodata <- cbind(expand.grid(list(echo = c("Negative", "Positive"),
                                     cad = c("No","Yes"))),
                    data.frame(pr=c(0.83,0.17,0.12,0.88)))
  expect_error(cpt(echo ~ cad, echodata, wt=c("wt", "pr")))
})

test_that("cpt with logical weights should return error",
{
  echodata <- cbind(expand.grid(list(echo = c("Negative", "Positive"),
                                     cad = c("No","Yes"))),
                    data.frame(pr=c(0.83,0.17,0.12,0.88)))
  expect_error(cpt(echo ~ cad, echodata, wt=c(TRUE, TRUE, FALSE, TRUE)))
})

test_that("cpt with inappropriate length for weights vector",
{
  echodata <- cbind(expand.grid(list(echo = c("Negative", "Positive"),
                                     cad = c("No","Yes"))),
                    data.frame(pr=c(0.83,0.17,0.12,0.88)))
  expect_error(cpt(echo ~ cad, echodata, wt=c(1, 2, 3)))
})

test_that("cpt with negative weights should cast error",
{
  echodata <- cbind(expand.grid(list(echo = c("Negative", "Positive"),
                                     cad = c("No","Yes"))),
                    data.frame(pr=c(0.83,0.17,0.12,0.88)))
  expect_error(cpt(echo ~ cad, echodata, wt=c(-1, 1, 1, 1)))
})

test_that("print.cpt succeeds",
{
  echodata <- cbind(expand.grid(list(echo = c("Negative", "Positive"),
                                     cad = c("No","Yes"))),
                    data.frame(pr=c(0.83,0.17,0.12,0.88)))
  x <- cpt(echo ~ cad, echodata, wt="pr")
  expect_that(print(x), not(throws_error()))
})
  