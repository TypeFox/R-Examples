# unit test of validating functions

context("Validating edge lists")

test_that("Valid edgelist passes validation", {
  df1 <- data.frame( ego = 1:5, alter = c(2,3,2,5,4))
  r1 <- intergraph:::validateEL(df1)
  expect_equal(df1, r1)
} )

test_that("Edgelist with single column throws error", {
  expect_error( intergraph:::validateEL( data.frame(ego=1:5)))
} )

test_that("Edgelist with NA (invalid) gives warning", {
  df2 <- data.frame( ego = 1:5, alter = c(2,NA,2,5,4))
  expect_warning( intergraph:::validateEL(df2) )
} )






context("Validating vertex data frames")


test_that("Valid vertex DB passes validation", {
  df1 <- data.frame(id=1:5, x=c(1,2,3,2,1), y=c(3,2,1,2,3))
  expect_equal(df1, validateVDB(df1))
} )


test_that("Empty vertex DB throws error", {
  df2 <- data.frame(id=numeric(0), x=character(0))
  expect_error( validateVDB(df2) )
} )

test_that("Vertex DB with duplicated ids throws error", {
  df3 <- data.frame(id=1:5, x=c(1,2,3,2,1), y=c(3,2,1,2,3))
  df3$id[3] <- 1
  expect_error( validateVDB(df3))
} )

test_that("NAs in vertex ids gives warning", {
  df4 <- data.frame(id=1:5, x=c(1,2,3,2,1), y=c(3,2,1,2,3))
  df4$id[2] <- NA
  expect_warning( validateVDB(df4))
} )








context("Validating edgelist vs vertex data frame")

test_that("Valid edgelist and vertex DB pass validation", {
  edb <- data.frame( ego = 1:5, alter = c(2,3,2,5,4))
  vdb <- data.frame(id=1:5, x=c(1,2,3,2,1), y=c(3,2,1,2,3))
  expect_true( validNetDB(edb, vdb))
} )


test_that("Edgelist with some ids not present in VDB throw error", {
  elist <- data.frame( ego = 1:5, alter = c(2,3,2,5,4))
  vdb <- data.frame(id=1:4, x=c(1,2,3,4), y=c(3,2,1,2))
  expect_error( validNetDB(elist, vdb))
} )
