context("verificaTipoColuna return")

test_that("verificaTipoColuna", {
  ID_plant <- c("ABC","ABC","ABC","ABC")
  ID_Re <- c(1,2,3,4)
  test <- data.frame(ID_Re,ID_plant)
  expect_equal(verificaTipoColuna(test$ID_Re), "as.numeric()")
})
