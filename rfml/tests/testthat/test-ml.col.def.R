context("ml.col.def")

test_that("ml.col.def methods", {
  skip_on_cran()
  myConn<-ml.connect(port = "8088")
  mlIris <- as.ml.data.frame(myConn, iris, "iris-test")
  expect_equal(mlIris$Petal.Width@.name, "Petal.Width")
  rm.ml.data.frame(mlIris)
})

