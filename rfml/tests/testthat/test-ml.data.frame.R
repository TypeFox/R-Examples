context("ml.data.frame")

test_that("can create and delete a ml.data.frame based on iris dataset using json format", {
   skip_on_cran()
   myConn <- ml.connect(port = "8088")
   mlIris <- as.ml.data.frame(myConn, iris, "iris-test", format = "json")
   expect_is(mlIris, "ml.data.frame")
   expect_true(is.ml.data.frame(mlIris))
   expect_equal(mlIris@.nrows, 150)
   expect_true(rm.ml.data.frame(mlIris))
 })

test_that("can create and delete a ml.data.frame based on iris dataset using xml format", {
  skip_on_cran()
  myConn <- ml.connect(port = "8088")
  mlIris <- as.ml.data.frame(myConn, iris, "iris-test", format = "XML")
  expect_is(mlIris, "ml.data.frame")
  expect_true(is.ml.data.frame(mlIris))
  expect_equal(mlIris@.nrows, 150)
  expect_true(rm.ml.data.frame(mlIris))
})

test_that("can create a ml.data.frame based on search", {
  skip_on_cran()
  myConn <- ml.connect(port = "8088")
   mlIris <- as.ml.data.frame(myConn, iris, "iris-test")
   mlIris2 <- ml.data.frame(myConn, collection = "iris-test")
   expect_is(mlIris2, "ml.data.frame")
   expect_true(is.ml.data.frame(mlIris2))
   expect_equal(mlIris2@.nrows, 150)
   expect_equal(nrow(mlIris2), 150)
   mlIris3 <- ml.data.frame(myConn, query = "setosa", collection = "iris-test")
   expect_equal(mlIris3@.nrows, 50)
   expect_equal(nrow(mlIris3), 50)
   mlIris4 <- ml.data.frame(myConn, query = "virginica", collection = "iris-test", directory = "/rfml/admin/iris-test/")
   expect_equal(mlIris4@.nrows, 50)
   expect_equal(nrow(mlIris4), 50)
   rm.ml.data.frame(mlIris)
})

test_that("can create a ml.data.frame based on fieldQuery", {
  skip_on_cran()
  myConn <- ml.connect(port = "8088")
  mlIris <- as.ml.data.frame(myConn, iris, "iris-test")
  mlIris2 <- ml.data.frame(myConn, fieldFilter = "Species == virginica", collection = "iris-test")
  expect_equal(mlIris2@.nrows, 50)
  expect_equal(nrow(mlIris2), 50)
  db <- "rfml"
  expect_message(ml.add.index(x = mlIris$Petal.Length, scalarType = "decimal", database =  db, conn = myConn), "Range element index created on Petal.Length")
  # We need to wait so that the index gets updated before using a function that leverage it
  Sys.sleep(10)
  mlIris3 <- ml.data.frame(myConn, fieldFilter = "Petal.Length > 4.5", collection = "iris-test")
  expect_equal(mlIris3@.nrows, 63)
  expect_equal(nrow(mlIris3), 63)
  mlIris4 <- ml.data.frame(myConn, query = "virginica", fieldFilter = "Petal.Length > 5", collection = "iris-test", directory = "/rfml/admin/iris-test/")
  expect_equal(mlIris4@.nrows, 41)
  expect_equal(nrow(mlIris4), 41)
  rm.ml.data.frame(mlIris)
})

test_that("can create new fields on a ml.data.frame", {
  skip_on_cran()
  myConn <- ml.connect(port = "8088")
  mlIris <- as.ml.data.frame(myConn,iris, "iris-test")
  mlIris$SepLength <- mlIris$Sepal.Length
  expect_is(mlIris$SepLength, "ml.col.def")
  expect_match(mlIris$SepLength@.expr, "rfmlResult['Sepal.Length']", fixed=TRUE)
  expect_equal(length(mlIris@.col.defs), 1)
  expect_equal(length(mlIris@.col.name), 6)
  mlIris$SepLength10 <- mlIris$Sepal.Length * 10
  expect_output(mlIris$SepLength10@.expr, "(rfmlResult['Sepal.Length']*10)", fixed=TRUE)
  expect_equal(length(mlIris@.col.defs), 2)
  expect_equal(length(mlIris@.col.name), 7)
  mlIris$SepRatio <- mlIris$Sepal.Length / mlIris$Sepal.Width
  expect_output(mlIris$SepRatio@.expr, "(rfmlResult['Sepal.Length']/rfmlResult['Sepal.Width'])", fixed=TRUE)
  expect_equal(length(mlIris@.col.defs), 3)
  expect_equal(length(mlIris@.col.name), 8)
  mlIris$SepLengthAbs <- abs(mlIris$Sepal.Length)
  expect_output(mlIris$SepLengthAbs@.expr, "fn.abs(rfmlResult['Sepal.Length'])", fixed=TRUE)
  expect_equal(length(mlIris@.col.defs), 4)
  expect_equal(length(mlIris@.col.name), 9)
  rm.ml.data.frame(mlIris)
})
test_that("sub select on a ml.data.frame", {
  skip_on_cran()
  myConn <- ml.connect(port = "8088")
  mlIris <- as.ml.data.frame(myConn, iris, "iris-test")
  mlIris2 <- mlIris[1:3]
  expect_equal(length(mlIris2@.col.name), 3)
  mlIris3 <- mlIris[,1:3]
  expect_equal(length(mlIris3@.col.name), 3)
  mlIris4 <- mlIris[,c("Sepal.Length","Sepal.Width","Petal.Length")]
  expect_equal(length(mlIris4@.col.name), 3)
  mlIris5 <- mlIris[mlIris$Species=="setosa", 1:3]
  expect_equal(nchar(mlIris5@.queryArgs$`rs:fieldQuery`), 99)
  expect_true(mlIris5@.extracted)
  expect_equal(length(mlIris5@.col.name), 3)
  mlIris6 <- mlIris[mlIris$Species=="setosa",]
  expect_equal(nchar(mlIris6@.queryArgs$`rs:fieldQuery`), 99)
  expect_equal(length(mlIris6@.col.name), 5)
  expect_false(mlIris6@.extracted)
  mlIris7 <- mlIris["setosa",]
  expect_output(mlIris7@.queryArgs$`rs:q`, "setosa")
  expect_equal(length(mlIris7@.col.name), 5)
  rm.ml.data.frame(mlIris)
})

test_that("can create data based on a ml.data.frame", {
  skip_on_cran()
  myConn <- ml.connect(port = "8088")
  mlIris <- as.ml.data.frame(myConn, iris, "iris-test")
  mlIris$SepLength <- mlIris$Sepal.Length
  mlIris$SepLength10 <- mlIris$Sepal.Length * 10
  mlIris$SepRatio <- mlIris$Sepal.Length / mlIris$Sepal.Width
  mlIris$SepLengthAbs <- abs(mlIris$Sepal.Length)
  newIris <- as.ml.data.frame(x = mlIris, name = "newIris-test" )
  expect_equal(nrow(newIris), 150)
  expect_equal(length(newIris@.col.name), 9)
  expect_equal(length(newIris@.col.defs), 0)
  rm.ml.data.frame(mlIris)
  rm.ml.data.frame(newIris)
})


