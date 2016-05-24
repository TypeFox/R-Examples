context("ODBCConnection")
#Common functions
make_test_connection <- function()
{
  USER <- 'sa'
  PASSWORD <- 'Password12!'
  dbConnect(RODBCDBI::ODBC(), dsn='test', user=USER, password=PASSWORD)
}

test_that("Connection should be established", {
  con <- make_test_connection()
  expect_true(!(is.null(con)))
  dbDisconnect(con)
})

test_that("Connection should not be established", {
  suppressWarnings(expect_error(dbConnect(RODBCDBI::ODBC(), dsn='test_nothing_connection')))
})


test_that("iris exists by dbExistsTable", {
  con <- make_test_connection()
  dbWriteTable(con, "iris", iris, overwrite=TRUE)
  expect_true(dbExistsTable(con, "iris"))
  dbRemoveTable(con, "iris")
  dbDisconnect(con)
})

test_that("iris exists by dbListTables", {
  con <- make_test_connection()
  dbWriteTable(con, "iris", iris, overwrite=TRUE)
  expect_true(tolower("iris") %in% tolower(dbListTables(con)))
  dbRemoveTable(con, "iris")
  dbDisconnect(con)
})

test_that("iris does not exists", {
  con <- make_test_connection()
  dbRemoveTable(con, "iris")
  expect_false(dbExistsTable(con, "iris"))
  dbDisconnect(con)
})

test_that("COlumn name should match", {
  con <- make_test_connection()
  dbWriteTable(con, "iris", iris, overwrite=TRUE, rownames=FALSE)
  expect_true(all(gsub("\\.", "", colnames(iris)) == dbListFields(con, "iris")))
  dbRemoveTable(con, "iris")
  dbDisconnect(con)
})

test_that("COlumn name should match including rownames", {
  con <- make_test_connection()
  dbWriteTable(con, "iris", iris, overwrite=TRUE)
  expect_true(all(c("rownames", gsub("\\.", "", colnames(iris))) == dbListFields(con, "iris")))
  dbRemoveTable(con, "iris")
  dbDisconnect(con)
})


test_that("it should not behave like bringing any error", {
  con <- make_test_connection()
  dbDisconnect(con)
  dbDisconnect(con)
  dbDisconnect(con)
})

test_that("iris table and raw iris data should be match", {
  con <- make_test_connection()
  dbWriteTable(con, "iris", iris, overwrite=TRUE, rownames=FALSE)
  suppressWarnings(expect_true(all(dbReadTable(con, "iris") == iris)))
  dbRemoveTable(con, "iris")
  dbDisconnect(con)
})

test_that("DB source name should be test", {
  con <- make_test_connection()
  expect_true(dbGetInfo(con)$sourcename=="test")
  dbDisconnect(con)
})
