context("ODBCResult")
#Common functions
make_test_connection <- function()
{
  USER <- 'sa'
  PASSWORD <- 'Password12!'
  dbConnect(RODBCDBI::ODBC(), dsn='test', user=USER, password=PASSWORD)
}

test_that("DB source name should be test", {
  con <- make_test_connection()
  dbWriteTable(con, "iris", iris, overwrite=TRUE)
  res <- dbSendQuery(con, "SELECT * FROM iris")
  expect_true(dbGetInfo(res)$sourcename=="test")
  dbRemoveTable(con, "iris")
  dbDisconnect(con)
})

test_that("All rows and columns should be returned", {
  con <- make_test_connection()
  dbWriteTable(con, "iris", iris, overwrite=TRUE, rownames=FALSE)
  res <- dbSendQuery(con, "SELECT * FROM iris")
  df <- dbFetch(res)
  expect_true(nrow(df)==nrow(iris))
  expect_true(ncol(df)==ncol(iris))
  dbRemoveTable(con, "iris")
  dbDisconnect(con)
})

test_that("dbColumnInfo should give the information about its result", {
  con <- make_test_connection()
  dbWriteTable(con, "iris", iris, overwrite=TRUE, rownames=FALSE)
  res <- dbSendQuery(con, "SELECT * FROM iris")
  column_info <- dbColumnInfo(res)
  expect_true(all(sapply(iris, class)==column_info$data.type))
  expect_true(all(gsub("\\.", "", colnames(iris))==column_info$name))
  dbRemoveTable(con, "iris")
  dbDisconnect(con)
})

test_that("dbGetRowCount function should give the ", {
  con <- make_test_connection()
  dbWriteTable(con, "iris", iris, overwrite=TRUE, rownames=FALSE)
  res <- dbSendQuery(con, "SELECT * FROM iris")
  size <- dbGetRowCount(res)
  expect_true(size == nrow(iris))
  dbRemoveTable(con, "iris")
  dbDisconnect(con)
})

test_that("dbFetch should return the result depending on n argument", {
  con <- make_test_connection()
  dbWriteTable(con, "iris", iris, overwrite=TRUE, rownames=FALSE)
  res <- dbSendQuery(con, "SELECT * FROM iris")
  df <- dbFetch(res, n=3)
  expect_true(nrow(df)==3)
  # -1 means all row
  df <- dbFetch(res, n=-1)
  expect_true(nrow(df)==nrow(iris))
  dbRemoveTable(con, "iris")
  dbDisconnect(con)
})
