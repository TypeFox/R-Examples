# Test the FocalPointQuery SQL builder
# Author: David Barker <support@causata.com>
#
library(testthat)
library(Causata)

equals <- testthat::equals

source("utils.R")

context("FocalPointQuery")

clean <- function(...) {
  strip(paste(...))
}

test_that("Simplest possible FocalPointQuery", {
  q <- FocalPointQuery("subscribe")

  expect_that(
    strip(as.character(q)),
    equals("SELECT * FROM Scenarios variable, `subscribe` WHERE variable.focal_point = `subscribe`.timestamp")
  )
})

test_that("FocalPointQuery on event attribute", {
  q <- FocalPointQuery("subscribe", event.attribute="pos-date")
  
  expect_that(
    strip(as.character(q)),
    equals("SELECT * FROM Scenarios variable, `subscribe` WHERE variable.focal_point = `subscribe`.`pos-date`")
  )
})

test_that("Extra events appear in the SELECT", {
  q <- FocalPointQuery("subscribe") + WithEvents("page-view", "purchase")
  
  expect_that(
    strip(as.character(q)),
    equals("SELECT * FROM Scenarios variable, `subscribe`, `page-view`, `purchase` WHERE variable.focal_point = `subscribe`.timestamp")
  )
})

test_that("Add variables", {
  original <- FocalPointQuery("subscribe")
  with.variables <- original + WithVariables("page-view-count", "total-spend")
  
  expect_that(strip(as.character(original)), equals(clean(
    "SELECT *",
    "FROM Scenarios variable, `subscribe`",
    "WHERE variable.focal_point = `subscribe`.timestamp"
  )))
  
  expect_that(strip(as.character(with.variables)), equals(clean(
    "SELECT `page-view-count`,`total-spend`",
    "FROM Scenarios variable, `subscribe`",
    "WHERE variable.focal_point = `subscribe`.timestamp"
  )))
})

test_that("Get and set variables", {
  original <- FocalPointQuery("subscribe")
  with.variables <- original + WithVariables("page-view-count", "total-spend")
  
  expect_that(Variables(original), equals(NULL))
  expect_that(Variables(with.variables), equals(c("page-view-count", "total-spend")))
  
  Variables(original) <- c("added-1", "added-2")
  expect_that(Variables(original), equals(c("added-1", "added-2")))
  expect_that(Variables(with.variables), equals(c("page-view-count", "total-spend")))
  
  expect_that(strip(as.character(original)), equals(clean(
    "SELECT `added-1`,`added-2`",
    "FROM Scenarios variable, `subscribe`",
    "WHERE variable.focal_point = `subscribe`.timestamp"
  )))
  
  expect_that(strip(as.character(with.variables)), equals(clean(
    "SELECT `page-view-count`,`total-spend`",
    "FROM Scenarios variable, `subscribe`",
    "WHERE variable.focal_point = `subscribe`.timestamp"
  )))
})

test_that("Add limit", {
  original <- FocalPointQuery("subscribe")
  with.limit <- original + Limit(1000)
  
  expect_that(strip(as.character(original)), equals(clean(
    "SELECT *",
    "FROM Scenarios variable, `subscribe`",
    "WHERE variable.focal_point = `subscribe`.timestamp"
  )))
  expect_that(strip(as.character(with.limit)), equals(clean(
    "SELECT *",
    "FROM Scenarios variable, `subscribe`",
    "WHERE variable.focal_point = `subscribe`.timestamp",
    "LIMIT 1000"
  )))
})

test_that("Get and set limit", {
  original <- FocalPointQuery("subscribe")
  with.limit <- original + Limit(1000)
  
  expect_that(Limit(original), equals(NULL))
  expect_that(Limit(with.limit), equals(1000))
  
  Limit(original) <- 123
  expect_that(Limit(original), equals(123))
  expect_that(Limit(with.limit), equals(1000))
    
  expect_that(strip(as.character(original)), equals(clean(
    "SELECT *",
    "FROM Scenarios variable, `subscribe`",
    "WHERE variable.focal_point = `subscribe`.timestamp",
    "LIMIT 123"
  )))
  expect_that(strip(as.character(with.limit)), equals(clean(
    "SELECT *",
    "FROM Scenarios variable, `subscribe`",
    "WHERE variable.focal_point = `subscribe`.timestamp",
    "LIMIT 1000"
  )))
})

test_that("Add WHERE clauses with functions", {
  q <- FocalPointQuery("subscribe") +
    WithVariables("page-view-count", "total-spend") +
    Where("total-spend", GreaterThan(100)) + 
    Where("page-view-count", GreaterThan(2))
  
  expect_that(
    strip(as.character(q)),
    equals(clean(
      "SELECT `page-view-count`,`total-spend`",
      "FROM Scenarios variable, `subscribe`",
      "WHERE variable.focal_point = `subscribe`.timestamp",
      "AND variable.`total-spend` > 100",
      "AND variable.`page-view-count` > 2"
    ))
  )
})

test_that("Add WHERE clauses with variable name", {
  q <- FocalPointQuery("subscribe") +
    WithVariables("page-view-count", "total-spend") +
    Where("total-spend", "<>123") + 
    Where("page-view-count", "IN (2, 3)")
  
  expect_that(
    strip(as.character(q)),
    equals(clean(
      "SELECT `page-view-count`,`total-spend`",
      "FROM Scenarios variable, `subscribe`",
      "WHERE variable.focal_point = `subscribe`.timestamp",
      "AND variable.`total-spend` <>123",
      "AND variable.`page-view-count` IN (2, 3)"
    ))
  )
})

test_that("Add raw WHERE clauses", {
  q <- FocalPointQuery("subscribe") +
       WithVariables("page-view-count", "total-spend") +
       Where("variable.`total-spend` > 100") + 
       Where("variable.`page-view-count` > 2")
  
  expect_that(
    strip(as.character(q)),
    equals(clean(
      "SELECT `page-view-count`,`total-spend`",
      "FROM Scenarios variable, `subscribe`",
      "WHERE variable.focal_point = `subscribe`.timestamp",
      "AND variable.`total-spend` > 100",
      "AND variable.`page-view-count` > 2"
    ))
  )
})

test_that("Oldest event can be specified", {
  q <- FocalPointQuery("event", "using.oldest.event")
  
  expect_that(
    strip(as.character(q)),
    equals(clean(
      "SELECT * FROM Scenarios variable, `event`",
      "WHERE variable.focal_point = `event`.timestamp",
      "AND IS_LAST(`event`.timestamp)"
    ))
  )
})

test_that("Newest event can be specified", {
  q <- FocalPointQuery("event", "using.newest.event")
  
  expect_that(
    strip(as.character(q)),
    equals(clean(
      "SELECT * FROM Scenarios variable, `event`",
      "WHERE variable.focal_point = `event`.timestamp",
      "AND IS_FIRST(`event`.timestamp)"
    ))
  )
})
