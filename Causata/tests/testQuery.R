# Test querying causata Customers table
# 
# Author: David Barker
###############################################################################
library(testthat)
library(Causata)

context("Query")

equals <- testthat::equals

source("utils.R")
has.causata.connection <- has.causata.connection.function()

test_that("Blank Query is SELECT * FROM Customers variable", {
  expect_that(as.character(Query()), equals("SELECT * FROM Customers variable"))
})

test_that("Limit can be added to the blank query", {
  q <- Query()
  limited <- q + Limit(100)
  expect_that(as.character(q), equals("SELECT * FROM Customers variable"))
  expect_that(as.character(limited), equals("SELECT * FROM Customers variable LIMIT 100"))
})

test_that("Variables can be specified", {
  q <- Query() + WithVariables("purchase", "page-view")
  expect_that(as.character(q), equals("SELECT `purchase`,`page-view` FROM Customers variable"))
})

test_that("Where clauses can be added", {
  q <- Query() + Where("purchase", EqualTo("letters"))
  expect_that(as.character(q), equals("SELECT * FROM Customers variable WHERE variable.`purchase` = 'letters'"))
})

test_that("Two Where clauses and a limit", {
  q <- Query() +
        Where("purchase", EqualTo("letters")) +
        Where("page-view-count", GreaterThan(10)) +
        Limit(1000)
  expect_that(as.character(q), equals("SELECT * FROM Customers variable WHERE variable.`purchase` = 'letters' AND variable.`page-view-count` > 10 LIMIT 1000"))
})

test_that("Variables, filters and limit specified", {
  q <- Query() +
    WithVariables("variable-one", "variable-two") +
    Where("purchase", EqualTo("letters")) +
    Where("page-view-count", GreaterThan(10)) +
    Limit(1000)
  expect_that(as.character(q), equals("SELECT `variable-one`,`variable-two` FROM Customers variable WHERE variable.`purchase` = 'letters' AND variable.`page-view-count` > 10 LIMIT 1000"))
})

test_that("Variables, filters and limit specified", {
  q <- Query() +
    WithVariables("variable-one", "variable-two") +
    Where("purchase", EqualTo("letters")) +
    Where("page-view-count", GreaterThan(10)) +
    Limit(1000)
  expect_that(as.character(q), equals("SELECT `variable-one`,`variable-two` FROM Customers variable WHERE variable.`purchase` = 'letters' AND variable.`page-view-count` > 10 LIMIT 1000"))
})

if (has.causata.connection){
  test_that("Selecting one customer by customer attribute", 
    with.local.connection(function(conn) {
      causata.variables <- list("customer-id"="SELECT PROPERTY customer-id")
      with.primary.variables(conn, causata.variables, {
        variables <- GetMetadata(conn)
        expect_that('customer-id' %in% variables$system.name, is_true())
      })
    })
  )
}

