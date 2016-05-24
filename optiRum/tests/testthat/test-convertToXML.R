
context("convertToXML")

test_that("data.frame works", {
    df <- data.frame(a = 1:10, b = letters[1:10])
    expect_true(inherits(convertToXML(df), "XMLInternalDOM"))
})

test_that("data.table works", {
    # library(data.table)
    dt <- data.table(a = 1:10, b = letters[1:10])
    expect_true(inherits(convertToXML(dt), "XMLInternalDOM"))
})

test_that("table works", {
    tab <- table(a = 1:10, b = letters[1:10])
    expect_true(inherits(convertToXML(tab), "XMLInternalDOM"))
})

test_that("0 rows data.frame works", {
    df <- data.frame(a = 1:10, b = letters[1:10])
    expect_true(inherits(convertToXML(df[0, ]), "XMLInternalDOM"))
})

test_that("0 rows data.table works", {
    # library(data.table)
    dt <- data.table(a = 1:10, b = letters[1:10])
    expect_true(inherits(convertToXML(dt[0, ]), "XMLInternalDOM"))
}) 
