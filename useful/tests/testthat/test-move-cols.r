context('Checking that moving columns to the fron or back of data.frame works')

library(dplyr)

# create testing data.frame
theDF <- data.frame(A=1:10, B=11:20, C=1:10, D=11:20)

test_that('cols to front/back returns characters', {
    expect_is(colsToFront(theDF, c('B', 'C')), 'character')
    expect_is(colsToFront(theDF, c('C', 'B')), 'character')
    expect_is(colsToFront(theDF, c('C', 'C')), 'character')
    expect_is(colsToBack(theDF, c('B', 'C')), 'character')
    expect_is(colsToBack(theDF, c('C', 'B')), 'character')
    expect_is(colsToBack(theDF, c('C', 'C')), 'character')
})

test_that('cols to front/back returns the right number', {
    expect_equal(length(colsToFront(theDF, c('B', 'C'))), ncol(theDF))
    expect_equal(length(colsToFront(theDF, c('C', 'B'))), ncol(theDF))
    expect_equal(length(colsToFront(theDF, c('C', 'C'))), ncol(theDF) + 1)
    expect_equal(length(colsToBack(theDF, c('B', 'C'))), ncol(theDF))
    expect_equal(length(colsToBack(theDF, c('C', 'B'))), ncol(theDF))
    expect_equal(length(colsToBack(theDF, c('C', 'C'))), ncol(theDF) + 1)
})

test_that('cols to front/back returns the right names', {
    expect_equal(colsToFront(theDF, c('B', 'C')), c('B', 'C', 'A', 'D'))
    expect_equal(colsToFront(theDF, c('C', 'B')), c('C', 'B', 'A', 'D'))
    expect_equal(colsToFront(theDF, c('C', 'C')), c('C', 'C', 'A', 'B', 'D'))
    expect_equal(colsToBack(theDF, c('B', 'C')), c('A', 'D', 'B', 'C'))
    expect_equal(colsToBack(theDF, c('C', 'B')), c('A', 'D', 'C', 'B'))
    expect_equal(colsToBack(theDF, c('C', 'C')), c('A', 'B', 'D', 'C', 'C'))
})


test_that('cols to front/back returns data.frame', {
    expect_is(moveToFront(theDF, c('B', 'C')), 'data.frame')
    expect_is(moveToFront(theDF, c('C', 'B')), 'data.frame')
    expect_is(moveToFront(theDF, c('C', 'C')), 'data.frame')
    expect_is(moveToBack(theDF, c('B', 'C')), 'data.frame')
    expect_is(moveToBack(theDF, c('C', 'B')), 'data.frame')
    expect_is(moveToBack(theDF, c('C', 'C')), 'data.frame')
})

test_that('cols to front/back returns the right dimension', {
    expect_equal(dim(moveToFront(theDF, c('B', 'C'))), dim(theDF))
    expect_equal(dim(moveToFront(theDF, c('C', 'B'))), dim(theDF))
    expect_equal(dim(moveToFront(theDF, c('C', 'C'))), dim(theDF))
    expect_equal(dim(moveToBack(theDF, c('B', 'C'))), dim(theDF))
    expect_equal(dim(moveToBack(theDF, c('C', 'B'))), dim(theDF))
    expect_equal(dim(moveToBack(theDF, c('C', 'C'))), dim(theDF))
})

test_that('cols to front/back returns the right names', {
    expect_equal(names(moveToFront(theDF, c('B', 'C'))), c('B', 'C', 'A', 'D'))
    expect_equal(names(moveToFront(theDF, c('C', 'B'))), c('C', 'B', 'A', 'D'))
    expect_equal(names(moveToFront(theDF, c('C', 'C'))), c('C', 'A', 'B', 'D'))
    expect_equal(names(moveToBack(theDF, c('B', 'C'))), c('A', 'D', 'B', 'C'))
    expect_equal(names(moveToBack(theDF, c('C', 'B'))), c('A', 'D', 'C', 'B'))
    expect_equal(names(moveToBack(theDF, c('C', 'C'))), c('A', 'B', 'D', 'C'))
})