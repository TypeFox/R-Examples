context("Learn")

test_that("Exceptions", { 
    expect_that(learn(), throws_error("Set is missing"))
})

test_that("Training set", { 
    data("norms", package="pnn")
    m <- learn(norms)
    expect_that(m$set[,-1], equals(norms[,-1]))
    expect_that(length(m$set[,1]), equals(400))
    expect_that(m$n, equals(400)) 
    n <- learn(norms, nn=m)
    expect_that(length(n$set[,1]), equals(800)) 
    expect_that(n$n, equals(800)) 
})

test_that("PNN parameters", { 
    data("norms", package="pnn")
    m <- learn(norms)
    expect_that(m$category.column, equals(1))
    expect_that(m$categories, equals(c("A","B")))
    expect_that(m$k, equals(2))
})
