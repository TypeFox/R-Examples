context("sparse.list")

test_that("check list structure", {
    ### Throws error when input is not a list
    expect_error(sparse.list(integer(0)), "Input must be a list!")

    ### Throws error when input is a list without a names attribute
    expect_error(sparse.list(list()), "coercable to an object of type sparse")

    ### Throws error when list has wrong names
    li <- list(not_rows = integer(0), not_cols = integer(0), not_vals = numeric(0), not_dim = integer(0), not_start = integer(0))
    expect_error(sparse.list(li), "coercable to an object of type sparse")

    ### Throws error when list has correct names but no values
    li <- list(rows = integer(0), cols = integer(0), vals = numeric(0), dim = integer(0), start = integer(0))
    expect_error(sparse.list(li), "dim attribute")
})

test_that("checks that dim is properly set", {
    ### Throws errors when dim is improperly set
    li <- list(rows = integer(0), cols = integer(0), vals = numeric(0), dim = integer(0), start = integer(0))
    li[["start"]] <- 0
    expect_error(sparse.list(li), "dim attribute")
    li[["start"]] <- 1
    expect_error(sparse.list(li), "dim attribute")
})

test_that("checks that start is properly set", {
    ### Throws errors when start is improperly set
    li <- list(rows = integer(0), cols = integer(0), vals = numeric(0), dim = c(1,1), start = integer(0))
    expect_error(sparse.list(li), "missing value where TRUE/FALSE needed")
    li[["start"]] <- -1
    expect_error(sparse.list(li), "start attribute")
    li[["start"]] <- 2
    expect_error(sparse.list(li), "start attribute")

})

test_that("checks that consistency of rows / cols / vals is properly set", {
    ### Accepts empty input for rows, cols, vals as long as dim & start are properly set
    li <- list(rows = integer(0), cols = integer(0), vals = numeric(0), dim = c(1,1), start = 1)
    expect_error(sparse.list(li), NA)

    ### Throws error when any of rows / cols / vals have unequal length
    li <- list(rows = c(0L), cols = integer(0), vals = numeric(0), dim = c(1,1), start = 0)
    expect_error(sparse.list(li), "elements have different sizes")
    li <- list(rows = integer(0), cols = c(0L), vals = numeric(0), dim = c(1,1), start = 0)
    expect_error(sparse.list(li), "elements have different sizes")
    li <- list(rows = integer(0), cols = integer(0), vals = c(0), dim = c(1,1), start = 0)
    expect_error(sparse.list(li), "elements have different sizes")
    li <- list(rows = c(0L), cols = integer(0), vals = c(0), dim = c(1,1), start = 0)
    expect_error(sparse.list(li), "elements have different sizes")
    li <- list(rows = integer(0), cols = c(0L), vals = c(0), dim = c(1,1), start = 0)
    expect_error(sparse.list(li), "elements have different sizes")
    li <- list(rows = c(0L), cols = c(0L), vals = numeric(0), dim = c(1,1), start = 0)
    expect_error(sparse.list(li), "elements have different sizes")

    ### Accepts input with equal lengths
    li <- list(rows = c(0L), cols = c(0L), vals = c(0), dim = c(1,1), start = 0)
    expect_error(sparse.list(li), NA)
})

test_that("checks types for rows / cols / vals", {
    ### Throws error if rows or cols are not integers
    li <- list(rows = c(pi), cols = c(0L), vals = c(0), dim = c(1,1), start = 0)
    expect_error(sparse.list(li))
    li <- list(rows = c(0L), cols = c(pi), vals = c(0), dim = c(1,1), start = 0)
    expect_error(sparse.list(li))

    ### Throws error if vals is not numeric
    li <- list(rows = c(0L), cols = c(0L), vals = c("a"), dim = c(1,1), start = 0)
    expect_error(sparse.list(li))
})

# test_that("constructor", {
#     expect_is(sparse(list()))
# })
