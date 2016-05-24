
###
###  phylobase.options
###

context("phylobase.options()")

test_that("test of match.arg", {
    op <- phylobase.options()
    ## test match.arg
    expect_error(phylobase.options(retic="test"))
    no <- phylobase.options(retic="f")
    expect_equal(no$retic, "fail")
    phylobase.options(op)
})

test_that("test of multiple arguments", {
    op <- phylobase.options()
    ## test multiple args
    no <- phylobase.options(retic="f", poly="f")
    expect_equal(no$retic, "fail")
    expect_equal(no$poly, "fail")
    phylobase.options(op)
})

test_that("test some failures", {
    op <- phylobase.options()
    ## check some failures
    expect_error(phylobase.options(1))
    expect_error(phylobase.options("foobar"="foo"))
    phylobase.options(op)
})
