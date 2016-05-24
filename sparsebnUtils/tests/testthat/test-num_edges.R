context("num.edges")

test_that("num.edges works on edgeLists", {
    ### Trivial case
    expect_equal(num.edges(generate_empty_edgeList()), 0)

    ### Non-trivial case
    edgeL <- generate_fixed_edgeList()
    expect_equal(num.edges(edgeL), 5)
})

test_that("num.edges works on sparsebnFit", {
    ### Trivial case
    cf <- generate_empty_sparsebnFit()
    expect_equal(num.edges(cf), 0)

    ### Non-trivial case
    cf <- generate_fixed_sparsebnFit()
    expect_equal(num.edges(cf), 5)
})

test_that("num.edges works on sparsebnPath", {
    ### Trivial case
    cp <- generate_empty_sparsebnPath()
    expect_equal(num.edges(cp), rep(0, length(cp)))

    ### Non-trivial case
    cp <- generate_fixed_sparsebnPath()
    expect_equal(num.edges(cp), rep(5, length(cp)))
})
