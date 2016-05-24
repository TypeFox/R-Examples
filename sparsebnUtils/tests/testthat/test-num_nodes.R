context("num.nodes")

test_that("num.nodes works on edgeLists", {
    ### Trivial case
    expect_equal(num.nodes(generate_empty_edgeList()), 1)

    ### Non-trivial case
    edgeL <- generate_fixed_edgeList()
    expect_equal(num.nodes(edgeL), 5)
})

test_that("num.nodes works on sparsebnFit", {
    ### Trivial case
    cf <- generate_empty_sparsebnFit()
    expect_equal(num.nodes(cf), 1)

    ### Non-trivial case
    cf <- generate_fixed_sparsebnFit()
    expect_equal(num.nodes(cf), 5)
})

test_that("num.nodes works on sparsebnPath", {
    ### Trivial case
    cp <- generate_empty_sparsebnPath()
    expect_equal(num.nodes(cp), 1)

    ### Non-trivial case
    cp <- generate_fixed_sparsebnPath()
    expect_equal(num.nodes.sparsebnPath(cp), 5)
})
