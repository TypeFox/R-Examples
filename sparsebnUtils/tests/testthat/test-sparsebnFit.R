context("sparsebnFit")

## True num.nodes = 5
## True num.edges = 5
edges <- generate_fixed_edgeList()

test_that("sparsebnFit correctly identifies when input number of nodes is not consistent", {
    li <- list(edges = edges, lambda = pi, nedge = num.edges(edges), pp = 1, nn = 10, time = exp(1))
    expect_error(cf <- sparsebnFit(li))
})

test_that("sparsebnFit correctly identifies when input number of edges is not consistent", {
    li <- list(edges = edges, lambda = pi, nedge = 20, pp = num.nodes(edges), nn = 10, time = exp(1))
    expect_error(cf <- sparsebnFit(li))
})

test_that("sparsebnFit is consistent with different ways of accessing nedge", {
    li <- list(edges = edges, lambda = pi, nedge = 5, pp = 5, nn = 10, time = exp(1))
    cf <- sparsebnFit(li) ### Should not generate an error anymore

    matrix.nedge <- Matrix::nnzero(get.adjacency.matrix(cf$edges))
    edgeL.nedge <- num.edges(edges)

    ### Note that sbm.nedge and cf$nedge are equal by construction since nedge is set manually above
    expect_equal(edgeL.nedge, matrix.nedge, cf$nedge)
})
