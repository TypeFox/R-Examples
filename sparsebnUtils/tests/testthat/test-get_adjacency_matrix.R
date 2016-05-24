context("get.adjacency.matrix")

test_that("get.adjacency.matrix works on empty graphs", {
    adj <- generate_empty_adjacency_matrix()

    edgeL <- generate_empty_edgeList()
    expect_equivalent(as.matrix(get.adjacency.matrix(edgeL)), adj) # Ignore attributes using expect_equivalent()

#     sbm <- generate_empty_SparseBlockMatrixR()
#     expect_equivalent(as.matrix(get.adjacency.matrix(sbm)), adj)

    cf <- generate_empty_sparsebnFit()
    expect_equivalent(as.matrix(get.adjacency.matrix(cf)), adj)

    cp <- generate_empty_sparsebnPath()
    expect_equivalent(lapply(get.adjacency.matrix(cp), as.matrix), list(adj, adj, adj, adj))
})

test_that("get.adjacency.matrix works on nontrivial graphs", {
    adj <- generate_fixed_adjacency_matrix()

    edgeL <- generate_fixed_edgeList()
    expect_equivalent(as.matrix(get.adjacency.matrix(edgeL)), adj) # Ignore attributes using expect_equivalent()

#     sbm <- generate_fixed_SparseBlockMatrixR()
#     expect_equivalent(as.matrix(get.adjacency.matrix(sbm)), adj)

    cf <- generate_fixed_sparsebnFit()
    expect_equivalent(as.matrix(get.adjacency.matrix(cf)), adj)

    cp <- generate_fixed_sparsebnPath()
    expect_equivalent(lapply(get.adjacency.matrix(cp), as.matrix), list(adj, adj, adj, adj))
})

