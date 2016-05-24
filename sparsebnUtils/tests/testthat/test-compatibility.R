context("compatibility")

sbf.empty <- generate_empty_sparsebnFit()
sbp.empty <- generate_empty_sparsebnPath()
sbp <- generate_fixed_sparsebnPath()
sbf <- generate_fixed_sparsebnFit()

test_that("Coercion works on empty graphs", {
    ### sparsebnFit
    expect_equivalent(sbf.empty, to_edgeList(to_graphNEL(sbf.empty)))
    expect_equivalent(sbf.empty, to_edgeList(to_igraph(sbf.empty)))
    expect_equivalent(sbf.empty, to_edgeList(to_network(sbf.empty)))

    ### sparsebnPath
    expect_equivalent(sbp.empty, to_edgeList(to_graphNEL(sbp.empty)))
    expect_equivalent(sbp.empty, to_edgeList(to_igraph(sbp.empty)))
    expect_equivalent(sbp.empty, to_edgeList(to_network(sbp.empty)))
})

test_that("Coercion to graphNEL works", {
    ### Check that output is a valid graphNEL object
    expect_is(to_graphNEL(sbf)$edges, "graphNEL")
    expect_true(check_list_class(lapply(to_graphNEL(sbp), function(x) x$edges), "graphNEL"))

    ### Check that nothing changes
    expect_equivalent(sbf, to_edgeList(to_graphNEL(sbf)))
    expect_equivalent(sbp, to_edgeList(to_graphNEL(sbp)))
})

test_that("Coercion to igraph works", {
    ### Check that output is a valid igraph object
    expect_is(to_igraph(sbf)$edges, "igraph")
    expect_true(check_list_class(lapply(to_igraph(sbp), function(x) x$edges), "igraph"))

    ### Check that nothing changes
    expect_equivalent(sbf, to_edgeList(to_igraph(sbf)))
    expect_equivalent(sbp, to_edgeList(to_igraph(sbp)))
})

test_that("Coercion to network works", {
    ### Check that output is a valid network object
    expect_is(to_network(sbf)$edges, "network")
    expect_true(check_list_class(lapply(to_network(sbp), function(x) x$edges), "network"))

    ### Check that nothing changes
    expect_equivalent(sbf, to_edgeList(to_network(sbf)))
    expect_equivalent(sbp, to_edgeList(to_network(sbp)))
})

test_that("Iterative coercing works", {

    ### Empty graphs: sparsebnFit
    expect_equivalent(sbf.empty, to_edgeList(to_graphNEL(to_igraph(to_network(sbf.empty)))))
    expect_equivalent(sbf.empty, to_edgeList(to_graphNEL(to_network(to_igraph(sbf.empty)))))

    expect_equivalent(sbf.empty, to_edgeList(to_igraph(to_graphNEL(to_network(sbf.empty)))))
    expect_equivalent(sbf.empty, to_edgeList(to_igraph(to_network(to_graphNEL(sbf.empty)))))

    expect_equivalent(sbf.empty, to_edgeList(to_network(to_igraph(to_graphNEL(sbf.empty)))))
    expect_equivalent(sbf.empty, to_edgeList(to_network(to_graphNEL(to_igraph(sbf.empty)))))

    ### Non-empty graphs: sparsebnFit
    expect_equivalent(sbf, to_edgeList(to_graphNEL(to_igraph(to_network(sbf)))))
    expect_equivalent(sbf, to_edgeList(to_graphNEL(to_network(to_igraph(sbf)))))

    expect_equivalent(sbf, to_edgeList(to_igraph(to_graphNEL(to_network(sbf)))))
    expect_equivalent(sbf, to_edgeList(to_igraph(to_network(to_graphNEL(sbf)))))

    expect_equivalent(sbf, to_edgeList(to_network(to_igraph(to_graphNEL(sbf)))))
    expect_equivalent(sbf, to_edgeList(to_network(to_graphNEL(to_igraph(sbf)))))

    ### Non-empty graphs: sparsebnPath
    expect_equivalent(sbp, to_edgeList(to_graphNEL(to_igraph(to_network(sbp)))))
    expect_equivalent(sbp, to_edgeList(to_graphNEL(to_network(to_igraph(sbp)))))

    expect_equivalent(sbp, to_edgeList(to_igraph(to_graphNEL(to_network(sbp)))))
    expect_equivalent(sbp, to_edgeList(to_igraph(to_network(to_graphNEL(sbp)))))

    expect_equivalent(sbp, to_edgeList(to_network(to_igraph(to_graphNEL(sbp)))))
    expect_equivalent(sbp, to_edgeList(to_network(to_graphNEL(to_igraph(sbp)))))

})
