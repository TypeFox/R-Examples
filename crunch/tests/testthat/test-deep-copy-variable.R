context("Deep copies of variables")

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(newDatasetFromFixture("apidocs")), {
            test_that("Can deep copy categorical", {
                ds$q1a <- copy(ds$q1, deep=TRUE)
                expect_identical(as.vector(ds$q1), as.vector(ds$q1a))
            })
            test_that("Can deep copy numeric", {
                ds$ndogsa <- copy(ds$ndogs, deep=TRUE)
                expect_identical(as.vector(ds$ndogs), as.vector(ds$ndogsa))
            })
            test_that("Can deep copy datetime", {
                ds$wavea <- copy(ds$wave, deep=TRUE)
                expect_identical(as.vector(ds$wave), as.vector(ds$wavea))
            })
            # with(temp.option(crunch.debug=TRUE), {
                test_that("Can deep copy multiple response", {
                    skip("WIP")
                    ds$allpetsa <- copy(ds$allpets, deep=TRUE)
                    expect_identical(as.vector(ds$allpets), as.vector(ds$allpetsa))
                    expect_equivalent(as.array(crtabs(~ allpets, data=ds)),
                        as.array(crtabs(~ allpetsa, data=ds)))
                })
                test_that("Can deep copy categorical array", {
                    skip("WIP")
                    ds$petloca <- copy(ds$petloc, deep=TRUE)
                    expect_identical(as.vector(ds$petloc), as.vector(ds$petloca))
                    expect_equivalent(as.array(crtabs(~ petloc, data=ds)),
                        as.array(crtabs(~ petloca, data=ds)))
                })
            # })

            with(test.dataset(newDatasetFromFixture("apidocs"), "part2"), {
                test_that("Deep copies don't get data when appending", {
                    out <- appendDataset(ds, part2)
                    ## Counts should be double in the original than in the copy
                    expect_equivalent(as.array(crtabs(~ q1, data=out)),
                        2 * as.array(crtabs(~ q1a, data=out)))
                    skip("WIP")
                    expect_equivalent(as.array(crtabs(~ petloc, data=out)),
                        2 * as.array(crtabs(~ petloca, data=out)))
                })
            })
        })
    })
}
