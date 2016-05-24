context("Debugging append")

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(NULL, "part0"), {
            with(test.dataset(newDatasetFromFixture("apidocs"), "part1"), {
                with(test.dataset(newDatasetFromFixture("apidocs"), "part2"), {
                    exclusion(part1) <- part1$q1 == "Dog"
                    exclusion(part2) <- part2$q1 == "Dog"
                    part0 <- suppressMessages(appendDataset(part0, part1))
                    part0 <- suppressMessages(appendDataset(part0, part2))
                    test_that("Appending happened, and rows were excluded", {
                        expect_identical(dim(part0),
                            c(nrow(part1)*2L, ncol(part1)))
                        expect_equivalent(table(part0$q1)["Dog"], 0)
                    })
                })
            })
        })
    })
}
