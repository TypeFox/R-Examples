context('Adding subvariables')

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(newDatasetFromFixture("apidocs")), {
            test_that("Adding a subvariable to an array", {
                ds$petloc_daycare <- VariableDefinition(factor(rep(c("Cat",
                    "Dog"), 10)), name="doggy daycare")
                expect_identical(c("Home", "Work"),
                    names(subvariables(ds$petloc)))
                expect_true("petloc_daycare" %in% aliases(variables(ds)))
                addSubvariable(ds$petloc, ds$petloc_daycare)
                ds <- refresh(ds)
                expect_identical(c("Home", "Work", "doggy daycare"),
                    names(subvariables(ds$petloc)))
            })
        })
    })
}
