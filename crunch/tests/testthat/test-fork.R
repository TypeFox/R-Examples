context("Fork and merge")

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            test_that("Fork catalog exists", {
                expect_true(inherits(forks(ds), "ForkCatalog"))
                expect_identical(length(forks(ds)), 0L)
            })

            f1 <- forkDataset(ds)
            f1.name <- paste("Fork of", name(ds))
            test_that("Can create a new fork", {
                expect_true(is.dataset(f1))
                expect_identical(name(f1), f1.name)
                expect_identical(names(forks(ds)), f1.name)
            })

            f2 <- forkDataset(ds, "Fork yeah!")
            test_that("Can create a fork with a given name", {
                expect_true(is.dataset(f2))
                expect_identical(name(f2), "Fork yeah!")
                expect_true("Fork yeah!" %in% names(forks(ds)))
            })

            f3 <- forkDataset(ds)
            f3.name <- paste("Fork #3 of", name(ds))
            test_that("Creating forks autonames with a fork number", {
                expect_true(is.dataset(f3))
                expect_identical(name(f3), f3.name)
                expect_true(setequal(names(forks(ds)),
                    c(f1.name, "Fork yeah!", f3.name)))
            })

            delete(f2)
            delete(f3)
            test_that("If you delete a fork, it disappears from upstream forks catalog", {
                expect_identical(names(refresh(forks(ds))), f1.name)
            })

            ## Make edits to fork #1. cf. test-versioning.R:
            # 1. Edit variable metadata
            names(categories(f1$v4))[1:2] <- c("d", "e")
            name(f1$v2) <- "Variable Two"
            description(f1$v3) <- "The third variable in the dataset"

            # 2. Edit dataset metadata
            description(f1) <- "A dataset for testing"

            # 3. Reorder variables
            ordering(f1) <- VariableOrder(VariableGroup("Even", f1[c(2,4,6)]),
                VariableGroup("Odd", f1[c(1,3,5)]))

            # 5. Add non-derived variable
            f1$v8 <- rep(1:5, 4)

            # 4. Derive variable
            f1$v7 <- f1$v3 - 6

            ## Assert those things
            test_that("The edits are made to the fork", {
                expect_identical(names(na.omit(categories(f1$v4))),
                    c("d", "e"))
                expect_identical(name(f1$v2), "Variable Two")
                expect_identical(description(f1$v3),
                    "The third variable in the dataset")
                expect_identical(description(f1), "A dataset for testing")
                expect_identical(as.vector(f1$v7), df$v3 - 6)
                expect_equivalent(as.vector(f1$v8), rep(1:5, 4))
                expect_identical(aliases(variables(f1)),
                    paste0("v", c(2,4,6,1,3,5,8,7)))
            })

            test_that("The upstream dataset is unaffected by edits to the fork", {
                validImport(ds)
            })

            ## So that f1 gets cleaned up even if merge fails
            with(test.dataset(f1, "f1"), {
                ## Now merge f1 back to ds
                ds <- mergeFork(ds, f1)
                test_that("The edits made to the fork are now upstream", {
                    expect_identical(names(na.omit(categories(ds$v4))),
                        c("d", "e"))
                    expect_identical(name(ds$v2), "Variable Two")
                    expect_identical(description(ds$v3),
                        "The third variable in the dataset")
                    expect_identical(as.vector(ds$v7), df$v3 - 6)
                    expect_equivalent(as.vector(ds$v8), rep(1:5, 4))
                    expect_identical(aliases(variables(ds)),
                        paste0("v", c(2,4,6,1,3,5,8,7)))
                    ## Extra checks for v7 and v8
                    expect_true("v7" %in% aliases(allVariables(ds)))
                    expect_true("v8" %in% aliases(allVariables(ds)))
                })
                test_that("Certain changes don't merge", {
                    expect_identical(description(ds), "")
                })
            })
        })
    })
}
