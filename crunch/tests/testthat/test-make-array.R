context("Categorical Array")

if (run.integration.tests) {
    with(test.authentication, {
        test_that("can make Categorical Array with Dataset subset", {
            with(test.dataset(mrdf), {
                var <- makeArray(ds[1:3], name="arrayVar")
                expect_true(inherits(var, "VariableDefinition"))
                ds$arrayVar <- var
                expect_equal(c("arrayVar", "v4"), names(ds))
                expect_true(is.CA(ds$arrayVar))
                name(ds$arrayVar) <- "TESTONE"
                ds <- refresh(ds)
                expect_equal(c("TESTONE", "v4"), names(variables(ds)))
            })
        })
        test_that("Can make and then delete an array (and its subvars delete too)", {
            with(test.dataset(mrdf), {
                ds$arrayVar <- makeArray(ds[c("mr_1", "mr_2", "mr_3")],
                    name="arrayVar")
                expect_equal(c("arrayVar", "v4"), names(ds))
                with(consent(), {
                    ## delete array variable
                    ds$arrayVar <- NULL
                })
                expect_identical(names(ds), "v4")
                expect_identical(ncol(ds), 1L)
            })
        })
        test_that("can make Categorical Array with pattern", {
            with(test.dataset(mrdf), {
                var <- makeArray(pattern="mr_[123]", dataset=ds,
                    name="arrayVar")
                expect_true(inherits(var, "VariableDefinition"))
                ds$arrayVar <- var
                expect_equal(c("arrayVar", "v4"), names(ds))
                ## unbind.
                u <- unbind(ds$arrayVar)
                ds <- refresh(ds)
                expect_true(setequal(names(ds), names(mrdf)))
                expect_identical(ncol(ds), 4L)
            })
        })

        with(test.dataset(mrdf), {
            test_that("makeArray error conditions", {
                no.name <- "Must provide the name for the new variable"
                no.match <- "Pattern did not match any variables"
                ds.mismatch <- "`list_of_variables` must be from `dataset`"
                expect_error(makeArray(), no.name)
                expect_error(makeArray(pattern="mr_[123]", dataset=ds),
                    no.name)
                expect_error(makeArray(pattern="rm_", dataset=ds,
                    name="foo"), no.match)
                expect_silent(ds$d <- is.CA(makeArray(c("mr_1", "mr_2", "mr_3"),
                    dataset=ds, name="foo")))
                skip("Errors, but with wrong error condition")
                with(test.dataset(df, "notds"), {
                    expect_error(makeArray(ds[1:3], dataset=notds,
                        name="arrayVar"), ds.mismatch)
                })
            })
        })

        with(test.dataset(mrdf), {
            ds$arrayVar <- makeArray(pattern="mr_[123]", dataset=ds,
                name="arrayVar")
            var <- ds$arrayVar
            test_that("setup to make MultipleResponse from CategoricalArray", {
                expect_true(is.CA(var))
            })

            test_that("can make MultipleResponse from CategoricalArray by editing category$selected", {
                categories(var)[[1]]$selected <- TRUE
                var <- refresh(var) ## Refresh required if changing type by editing categories
                expect_true(is.Multiple(var))
                categories(var)[[1]]$selected <- FALSE
                var <- refresh(var) ## Refresh required if changing type by editing categories
                expect_true(is.CA(var))
            })

            test_that("can make MultipleResponse from CategoricalArray by dichotomizing categories (and back by undichotomize)", {
                categories(var) <- dichotomize(categories(var), 1)
                var <- refresh(var) ## Refresh required if changing type by editing categories
                expect_true(is.Multiple(var))
                categories(var) <- undichotomize(categories(var))
                var <- refresh(var) ## Refresh required if changing type by editing categories
                expect_true(is.CA(var))
            })
            test_that("can (un)dichotomize directly on the variable", {
                var <- dichotomize(var, 1)
                expect_true(is.Multiple(var))
                expect_true(is.Multiple(refresh(var)))
                var <- undichotomize(var)
                expect_true(is.CA(var))
                expect_true(is.CA(refresh(var)))
            })
            test_that("can (un)dichotomize on var in dataset", {
                ds <- refresh(ds)
                ds$arrayVar <- dichotomize(ds$arrayVar, 1)
                expect_true(is.Multiple(ds$arrayVar))
                expect_true(is.Multiple(refresh(ds)$arrayVar))
                ds$arrayVar <- undichotomize(ds$arrayVar)
                expect_true(is.CA(ds$arrayVar))
                expect_true(is.CA(refresh(ds)$arrayVar))
            })
        })

        with(test.dataset(mrdf), {
            test_that("can make MultipleResponse directly", {
                cast.these <- grep("mr_", names(ds))
                ds[cast.these] <- lapply(ds[cast.these],
                    castVariable, "categorical")
                ds$arrayVar <- makeMR(pattern="mr_[123]", dataset=ds,
                    name="arrayVar", selections="1.0")
                var <- ds$arrayVar
                expect_true(is.Multiple(var))

                var <- undichotomize(var)
                expect_true(is.CA(var))

                ## unbind.
                u <- unbind(var)
                ds <- refresh(ds)
                expect_true(setequal(names(ds), names(mrdf)))
                expect_identical(ncol(ds), 4L)
            })

            with(test.dataset(mrdf), {
                test_that("makeMR error conditions", {
                    no.name <- "Must provide the name for the new variable"
                    no.match <- "Pattern did not match any variables"
                    need.variables <- "Invalid list of Variables to combine"
                    ds.mismatch <- "`list_of_variables` must be from `dataset`"
                    no.selections <- paste("Must provide the names of the",
                        "category or categories that indicate the dichotomous",
                        "selection")
                    invalid.selection <- paste("not found in",
                        "variable's categories")
                    not.categorical <- "are not Categorical"
                    expect_error(makeMR(), no.name)
                    expect_error(makeMR(pattern="mr_[123]", dataset=ds),
                        no.name)
                    expect_error(makeMR(pattern="rm_", dataset=ds,
                        name="foo", selections="foo"), no.match)
                    expect_error(makeMR(c("mr_1", "mr_2", "mr_3"),
                        dataset=ds, name="foo", selections="foo"),
                        not.categorical)
                    expect_error(makeMR(pattern="mr_[123]", dataset=ds,
                        name="arrayVar", selections="Not a Selection!"),
                        not.categorical)
                    cast.these <- grep("mr_", names(ds))
                    ds[cast.these] <- lapply(ds[cast.these],
                        castVariable, "categorical")
                    expect_error(makeMR(pattern="mr_[123]", dataset=ds,
                        name="arrayVar", selections="Not a Selection!"),
                        invalid.selection)
                    expect_error(makeMR(pattern="mr_[123]", dataset=ds,
                        name="arrayVar"), no.selections)
                    skip("Errors, but with Error in match(x, table, nomatch = 0L) : \n  'match' requires vector arguments\n")
                    with(test.dataset(df, "notds"), {
                        expect_error(makeMR(ds[1:3], dataset=notds,
                            name="arrayVar", selections="1.0"), ds.mismatch)
                    })
                })
            })
        })
    })
}
