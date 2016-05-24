context("Hiding variables")

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            var1 <- ds[[1]]
            test_that("Hide and unhide method for variables", {
                expect_true(name(var1) %in% findVariables(ds, key="name",
                    value=TRUE))
                var1 <- hide(var1)
                ds <- refresh(ds)
                expect_false(name(var1) %in% findVariables(ds, key="name",
                    value=TRUE))

                var1 <- unhide(var1)
                ds <- refresh(ds)
                expect_true(name(var1) %in% findVariables(ds, key="name",
                    value=TRUE))
            })
        })

        with(test.dataset(df), {
            test_that("There are no hidden variables to start", {
                expect_equivalent(index(hidden(ds)), list())
                expect_identical(hiddenVariables(ds), c())
                expect_identical(dim(ds), dim(df))
            })

            try(ds <- hideVariables(ds, c("v2", "v3")))
            test_that("hideVariables hides by alias", {
                expect_identical(names(ds)[1:2], c("v1", "v4"))
                expect_identical(hiddenVariables(ds), c("v2", "v3"))
                expect_identical(length(hidden(ds)), 2L)
                expect_identical(length(variables(ds)), ncol(df) - 2L)
                expect_identical(dim(ds), c(nrow(df), ncol(df) - 2L))
            })

            try(hiddenVariables(ds) <- "v3")
            ## work like is.na<-, i.e. adds but doesn't unhide by omitting
            test_that("hiddenVariables<- does nothing if already hidden", {
                expect_identical(hiddenVariables(ds), c("v2", "v3"))
                expect_identical(names(ds)[1:2], c("v1", "v4"))
                expect_identical(dim(ds), c(nrow(df), ncol(df) - 2L))
            })

            try(hiddenVariables(ds) <- "v4")
            test_that("hiddenVariables<- adds variables", {
                expect_identical(names(ds)[1:2], c("v1", "v5"))
                expect_identical(hiddenVariables(ds), c("v2", "v3", "v4"))
                expect_identical(dim(ds), c(nrow(df), ncol(df)-3L))
            })

            test_that("hidden variables can be accessed with $", {
                expect_warning(ds$v2, "hidden")
                expect_true(is.Text(suppressWarnings(ds$v2)))
            })

            try(ds <- unhideVariables(ds, c("v2", "v3", "v4")))

            test_that("unhideVariables by alias", {
                expect_identical(hiddenVariables(ds), c())
                expect_identical(dim(ds), dim(df))
                expect_warning(ds$v2, NA)
                expect_true(is.Text(ds$v2))
            })
        })

        with(test.dataset(df), {
            test_that("hideVariables with grep (and by index)", {
                ds <- hideVariables(ds, pattern="v[23]")
                expect_identical(names(ds)[1:2], c("v1", "v4"))

                ds <- unhideVariables(ds, pattern="v[23]")
                expect_identical(hiddenVariables(ds), c())
            })

            test_that("Error handling", {
                expect_identical(hiddenVariables(ds), c()) # To be clear
                ## Need something better than subscript out of bounds, probably

            })
        })

        with(test.dataset(df), {
            test_that("can hide variables by group", {
                ordering(ds) <- VariableOrder(
                    VariableGroup(name="g1", variables=list(ds$v1)),
                    VariableGroup(name="group2", variables=ds[c("v3", "v4")])
                )
                expect_identical(length(grouped(ordering(ds))), 2L)
                ds <- hideVariables(ds, ungrouped(ordering(ds)))
                expect_identical(length(hiddenVariables(ds)), ncol(df) - 3L)
                expect_true(all(c("v2", "v5") %in% hiddenVariables(ds)))
            })
        })

        with(test.dataset(mrdf), {
            ds <- mrdf.setup(ds)
            test_that("Can hide array variables", {
                expect_true("CA" %in% names(ds))
                try(hiddenVariables(ds) <- c("CA", "v4"))
                expect_false("CA" %in% names(ds))
            })
        })
        with(test.dataset(mrdf), {
            ds <- mrdf.setup(ds, selections="1.0")
            test_that("Can hide MR variables", {
                expect_true("MR" %in% names(ds))
                try(hiddenVariables(ds) <- "MR")
                expect_false("MR" %in% names(ds))
            })
        })

        with(test.dataset(mrdf), {
            ds <- mrdf.setup(ds, pattern="mr_1")
            test_that("Can hide array variables even if they only have one subvar", {
                expect_identical(names(ds), c("CA", "mr_2", "mr_3", "v4"))
                expect_identical(length(subvariables(ds$CA)), 1L)
                try(hiddenVariables(ds) <- "CA")
                expect_false("CA" %in% names(ds))
            })
        })
    })
}
