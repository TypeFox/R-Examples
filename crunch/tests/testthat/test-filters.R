context("Filters")

test_that("show method exists", {
    expect_true(is.character(capture.output(print(CrunchFilter()))))
})

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            test_that("We have an empty filter catalog", {
                expect_true(inherits(filters(ds), "FilterCatalog"))
                expect_identical(length(filters(ds)), 0L)
                expect_output(filters(ds), capture.output(print(data.frame())))
            })

            filters(ds)[["Test filter"]] <- ds$v4 == "B"
            test_that("We can create a filter", {
                expect_identical(length(filters(ds)), 1L)
                expect_identical(names(filters(ds)), "Test filter")
                expect_identical(name(filters(ds)[[1]]), "Test filter")
            })
            test_that("Show methods for filter and filter catalog", {
                expect_output(filters(ds),
                    paste(capture.output(print(data.frame(
                            name="Test filter",
                            id=filters(ds)[["Test filter"]]@body$id,
                            is_public=FALSE
                        ))), collapse="\n"), fixed=TRUE)
                expect_output(filters(ds)[["Test filter"]],
                    paste0('Crunch filter ', dQuote("Test filter"),
                        '\nExpression: v4 == "B"'))
            })

            test_that("We can make it public/private", {
                expect_false(is.public(filters(ds)[["Test filter"]]))
                is.public(filters(ds)[["Test filter"]]) <- TRUE
                expect_true(is.public(filters(ds)[["Test filter"]]))
                is.public(filters(ds)[["Test filter"]]) <- FALSE
                expect_false(is.public(filters(ds)[["Test filter"]]))
            })

            test_that("Setter/getter by index", {
                expect_false(is.public(filters(ds)[[1]]))
                is.public(filters(ds)[[1]]) <- TRUE
                expect_true(is.public(filters(ds)[[1]]))
                is.public(filters(ds)[[1]]) <- FALSE
                expect_false(is.public(filters(ds)[[1]]))
            })

            test_that("Can update a filter's expression by name", {
                expect_json_equivalent(zcl(expr(filters(ds)[[1]])),
                    zcl(ds$v4 == "B"))
                filters(ds)[["Test filter"]] <- ds$v4 == "C"
                expect_json_equivalent(zcl(expr(filters(ds)[[1]])),
                    zcl(ds$v4 == "C"))
            })
            test_that("Can update a filter's expression by index", {
                expect_json_equivalent(zcl(expr(filters(ds)[[1]])),
                    zcl(ds$v4 == "C"))
                filters(ds)[[1]] <- ds$v4 == "B"
                expect_json_equivalent(zcl(expr(filters(ds)[[1]])),
                    zcl(ds$v4 == "B"))
            })
            test_that("Error handling for [[<-", {
                expect_error(filters(ds)[[6]] <- ds$v4 == "B",
                    "Subscript out of bounds: 6")
            })

            test_that("We have an applied filters view", {
                expect_identical(length(appliedFilters(ds)), 0L)
            })

            test_that("We can 'apply' a filter", {
                appliedFilters(ds) <- filters(ds)[["Test filter"]]
                expect_identical(length(appliedFilters(ds)), 1L)
            })

            test_that("'applied filters' for the UI don't affect R", {
                expect_identical(length(appliedFilters(ds)), 1L)
                validImport(ds)
            })

            test_that("We also have 'active filter' for the R object", {
                expect_identical(activeFilter(ds), NULL)
            })

            test_that("We can set 'active filter'", {
                activeFilter(ds) <- ds$v4 == "C"
                expect_identical(zcl(activeFilter(ds)), zcl(ds$v4 == "C"))
            })

            test_that("If we set an active filter, cubes will be filtered by it (and not UI filters)", {
                activeFilter(ds) <- ds$v4 == "C"
                expect_equivalent(as.array(crtabs(~ v4, data=ds)),
                    array(c(0, 10), dim=2L, dimnames=list(v4=c("B", "C"))))
            })
        })
    })
}
