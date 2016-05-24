context("VariableDefinition")

test_that("VariableDefinition constructs", {
    expect_identical(class(VariableDefinition()), "VariableDefinition")
    expect_identical(class(VarDef(name="foo")), "VariableDefinition")
    expect_equivalent(VarDef(name="Short", description="More verbose"),
        list(name="Short", description="More verbose"))
    expect_true(inherits(VarDef(), "VariableDefinition"))
})

test_that("VarDef takes 'data' and does toVariable on it", {
    ## See test-add-variable.R
    expect_equivalent(VarDef(data=as.factor(rep(LETTERS[2:3], 3)),
        name="Test Cats"),
        structure(list(values=rep(1:2, 3),
            type="categorical",
            categories=list(
                list(id=1L, name="B", numeric_value=1L, missing=FALSE),
                list(id=2L, name="C", numeric_value=2L, missing=FALSE),
                list(id=-1L, name="No Data", numeric_value=NULL, missing=TRUE)
            ),
            name="Test Cats"), class="VariableDefinition"))
})

test_that("VarDef takes 'values', whatever they are", {
    ## Do any validation here? or just let the API reject invalid requests?
    expect_equivalent(VarDef(name="var", values=1:5),
        structure(list(name="var", values=1:5), class="VariableDefinition"))
})

test_that("toVariable returns VarDef", {
    expect_identical(class(toVariable(letters[1:3])), "VariableDefinition")
})

test_that("VarDef(data=VarDef, ...)", {
    expect_identical(VarDef(VarDef(as.factor(rep(LETTERS[2:3], 3))), name="x"),
        VarDef(data=as.factor(rep(LETTERS[2:3], 3)), name="x"))
})

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            test_that("Wrapping VarDef has same result as just ds<-", {
                ds$newvar <- VarDef(df$v4)
                expect_true("newvar" %in% names(ds))
                expect_true(is.Categorical(ds$newvar))
                expect_identical(name(ds$newvar), "newvar")
                expect_identical(as.vector(ds$newvar), as.vector(ds$v4))
            })

            test_that("Extra attrs come in", {
                ds$newvar2 <- VarDef(df$v3, name="New var 2",
                    description="Second new variable")
                expect_true("newvar2" %in% names(ds))
                expect_true(is.Numeric(ds$newvar2))
                expect_identical(name(ds$newvar2), "New var 2")
                expect_identical(description(ds$newvar2), "Second new variable")
                expect_identical(as.vector(ds$newvar2), as.vector(ds$v3))
            })

            test_that("Can insert VarDef with no values", {
                expect_warning(ds$newvar3 <- VarDef(name="Empty", type="numeric"),
                    "Adding variable with no rows of data")
                expect_identical(as.vector(ds$newvar3), rep(NA_real_, 20L))
            })

            dropCache(self(ds)) ## Just so whether we have caching on doesn't affect the log we collect
            unifs <- runif(20)
            with(temp.option(httpcache.log=""), {
                avlog <- capture.output(ds <- addVariables(ds, list(
                    VarDef(1, name="One", description="the loneliest"),
                    VarDef(unifs, name="Some random stuff", alias="runif")
                )))
            })
            logdf <- loadLogfile(textConnection(avlog))
            test_that("token tests of log parsing while we're here", {
                expect_true(is.list(cacheLogSummary(logdf)))
                expect_true(is.list(requestLogSummary(logdf)))
                expect_true(is.numeric(blockTimings(logdf)))
                expect_error(getLogBlock(logdf, "NOBLOCK"),
                    "Block NOBLOCK not found")
            })
            test_that("addVariables can add multiple VarDefs", {
                expect_equivalent(as.vector(ds$One), rep(1, 20))
                expect_equivalent(as.vector(ds$runif), unifs)
            })
            test_that("addVariables doesn't refresh between each POST", {
                ## Parse avlog (and thus test the log parsing here)
                reqdf <- logdf[logdf$scope == "HTTP",]
                ## GET summary (nrows, to validate); POST var, POST var,
                ## with no GETs between the POSTs
                expect_identical(reqdf$verb[1:3], c("GET", "POST", "POST"))
                expect_identical(reqdf$url[2:3], rep(variableCatalogURL(ds), 2))
            })

            test_that("addVariables row validation happens before POSTing", {
                expect_error(addVariables(ds,
                    VarDef(2, name="Two"),
                    VarDef(1:3, name="Wrong num rows", alias="whatever")
                ), "replacement has 3 rows, data has 20")
                ## Confirm that refresh(ds) is unchanged
                expect_true(is.null(refresh(ds)$Two))
            })
            test_that("General input validation on addVariables", {
                expect_error(addVariables(ds,
                    VarDef(2, name="Two"),
                    5 ## Not a variable
                ), "Must supply VariableDefinitions")
                ## Confirm that refresh(ds) is unchanged
                expect_true(is.null(refresh(ds)$Two))
            })
            test_that("addVariables server error handling", {
                without_internet({
                    ## Add two expr vars (no GET on rows first)
                    expect_error(addVariables(ds,
                        VarDef(ds$v3 + 4, name="v3plus4"),
                        VarDef(ds$v3 + 5, name="v3plus5")
                    ), "The following variable definition\\(s\\) errored on upload: 1, 2")
                })
                ## Confirm that refresh(ds) is unchanged
                ds <- refresh(ds)
                expect_true(is.null(ds$v3plus4))
                expect_true(is.null(ds$v3plus5))
            })
        })
    })
}
