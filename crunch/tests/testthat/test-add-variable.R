context("Add a variable to a dataset")

test_that("toVariable parses R numerics", {
    expect_identical(toVariable(2L:4L, name="Numbers!", alias="num"),
        structure(list(values=2L:4L, type="numeric", name="Numbers!", alias="num"), class="VariableDefinition"))
    expect_equivalent(toVariable(2L:4L, name="Numbers!", alias="num"),
        list(values=2L:4L, type="numeric", name="Numbers!", alias="num"))
})
test_that("toVariable parses R characters", {
    expect_identical(toVariable(letters[1:3]),
        structure(list(values=c("a", "b", "c"), type="text"),
        class="VariableDefinition"))
})
test_that("toVariable parses factors", {
    expect_equivalent(toVariable(as.factor(rep(LETTERS[2:3], 3))),
        list(values=rep(1:2, 3), type="categorical", categories=list(
            list(id=1L, name="B", numeric_value=1L, missing=FALSE),
            list(id=2L, name="C", numeric_value=2L, missing=FALSE),
            list(id=-1L, name="No Data", numeric_value=NULL, missing=TRUE)
        ))) ## unclear why these aren't identical
    with(temp.options(crunch.max.categories=4), {
        expect_equivalent(toVariable(as.factor(letters[1:5])),
            list(values=c("a", "b", "c", "d", "e"), type="text"))
        expect_equivalent(toVariable(as.factor(letters[1:5]), name="v1"),
            list(values=c("a", "b", "c", "d", "e"), type="text", name="v1"))
    })
})

test_that("toVariable parses R Date class", {
    expect_equivalent(toVariable(as.Date(c("2014-12-16", "2014-12-17"))),
        list(values=c("2014-12-16", "2014-12-17"), type="datetime",
            resolution="D"))
})

test_that("toVariable handles POSIX datetimes", {
    skip("Investigate precision")
    numtime <- 1454238117.123 ## Note that it's off by 1ms below...
    expect_equivalent(toVariable(as.POSIXct(numtime, origin="1970-01-01", tz="UTC")),
        list(values="2016-01-31T03:01:57.122", type="datetime",
            resolution="ms"))
})

test_that("POSTNewVariable rejects invalid categories", {
    expect_error(POSTNewVariable("",
        list(type="categorical", name="bad ids",
            categories=list(
                list(id=-1L, name="B", numeric_value=1L, missing=FALSE),
                list(id=2L, name="C", numeric_value=2L, missing=FALSE),
                list(id=-1L, name="No Data", numeric_value=NULL, missing=TRUE)
            ))),
        "Invalid category ids: must be unique")
    expect_error(POSTNewVariable("",
        list(type="categorical", name="bad names",
            categories=list(
                list(id=1L, name="Name 1", numeric_value=1L, missing=FALSE),
                list(id=2L, name="Name 1", numeric_value=2L, missing=FALSE),
                list(id=-1L, name="No Data", numeric_value=NULL, missing=TRUE)
            ))),
        "Invalid category names: must be unique")
})

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            test_that("addVariable creates a new remote numeric variable", {
                ds <- addVariables(ds,
                    VariableDefinition(df$v3, name="New var", alias="newVar"))
                expect_true("newVar" %in% names(ds))
                nv <- ds$newVar
                expect_true(is.Numeric(nv))
                expect_true(is.Numeric(ds[['v3']]))
                expect_identical(as.vector(nv), as.vector(ds$v3))
            })
            test_that("addVariable creates text variables from character", {
                ds <- addVariables(ds,
                    VariableDefinition(df$v2, name="New var2", alias="newVar2"))
                expect_true("newVar2" %in% names(ds))
                nv <- ds$newVar2
                expect_true(is.Text(nv))
                expect_identical(as.vector(nv)[1:15],
                    as.vector(ds$v2)[1:15])
                    ## note that NAs aren't getting caught in the CSV importer
                    ## anymore, but they're right in the addVariable method
            })
            test_that("addVariable creates categorical from factor", {
                ds <- addVariables(ds,
                    VariableDefinition(df$v4, name="New var3", alias="newVar3"))
                expect_true("newVar3" %in% names(ds))
                nv <- ds$newVar3
                expect_true(is.Categorical(nv))
                expect_identical(as.vector(nv), as.vector(ds$v4))
            })
            test_that("addVariable creates datetime from Date", {
                ds <- addVariables(ds,
                    VariableDefinition(df$v5, name="New var4", alias="newVar4"))
                expect_true("newVar4" %in% names(ds))
                nv <- ds$newVar4
                expect_true(is.Datetime(nv))
                expect_identical(as.vector(nv), as.vector(ds$v5))
            })
            test_that("addVariable creates datetime from POSIXct", {
                skip("Can't support POSIXt until the app supports timezones")
                ds <- addVariables(ds, VariableDefinition(as.POSIXct(df$v5),
                    name="New var 5", alias="newVar5"))
                expect_true("newVar5" %in% names(ds))
                nv <- ds$newVar5
                expect_true(is.Datetime(nv))
                expect_identical(as.vector(nv), as.vector(ds$v5))
            })
            test_that("adding variable with duplicate name fails", {
                expect_error(addVariables(ds, VariableDefinition(df$v5,
                    name="New var4", alias="newVar4")),
                    "Variable with name: New var4 already exists")
            })
        })

        with(test.dataset(df), {
            test_that("assignment restrictions", {
                expect_error(ds[[2]] <- 1:20,
                    "Only character \\(name\\) indexing supported")
            })
            test_that("[[<- adds variables", {
                ds$newvariable <- 20:1
                expect_true(is.Numeric(ds$newvariable))
                expect_identical(mean(ds$newvariable), 10.5)
            })
            test_that("Variable lengths must match, in an R way", {
                expect_error(ds[['not valid']] <- 1:7,
                    "replacement has 7 rows, data has 20")
                ds$ok <- 1
                expect_identical(as.vector(ds$ok), rep(1, 20))
            })
        })

        with(test.dataset(df), {
            test_that("Categorical to R and back", {
                v4 <- as.vector(ds$v4)
                expect_identical(levels(v4), c("B", "C"))
                ds$v4a <- v4
                expect_equivalent(as.vector(ds$v4), as.vector(ds$v4a))
            })

            exclusion(ds) <- ds$v3 == 10
            test_that("Categorical to R and back with an exclusion", {
                v4b <- as.vector(ds$v4)
                expect_identical(levels(v4b), c("B", "C"))
                expect_identical(length(v4b), 19L)
                ds$v4b <- v4b
                expect_equivalent(as.vector(ds$v4b), as.vector(ds$v4a))
            })
        })
    })
}
