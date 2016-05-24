context("Variable catalog")

with_mock_HTTP({
    variables.catalog.url <- "/api/datasets/dataset1/variables.json"
    varblob <- crGET(variables.catalog.url)

    test_that("VariableCatalog instantiates from Shoji", {
        expect_true(inherits(VariableCatalog(varblob),
            "VariableCatalog"))
    })

    varcat <- VariableCatalog(varblob)
    order.url <- "/api/datasets/dataset1/variables/hierarchical.json"
    varorder <- VariableOrder(crGET(order.url))

    test_that("VariableCatalog index method", {
        expect_identical(names(index(varcat)), names(varcat@index))
        expect_identical(names(index(varcat)), names(varblob$index))
    })

    test_that("VariableCatalog has the right contents", {
        expect_true(all(grepl("/api/datasets/dataset1/variables",
            urls(varcat))))
        expect_identical(self(varcat), variables.catalog.url)
        expect_identical(entities(ordering(varcat)), entities(varorder))
    })

    test_that("active/hidden getters", {
        expect_identical(index(active(varcat)),
            index(varcat)[urls(ordering(varcat))])
        expect_equivalent(index(hidden(varcat)), list())
        index(varcat)[[2]]$discarded <- TRUE
        expect_true(inherits(active(varcat), "VariableCatalog"))
        expect_true(inherits(hidden(varcat), "VariableCatalog"))
        expect_identical(urls(active(varcat)),
            c("/api/datasets/dataset1/variables/gender.json",
            "/api/datasets/dataset1/variables/mymrset.json",
            "/api/datasets/dataset1/variables/textVar.json",
            "/api/datasets/dataset1/variables/starttime.json",
            "/api/datasets/dataset1/variables/catarray.json"))
        expect_identical(length(active(varcat)), 5L)
        expect_identical(urls(hidden(varcat)),
            "/api/datasets/dataset1/variables/birthyr.json")
        expect_identical(length(hidden(varcat)), 1L)
        expect_identical(length(varcat), 6L)
        expect_identical(active(hidden(varcat)), hidden(active(varcat)))
    })

    gender.url <- "/api/datasets/dataset1/variables/gender.json"
    test_that("Extract methods: character and numeric", {
        expect_true(inherits(varcat[[gender.url]], "VariableTuple"))
        expect_identical(varcat[[gender.url]]@body,
            index(varcat)[[gender.url]])
        expect_identical(index(varcat[2:3]), index(varcat)[2:3])
    })

    test_that("Extract methods: invalid input", {
        expect_error(varcat[[999]], "subscript out of bounds") ## base R
        expect_error(varcat[["asdf"]], "Subscript out of bounds: asdf")
        expect_error(varcat[[NA]], "Subscript out of bounds: NA")
        expect_error(varcat[999:1000], "Subscript out of bounds: 999:1000")
    })

    test_that("Extract methods: VariableOrder/Group", {
        ents <- c("/api/datasets/dataset1/variables/gender.json",
            "/api/datasets/dataset1/variables/mymrset.json")
        ord <- VariableOrder(VariableGroup("G1", entities=ents))
        expect_identical(names(varcat[ents]), c("Gender", "mymrset"))
        expect_identical(varcat[ord[[1]]], varcat[ents])
        expect_identical(varcat[ord], varcat[ents])
    })

    test_that("Construct Variable from Tuple", {
        expect_true(is.Categorical(CrunchVariable(varcat[[gender.url]])))
    })

    test_that("name and alias getters", {
        expect_identical(names(varcat)[1:3],
            c("Gender", "Birth Year", "starttime"))
        expect_identical(aliases(varcat)[1:2], c("gender", "birthyr"))
    })

    test_that("types getter", {
        expect_identical(types(varcat)[1:3],
            c("categorical", "numeric", "datetime"))
    })

    test_that("show method", {
        expect_identical(capture.output(print(varcat[1:3])),
            capture.output(print(data.frame(
                alias=c("gender", "birthyr", "starttime"),
                name=c("Gender", "Birth Year", "starttime"),
                type=c("categorical", "numeric", "datetime")
            ))))
    })
})

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            test_that("Can set descriptions (and doing so doesn't PUT order)", {
                with(temp.options(httpcache.log=""), {
                    expect_identical(descriptions(variables(ds)),
                        rep("", ncol(ds)))
                    logs <- capture.output(descriptions(variables(ds))[2:3] <- c("Des 1", "Des 2"))
                    expect_identical(descriptions(variables(ds))[1:4],
                        c("", "Des 1", "Des 2", ""))
                })
                expect_identical(length(logs), 2L) ## PATCH, DROP
            })
            test_that("Can set names and aliases", {
                n <- names(df)
                expect_identical(names(variables(ds)), n)
                expect_identical(aliases(variables(ds)), n)
                names(variables(ds))[2:3] <- c("two", "three")
                n2 <- n
                n2[2:3] <- c("two", "three")
                expect_identical(names(variables(ds)), n2)
                expect_identical(names(variables(refresh(ds))), n2)
                n3 <- n
                n3[c(2,4)] <- c("due", "quattro")
                aliases(variables(ds))[c(2,4)] <- c("due", "quattro")
                expect_identical(aliases(variables(ds)), n3)
                expect_identical(aliases(variables(refresh(ds))), n3)
            })

            test_that("Can [<- with VariableGroup/Order", {
                names(variables(ds))[2:3] <- c("two", "three")
                ord <- VariableOrder(VariableGroup("a group", entities=ds[2:3]))
                expect_identical(names(variables(ds)[ord]),
                    c("two", "three"))
                try(names(variables(ds)[ord[[1]]]) <- c("TWO", "Three"))
                expect_identical(names(variables(ds)[ord]),
                    c("TWO", "Three"))
                try(names(variables(ds)[ord]) <- c("2", "3"))
                expect_identical(names(variables(ds)[ord]),
                    c("2", "3"))
            })
        })
    })
}
