context("Filtering datasets and variables in the R session")

with_mock_HTTP({
    ds <- loadDataset("test ds")
    ds2 <- ds[ds$gender == "Male", ]
    ds3 <- ds2[ds$birthyr > 1981, ]

    test_that("A clean dataset has NULL activeFilter", {
        expect_identical(activeFilter(ds), NULL)
    })
    test_that("Getting a variable from a clean dataset has NULL activeFilter", {
        expect_identical(activeFilter(ds$gender), NULL)
    })

    test_that("[ method on dataset adds an active filter", {
        expect_identical(activeFilter(ds2), ds$gender == "Male")
    })
    test_that("Active filter persists on refreshing dataset", {
        expect_identical(activeFilter(refresh(ds2)), ds$gender == "Male")
    })
    test_that("Further [ on a filtered dataset ands the filters together", {
        expect_identical(activeFilter(ds3),
            ds$gender == "Male" & ds$birthyr > 1981)
    })

    test_that("subset method for dataset does the same", {
        expect_identical(subset(ds, ds$gender == "Male"), ds2)
    })

    test_that("Variables extracted from a filtered dataset are also filtered", {
        expect_identical(activeFilter(ds2$birthyr), ds$gender == "Male")
        skip("Not really supporting this correctly")
        expect_identical(ds[ds$gender == "Male", "birthyr"],
            ds$gender == "Male")
    })

    test_that("Subvariables extracted from a filtered array are also filtered", {
        mr <- ds2$mymrset
        expect_true(is.Multiple(mr))
        expect_identical(activeFilter(mr), ds$gender == "Male")
        sv1 <- subvariables(mr)[[1]]
        expect_true(is.Categorical(sv1))
        expect_identical(activeFilter(sv1), ds$gender == "Male")
    })

    test_that("Active filter persists on refreshing variable", {
        expect_identical(activeFilter(refresh(ds2$birthyr)),
            ds$gender == "Male")
    })

    test_that("Getting weight variable from filtered dataset is filtered", {
        ds4 <- ds2
        ds4@body$weight <- "/api/datasets/dataset1/variables/starttime.json"
        expect_identical(weight(ds4), ds4$starttime)
        expect_identical(activeFilter(weight(ds4)), ds$gender == "Male")
    })

    test_that("activeFilter from filtered CrunchVariable", {
        expect_identical(activeFilter(ds$birthyr), NULL)
        expect_identical(activeFilter(ds2$birthyr),
            activeFilter(ds$birthyr[ds$gender == "Male"]))
        expect_identical(activeFilter(ds$birthyr[ds$gender == "Male"]),
            ds$gender == "Male")
    })

    test_that("activeFilter from CrunchExpr", {
        ## Need the @expression here because CrunchExpr@filter is list
        expect_identical(activeFilter((ds$birthyr - ds$starttime)[ds$gender == "Male"]),
            (ds$gender == "Male")@expression)
    })
    test_that("activeFilter passes across operations among vars/exprs", {
        skip("TODO: filtered variables/exprs should pass their filters along")
        expect_identical(activeFilter(ds$birthyr[ds$gender == "Male"] - ds$starttime[ds$gender == "Male"]),
            ds$gender == "Male")
    })
})

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            ds2 <- ds[ds$v4 == "C",]
            ds2b <- ds[ds$v4 != "B",]
            ds3 <- ds[ds$v3 > 11,]
            ds4 <- ds[is.na(ds$v1),]

            test_that("filtered dim", {
                expect_identical(dim(ds2), c(10L, 6L))
                expect_identical(dim(ds2b), c(10L, 6L))
                expect_identical(dim(ds3), c(16L, 6L))
                expect_identical(dim(ds4), c(5L, 6L))
            })

            test_that("activeFilter appears in print method for dataset", {
                expect_output(ds3, "Filtered by v3 > 11")
                expect_false(any(grepl("Filtered by",
                    capture.output(print(ds)))))
            })

            test_that("Filtered variables return filtered values from as.vector", {
                expect_identical(as.vector(ds2$v3),
                    c(9, 11, 13, 15, 17, 19, 21, 23, 25, 27))
                expect_identical(as.vector(ds2b$v3),
                    c(9, 11, 13, 15, 17, 19, 21, 23, 25, 27))
                expect_identical(as.vector(ds3$v3),
                    as.numeric(12:27))
                expect_identical(as.vector(ds4$v3),
                    as.numeric(8:12))
            })

            test_that("activeFilter appears in print method for variables", {
                expect_output(ds3$v3, "Filtered by v3 > 11")
            })

            test_that("as.data.frame when filtered", {
                df2 <- as.data.frame(ds2)
                expect_identical(df2$v3,
                    c(9, 11, 13, 15, 17, 19, 21, 23, 25, 27))
                expect_equivalent(as.data.frame(ds2[,c("v3", "v4")],
                    force=TRUE),
                    df[df$v4 == "C", c("v3", "v4")])
                df3 <- as.data.frame(ds3)
                expect_equivalent(df3$v3, 12:27)
            })

            test_that("filtered cubing", {
                expect_equivalent(as.array(crtabs(~ v4,
                    data=ds[ds$v4 == "C",])),
                    array(c(0, 10), dim=2L, dimnames=list(v4=c("B", "C"))))
                expect_equivalent(as.array(crtabs(~ v4,
                    data=ds4)),
                    array(c(3, 2), dim=2L, dimnames=list(v4=c("B", "C"))))
            })

            test_that("filtered updating", {
                skip("TODO")
            })
        })
    })
}
