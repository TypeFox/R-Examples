context("Crosstabbing")

cubedf <- df
cubedf$v7 <- as.factor(c(rep("C", 10), rep("D", 5), rep("E", 5)))
cubedf$v8 <- as.Date(0:1, origin="1955-11-05")

test_that("bin CrunchExpr", {
    x <- list(variable="test") ## "ZCL"
    expect_true(inherits(bin(x), "CrunchExpr"))
    expect_identical(zcl(bin(x)),
        list(`function`="bin", args=list(list(variable="test"))))
})

test_that("rollup CrunchExpr from zcl variable", {
    x <- list(variable="test") ## "ZCL"
    expect_true(inherits(rollup(x), "CrunchExpr"))
    expect_identical(zcl(rollup(x)),
        list(`function`="rollup", args=list(list(variable="test"),
            list(value=NULL))))

    expect_true(inherits(rollup(x, resolution="Y"), "CrunchExpr"))
    expect_identical(zcl(rollup(x, resolution="Y")),
        list(`function`="rollup", args=list(list(variable="test"),
            list(value="Y"))))
})

with_mock_HTTP({
    ds <- loadDataset("test ds")
    v <- ds$starttime

    test_that("rollup CrunchExpr from DatetimeVariable", {
        expect_true(inherits(rollup(v), "CrunchExpr"))
        expect_identical(zcl(rollup(v)),
            list(`function`="rollup",
                args=list(list(variable="/api/datasets/dataset1/variables/starttime.json"),
                list(value="s"))))
        expect_true(inherits(rollup(v, resolution="Y"), "CrunchExpr"))
        expect_identical(zcl(rollup(v, resolution="Y")),
            list(`function`="rollup",
                args=list(list(variable="/api/datasets/dataset1/variables/starttime.json"),
                list(value="Y"))))
        expect_true(inherits(rollup(v, resolution=NULL), "CrunchExpr"))
        expect_identical(zcl(rollup(v, resolution=NULL)),
            list(`function`="rollup",
                args=list(list(variable="/api/datasets/dataset1/variables/starttime.json"),
                list(value=NULL))))
    })
})

test_that("rollup resolution validation", {
    expect_error(rollup("a", resolution="Invalid"),
        " is invalid. Valid values are NULL, ")
    expect_error(rollup("a", resolution=42),
        " is invalid. Valid values are NULL, ")
})

adims <- CubeDims(list(
    v4=list(name=c("B", "C"), any.or.none=rep(FALSE, 2),
        missing=rep(FALSE, 2)),
    v7=list(name=c("C", "D", "E", "No Data"), any.or.none=rep(FALSE, 4),
        missing=c(rep(FALSE, 3), TRUE))))
a1 <- CrunchCube(arrays=list("count"=array(c(
        8, 6,
        3, 2,
        2, 3,
        0, 0), dim=c(2L, 4L))),
    dims=adims)
#    v7
# v4  C D E No Data
#   B 8 3 2 0
#   C 6 2 3 0

test_that("simple margin.table", {
    expect_equivalent(margin.table(a1, 1), margin.table(a1@arrays[[1]], 1))
    expect_identical(margin.table(a1, 1),
        array(c(13, 11), dim=2L, dimnames=list(v4=c("B", "C"))))
    expect_identical(margin.table(a1, 2),
        array(c(14, 5, 5), dim=3L, dimnames=list(v7=LETTERS[3:5])))
    expect_equivalent(margin.table(a1), margin.table(a1@arrays[[1]]))
    expect_identical(margin.table(a1), 24)
})

test_that("margin.table with any/none", {
    a2 <- a1
    a2@dims[[2]]$any.or.none[1] <- TRUE ## "C"
    expect_identical(margin.table(a2, 1),
        array(c(8, 6), dim=2L, dimnames=list(v4=c("B", "C"))))
    expect_identical(margin.table(a2, 2),
        array(c(5, 5), dim=2L, dimnames=list(v7=LETTERS[4:5])))
    expect_identical(margin.table(a2), 14)
})

test_that("margin.table with missing", {
    a2 <- a1
    a2@dims[[2]]$missing[2] <- TRUE ## "D"
    expect_identical(a2@useNA, "no") ## The default.
    expect_identical(margin.table(a2, 1),
        array(c(10, 9), dim=2L, dimnames=list(v4=c("B", "C"))))
    expect_identical(margin.table(a2, 2),
        array(c(14, 5), dim=2L, dimnames=list(v7=c("C", "E"))))
    expect_identical(margin.table(a2), 19)

    a2@useNA <- "ifany"
    ## Should be the same as first tests
    expect_identical(margin.table(a2, 1),
        array(c(13, 11), dim=2L, dimnames=list(v4=c("B", "C"))))
    expect_identical(margin.table(a2, 2),
        array(c(14, 5, 5), dim=3L, dimnames=list(v7=LETTERS[3:5])))
    expect_identical(margin.table(a2), 24)

    a2@useNA <- "always"
    expect_identical(margin.table(a2, 1),
        array(c(13, 11), dim=2L, dimnames=list(v4=c("B", "C"))))
    expect_identical(margin.table(a2, 2),
        array(c(14, 5, 5, 0), dim=4L,
            dimnames=list(v7=c(LETTERS[3:5], "No Data"))))
    expect_identical(margin.table(a2), 24)
})

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(cubedf), {
            test_that("cubedf setup", {
                expect_identical(names(categories(ds$v7)),
                    c("C", "D", "E", "No Data"))
            })
            test_that("We can get a univariate categorical cube", {
                kube <- try(crtabs(~ v7, data=ds))
                expect_true(inherits(kube, "CrunchCube"))
                expect_equivalent(as.array(kube),
                    array(c(10, 5, 5), dim=c(3L),
                        dimnames=list(v7=LETTERS[3:5])))
                ## Not sure why not identical, str makes them look the same
            })

            test_that("We can get a bivariate categorical cube", {
                kube <- try(crtabs(~ v4 + v7, data=ds))
                expect_true(inherits(kube, "CrunchCube"))
                expect_identical(as.array(kube),
                    array(c(5, 5, 3, 2, 2, 3), dim=c(2L, 3L),
                        dimnames=list(v4=c("B", "C"), v7=LETTERS[3:5])))
            })

            is.na(categories(ds$v7)) <- "D"
            test_that("useNA on univariate cube", {
                expect_equivalent(as.array(crtabs(~ v7, data=ds)),
                    array(c(10, 5), dim=2L, dimnames=list(v7=c("C", "E"))))
                expect_equivalent(as.array(crtabs(~ v7, data=ds,
                    useNA="ifany")),
                    array(c(10, 5, 5), dim=c(3L),
                        dimnames=list(v7=LETTERS[3:5])))
                expect_equivalent(as.array(crtabs(~ v7, data=ds,
                    useNA="always")),
                    array(c(10, 5, 5, 0), dim=c(4L),
                        dimnames=list(v7=c(LETTERS[3:5], "No Data"))))
            })
            test_that("useNA on bivariate cube", {
                expect_equivalent(as.array(crtabs(~ v4 + v7, data=ds)),
                    array(c(5, 5, 2, 3), dim=c(2L, 2L),
                        dimnames=list(v4=c("B", "C"), v7=c("C", "E"))))
                expect_equivalent(as.array(crtabs(~ v4 + v7, data=ds,
                    useNA="ifany")),
                    array(c(5, 5, 3, 2, 2, 3), dim=c(2L, 3L),
                        dimnames=list(v4=c("B", "C"), v7=LETTERS[3:5])))
                expect_equivalent(as.array(crtabs(~ v4 + v7, data=ds,
                    useNA="always")),
                    array(c(5, 5, 0,
                            3, 2, 0,
                            2, 3, 0,
                            0, 0, 0), dim=c(3L, 4L),
                        dimnames=list(v4=c("B", "C", "No Data"),
                                    v7=c(LETTERS[3:5], "No Data"))))
            })

            test_that("univariate datetime cube", {
                kube <- try(crtabs(~ v8, data=ds))
                expect_true(inherits(kube, "CrunchCube"))
                expect_equivalent(as.array(kube),
                    array(c(10, 10), dim=c(2L),
                        dimnames=list(v8=c("1955-11-05", "1955-11-06"))))
            })
            test_that("bivariate cube with datetime", {
                expect_equivalent(as.array(crtabs(~ v8 + v7, data=ds)),
                    array(c(5, 5, 2, 3), dim=c(2L, 2L),
                        dimnames=list(v8=c("1955-11-05", "1955-11-06"),
                        v7=c("C", "E"))))
                expect_equivalent(as.array(crtabs(~ v8 + v7, data=ds,
                    useNA="ifany")),
                    array(c(5, 5, 3, 2, 2, 3), dim=c(2L, 3L),
                        dimnames=list(v8=c("1955-11-05", "1955-11-06"),
                        v7=LETTERS[3:5])))
                expect_equivalent(as.array(crtabs(~ v8 + v7, data=ds,
                    useNA="always")),
                    array(c(5, 5,
                            3, 2,
                            2, 3,
                            0, 0), dim=c(2L, 4L),
                        dimnames=list(v8=c("1955-11-05", "1955-11-06"),
                        v7=c(LETTERS[3:5], "No Data"))))
            })

            test_that("datetime rollup cubes", {
                ## Default rollup resolution for this should be same as
                ## its resolution, given the date range
                expect_equivalent(as.array(crtabs(~ rollup(v8) + v7,
                    data=ds)),
                    as.array(crtabs(~ v8 + v7, data=ds)))
                expect_equivalent(as.array(crtabs(~ rollup(v8, "M") + v7,
                    data=ds)),
                    array(c(10, 5), dim=c(1L, 2L),
                        dimnames=list(v8=c("1955-11"),
                        v7=c("C", "E"))))
                expect_equivalent(as.array(crtabs(~ rollup(v8, "Y") + v7,
                    data=ds)),
                    array(c(10, 5), dim=c(1L, 2L),
                        dimnames=list(v8=c("1955"),
                        v7=c("C", "E"))))
            })

            test_that("univariate cube with binned numeric", {
                kube <- try(crtabs(~ bin(v3), data=ds))
                expect_true(inherits(kube, "CrunchCube"))
                expect_equivalent(as.array(kube),
                    array(c(2, 5, 5, 5, 3), dim=c(5L),
                        dimnames=list(v3=c("5-10", "10-15", "15-20", "20-25",
                        "25-30"))))
            })
            test_that("bivariate cube with binned numeric", {
                expect_equivalent(as.array(crtabs(~ bin(v3) + v7, data=ds)),
                    array(c(2, 5, 3, 0, 0,
                            0, 0, 0, 2, 3), dim=c(5L, 2L),
                        dimnames=list(v3=c("5-10", "10-15", "15-20", "20-25",
                        "25-30"),
                        v7=c("C", "E"))))
                expect_equivalent(as.array(crtabs(~ bin(v3) + v7, data=ds,
                    useNA="ifany")),
                    array(c(2, 5, 3, 0, 0,
                            0, 0, 2, 3, 0,
                            0, 0, 0, 2, 3), dim=c(5L, 3L),
                        dimnames=list(v3=c("5-10", "10-15", "15-20", "20-25",
                        "25-30"),
                        v7=LETTERS[3:5])))
                expect_equivalent(as.array(crtabs(~ bin(v3) + v7, data=ds,
                    useNA="always")),
                    array(c(2, 5, 3, 0, 0,
                            0, 0, 2, 3, 0,
                            0, 0, 0, 2, 3,
                            0, 0, 0, 0, 0), dim=c(5L, 4L),
                        dimnames=list(v3=c("5-10", "10-15", "15-20", "20-25",
                        "25-30"),
                        v7=c(LETTERS[3:5], "No Data"))))
            })
            test_that("unbinned numeric", {
                expect_equivalent(as.array(crtabs(~ v1, data=ds)),
                    array(rep(1, 15), dim=15L, dimnames=list(v1=df$v1[6:20])))
                expect_equivalent(as.array(crtabs(~ v1, data=ds,
                    useNA="ifany")),
                    array(c(rep(1, 15), 5), dim=16L,
                        dimnames=list(v1=c(df$v1[6:20], "<NA>"))))
            })

            test_that("Weighted cubes", {
                weight(ds) <- ds$v3
                expect_equivalent(as.array(crtabs(~ v8 + v7, data=ds)),
                    array(c(60, 65, 50, 75), dim=c(2L, 2L),
                        dimnames=list(v8=c("1955-11-05", "1955-11-06"),
                        v7=c("C", "E"))))
                expect_equivalent(as.array(crtabs(~ v8 + v7, data=ds,
                    weight=NULL)),
                    array(c(5, 5, 2, 3), dim=c(2L, 2L),
                        dimnames=list(v8=c("1955-11-05", "1955-11-06"),
                        v7=c("C", "E"))))
                weight(ds) <- NULL
                expect_equivalent(as.array(crtabs(~ v8 + v7, data=ds)),
                    array(c(5, 5, 2, 3), dim=c(2L, 2L),
                        dimnames=list(v8=c("1955-11-05", "1955-11-06"),
                        v7=c("C", "E"))))
                expect_equivalent(as.array(crtabs(~ v8 + v7, data=ds,
                    weight=ds$v3)),
                    array(c(60, 65, 50, 75), dim=c(2L, 2L),
                        dimnames=list(v8=c("1955-11-05", "1955-11-06"),
                        v7=c("C", "E"))))
            })

            test_that("Numeric aggregates", {
                expect_equivalent(as.array(crtabs(mean(v3) ~ v8 + v7,
                    data=ds)),
                    array(c(12, 13, 25, 25), dim=c(2L, 2L),
                        dimnames=list(v8=c("1955-11-05", "1955-11-06"),
                        v7=c("C", "E"))))
                expect_equivalent(as.array(crtabs(sum(v3) ~ v8 + v7,
                    data=ds)),
                    array(c(60, 65, 50, 75), dim=c(2L, 2L),
                        dimnames=list(v8=c("1955-11-05", "1955-11-06"),
                        v7=c("C", "E"))))
                expect_equivalent(as.array(crtabs(min(v3) ~ v8 + v7,
                    data=ds)),
                    array(c(8, 9, 24, 23), dim=c(2L, 2L),
                        dimnames=list(v8=c("1955-11-05", "1955-11-06"),
                        v7=c("C", "E"))))
            })

            test_that("Missing values in cubes", {
                expect_equivalent(round(as.array(crtabs(sd(v3) ~ bin(v3) + v7,
                    data=ds)), 3),
                    array(c(0.707, 1.581, 1, NaN, NaN,
                            NaN, NaN, NaN, 0.707, 1), dim=c(5L, 2L),
                        dimnames=list(v3=c("5-10", "10-15", "15-20", "20-25",
                        "25-30"),
                        v7=c("C", "E"))))
            })

            test_that("round cubes", {
                expect_equivalent(round(crtabs(sd(v3) ~ bin(v3) + v7,
                    data=ds), 3),
                    array(c(0.707, 1.581, 1, NaN, NaN,
                            NaN, NaN, NaN, 0.707, 1), dim=c(5L, 2L),
                        dimnames=list(v3=c("5-10", "10-15", "15-20", "20-25",
                        "25-30"),
                        v7=c("C", "E"))))
            })

            test_that("Cube with variables and R objects", {
                skip("(400) Bad Request: No such category id: '1'.")
                d4 <- cubedf$v4
                expect_equivalent(as.array(crtabs(~ d4 + v7, data=ds)),
                    array(c(5, 5, 2, 3), dim=c(2L, 2L),
                        dimnames=list(v4=c("B", "C"), v7=c("C", "E"))))
            })

            test_that("Cube with transformations", {
                expect_equivalent(as.array(crtabs(~ bin(v3 + 5), data=ds)),
                    array(c(2, 5, 5, 5, 3), dim=c(5L),
                        dimnames=list(v3=c("10-15", "15-20", "20-25",
                        "25-30", "30-35"))))
            })

            test_that("prop.table on univariate cube", {
                expect_equivalent(prop.table(crtabs(~ bin(v3 + 5), data=ds)),
                    array(c(2, 5, 5, 5, 3)/20, dim=c(5L),
                        dimnames=list(v3=c("10-15", "15-20", "20-25",
                        "25-30", "30-35"))))
            })

            test_that("prop.table on crosstab", {
                expect_equivalent(prop.table(crtabs(~ bin(v3) + v7, data=ds)),
                    array(c(2, 5, 3, 0, 0,
                            0, 0, 0, 2, 3)/15, dim=c(5L, 2L),
                        dimnames=list(v3=c("5-10", "10-15", "15-20", "20-25",
                        "25-30"),
                        v7=c("C", "E"))))
                expect_equivalent(prop.table(crtabs(~ bin(v3) + v7, data=ds),
                    margin=1),
                    array(c(1, 1, 1, 0, 0,
                            0, 0, 0, 1, 1), dim=c(5L, 2L),
                        dimnames=list(v3=c("5-10", "10-15", "15-20", "20-25",
                        "25-30"),
                        v7=c("C", "E"))))
                expect_equivalent(prop.table(crtabs(~ bin(v3) + v7, data=ds),
                    margin=2),
                    array(c(.2, .5, .3, 0, 0,
                            0, 0, 0, .4, .6), dim=c(5L, 2L),
                        dimnames=list(v3=c("5-10", "10-15", "15-20", "20-25",
                        "25-30"),
                        v7=c("C", "E"))))
            })

            test_that("Univariate stats", {
                expect_equivalent(as.array(crtabs(mean(v3) ~ 1, data=ds)),
                    17.5)
            })
        })
    })
}
