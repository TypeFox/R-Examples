context("Update variables with NAs")

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(newDatasetFromFixture("apidocs")), {
            test_that("Insert NA into numeric", {
                expect_equivalent(as.vector(ds$ndogs[1:5]),
                    c(1, NA, 2, 3, 1))
                ds$ndogs[ds$ndogs > 2] <- NA
                expect_equivalent(as.vector(ds$ndogs[1:5]),
                    c(1, NA, 2, NA, 1))
            })
            test_that("Insert NA into categorical", {
                expect_equivalent(as.character(as.vector(ds$q1[1:5])),
                    c(NA, "Cat", NA, "Dog", "Dog"))
                expect_equivalent(as.vector(ds$q1[1:5], mode="id"),
                    c(8, 1, 9, 2, 2))
                ds$q1[4] <- NA
                expect_equivalent(as.character(as.vector(ds$q1[1:5])),
                    c(NA, "Cat", NA, NA, "Dog"))
                expect_equivalent(as.vector(ds$q1[1:5], mode="id"),
                    c(8, 1, 9, -1, 2))
                expect_true(-1 %in% ids(categories(ds$q1)))
            })
            test_that("Insert NA into datetime", {
                expect_equivalent(as.vector(ds$wave[1:5]),
                    rep(as.Date("2014-12-01"), 5))
                skip("See https://www.pivotaltracker.com/story/show/103531536")
                ds$wave[4] <- NA
                expect_equivalent(as.vector(ds$wave[1:5]),
                    c(rep(as.Date("2014-12-01"), 3), NA, as.Date("2014-12-01")))
            })
            test_that("Insert NA into text", {
                expect_equivalent(as.vector(ds$q3[1:3]),
                    c("Jasmine", "Clyde", "Geoffrey"))
                ds$q3[2] <- NA
                expect_equivalent(as.vector(ds$q3[1:3]),
                    c("Jasmine", NA, "Geoffrey"))
            })
            test_that("Insert NA into multiple response", {
                expect_equivalent(as.vector(ds$allpets$Cat, mode="id")[1:5],
                    c(1, 9, 1, 1, 9))
                expect_equivalent(as.vector(ds$allpets$Dog, mode="id")[1:5],
                    c(8, 1, 9, 9, 9))
                expect_equivalent(as.vector(ds$allpets$Bird, mode="id")[1:5],
                    c(8, 2, 8, 8, 8))
                ds$allpets[2] <- NA
                expect_equivalent(as.vector(ds$allpets$Cat, mode="id")[1:5],
                    c(1, -1, 1, 1, 9))
                expect_equivalent(as.vector(ds$allpets$Dog, mode="id")[1:5],
                    c(8, -1, 9, 9, 9))
                expect_equivalent(as.vector(ds$allpets$Bird, mode="id")[1:5],
                    c(8, -1, 8, 8, 8))
            })
            test_that("Insert NA into categorical array", {
                expect_equivalent(as.vector(ds$petloc$Home, mode="id")[1:5],
                    c(8, 2, 9, 9, 1))
                expect_equivalent(as.vector(ds$petloc$Work, mode="id")[1:5],
                    c(9, 3, 3, 2, 2))
                ds$petloc[3] <- NA
                expect_equivalent(as.vector(ds$petloc$Home, mode="id")[1:5],
                    c(8, 2, -1, 9, 1))
                expect_equivalent(as.vector(ds$petloc$Work, mode="id")[1:5],
                    c(9, 3, -1, 2, 2))
            })
        })
        with(test.dataset(newDatasetFromFixture("apidocs")), {
            test_that("Insert values including NA into numeric", {
                expect_equivalent(as.vector(ds$ndogs[1:5]),
                    c(1, NA, 2, 3, 1))
                ds$ndogs[2:4] <- c(1, NA, 2)
                expect_equivalent(as.vector(ds$ndogs[1:5]),
                    c(1, 1, NA, 2, 1))
            })
            test_that("Insert values including NA into categorical", {
                expect_equivalent(as.character(as.vector(ds$q1[1:5])),
                    c(NA, "Cat", NA, "Dog", "Dog"))
                ds$q1[2:4] <- c(NA, "Cat", "Cat")
                expect_equivalent(as.character(as.vector(ds$q1[1:5])),
                    c(NA, NA, "Cat", "Cat", "Dog"))
            })
            test_that("Insert values including NA into datetime", {
                expect_equivalent(as.vector(ds$wave[1:5]),
                    rep(as.Date("2014-12-01"), 5))
                ds$wave[2:4] <- as.Date(c("2014-12-15", NA, "2014-11-01"))
                expect_equivalent(as.Date(as.vector(ds$wave[1:5])), ## it's POSIXt
                    as.Date(c("2014-12-01", "2014-12-15", NA, "2014-11-01", "2014-12-01")))
            })
            test_that("Insert values including NA into text", {
                expect_equivalent(as.vector(ds$q3[1:3]),
                    c("Jasmine", "Clyde", "Geoffrey"))
                ds$q3[2:3] <- c(NA, "Jeff")
                expect_equivalent(as.vector(ds$q3[1:3]),
                    c("Jasmine", NA, "Jeff"))
            })
            test_that("Insert values including NA into multiple response", {

            })
            test_that("Insert values including NA into categorical array", {

            })
        })

        with(test.dataset(df), {
            test_that("Can set missing rules", {
                expect_error(is.na(ds$v5) <- ds$v4 == "B",
                    "is.na<- not yet supported")
                skip("mark missing as function seems not to work")
                try(is.na(ds$v3) <- ds$v3 >= 10 & ds$v3 < 13)
                expect_equivalent(as.vector(ds$v3),
                    c(8, 9, rep(NA, 3), 13:27))
                try(is.na(ds$v5) <- ds$v4 == "B")
                expect_identical(sum(is.na(as.vector(ds$v5))), 10L)
                try(is.na(ds$v3) <- ds$v4 %in% "C")
                print(as.vector(ds$v3))
                expect_identical(sum(is.na(as.vector(ds$v3))), 10L)
            })
        })
    })
}
