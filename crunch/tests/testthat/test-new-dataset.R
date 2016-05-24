context("Making a new dataset")

test_that("fake.csv is what we expect", {
    expect_identical(dim(testfile.df), c(20L, 6L))
})

if (run.integration.tests) {
    test_that("Source file cannot be uploaded if not logged in", {
        logout()
        expect_error(createSource(testfile.csv),
            "You must authenticate before making this request")
    })
    test_that("Dataset container object cannot be created if not logged in", {
        logout()
        expect_error(createDataset("testfile.csv"),
            "You must authenticate before making this request")
    })

    with(test.authentication, {
        ## New dataset by file upload method
        test_that("Source file can be uploaded if logged in", {
            expect_true(isTRUE(createSource(testfile.csv,
                status.handlers=list(`201`=function (response) TRUE))))
        })
        test_that("Dataset container object can be created if logged in", {
            with(test.dataset(), {
                expect_true(is.dataset(ds))
            })
        })
        test_that("Source can be added to Dataset", {
            source <- createSource(testfile.csv)
            with(test.dataset(), {
                ds <- try(addSourceToDataset(ds, source))
                expect_true(is.dataset(ds))
                expect_identical(nrow(ds), 20L)
                expect_identical(ncol(ds), 6L)
            })
        })
        test_that("newDatasetFromFile creates a dataset", {
            with(test.dataset(newDatasetFromFile(testfile.csv,
                                                name=uniqueDatasetName())), {
                expect_true(is.dataset(ds))
                expect_identical(nrow(ds), 20L)
                expect_identical(ncol(ds), 6L)
                expect_equivalent(mean(ds[[2]]), mean(testfile.df[[2]]))
            })
        })

        test_that("newDataset input validation", {
            expect_error(newDataset(NULL),
                "Can only make a Crunch dataset from a two-dimensional data")
            expect_error(newDataset(1:5),
                "Can only make a Crunch dataset from a two-dimensional data")
        })

        dsname <- uniqueDatasetName()
        with(test.dataset(newDatasetByColumn(df, name=dsname,
                                    description="a description")), {
            test_that("newDataset by addVariables", {
                expect_true(dsname %in% listDatasets())
                expect_true(is.dataset(ds))
                expect_identical(description(ds), "a description")
                expect_equivalent(mean(ds$v3), mean(df$v3))
                expect_identical(dim(ds), dim(df))
            })
        })

        with(test.dataset(newDatasetByColumn(df, name=uniqueDatasetName())), {
            test_that("Dataset-by-column variable types get set correctly",
                validImport(ds)
            )

            with(test.dataset(mrdf, "testmrdf"), {
                test_that("names() are the same and in the right order", {
                    expect_true(setequal(names(df), names(ds)))
                    expect_identical(names(df), names(ds))
                    expect_true(setequal(names(mrdf), names(testmrdf)))
                    expect_identical(names(mrdf), names(testmrdf))
                })
            })
        })

        with(test.dataset(suppressMessages(newDatasetByCSV(df,
                                            name=uniqueDatasetName()))), {
            test_that("newDataset via CSV + JSON", validImport(ds))
        })

        test_that("createWithMetadataAndFile using docs example", {
            with(test.dataset(newDatasetFromFixture("apidocs")), {
                expect_true(is.dataset(ds))
                expect_identical(name(ds), "Example dataset")
                expect_identical(names(categories(ds$q1)),
                    c("Cat", "Dog", "Bird", "Skipped", "Not Asked"))
            })
        })

        m <- fromJSON(file.path("dataset-fixtures", "apidocs.json"),
            simplifyVector=FALSE)

        test_that("Can create dataset with data in S3", {
            ds <- try(createWithMetadataAndFile(m,
                file="s3://public.testing.crunch.io/example-dataset.csv"))
            expect_true(is.dataset(ds))
            with(test.dataset(newDatasetFromFixture("apidocs")), as="ds2", {
                ## Compare to dataset imported from local file upload
                expect_identical(dim(ds), dim(ds2))
                expect_identical(as.vector(ds$q1), as.vector(ds2$q1))
                ## Could add more assertions
            })
            delete(ds)
        })

        test_that("Duplicate subvariables are forbidden", {
            m2 <- m
            ## Add a duplicate subvariable
            m2$body$table$metadata$allpets$subvariables[[4]] <- list(name="Another", alias="allpets_1")
            expect_error(createWithMetadataAndFile(m2,
                file.path("dataset-fixtures", "apidocs.csv")))
        })

        dsz <- try(suppressMessages(newDataset(df)))
        test_that("newDataset without specifying name grabs object name", {
            expect_true(is.dataset(dsz))
            expect_identical(name(dsz), "df")
            with(test.dataset(dsz), validImport(dsz))
        })

        test_that("Datasets can be deleted", {
            dsname <- uniqueDatasetName()
            testdf <- suppressMessages(newDataset(df, name=dsname))
            expect_true(dsname %in% listDatasets())
            expect_true(isTRUE(crDELETE(self(testdf),
                status.handlers=list(`204`=function (response) TRUE))))
            expect_false(dsname %in% listDatasets(refresh=TRUE))
        })
        test_that("Datasets can be deleted by S4 method", {
            dsname <- uniqueDatasetName()
            testdf <- suppressMessages(newDataset(df, name=dsname))
            expect_true(dsname %in% listDatasets())
            delete(testdf)
            expect_false(dsname %in% listDatasets())
        })
    })
}
