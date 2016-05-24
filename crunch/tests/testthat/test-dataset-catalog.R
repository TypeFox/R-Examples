context("Dataset catalog")

with_mock_HTTP({
    cr <- session()
    datcat <- cr$datasets

    test_that("DatasetCatalog instantiates from Shoji", {
        expect_true(inherits(datcat, "DatasetCatalog"))
    })

    test_that("DatasetCatalog has the right contents", {
        expect_identical(urls(datcat),
            c("/api/datasets/dataset3.json",
              "/api/datasets/dataset2.json",
              "/api/datasets/dataset1.json")) ## C sorting on names
        expect_identical(self(datcat),
            "/api/datasets.json")
    })

    test_that("active/archived getters", {
        ## NOTE: deferring the "shared" collection
        expect_true(inherits(active(datcat), "DatasetCatalog"))
        expect_true(inherits(archived(datcat), "DatasetCatalog"))
        expect_identical(index(active(datcat)), index(datcat))
        expect_equivalent(index(archived(datcat)), list())
        index(datcat)[[which(names(datcat) == "an archived dataset")]]$archived <- TRUE
        expect_true(inherits(active(datcat), "DatasetCatalog"))
        expect_true(inherits(archived(datcat), "DatasetCatalog"))
        expect_identical(urls(active(datcat)),
            c("/api/datasets/dataset3.json", "/api/datasets/dataset1.json"))
        expect_identical(length(active(datcat)), 2L)
        expect_identical(urls(archived(datcat)),
            "/api/datasets/dataset2.json")
        expect_identical(length(archived(datcat)), 1L)
        expect_identical(length(datcat), 3L)
        expect_identical(active(archived(datcat)), archived(active(datcat)))
    })

    test_that("Extract methods", {
        expect_true(inherits(datcat[["/api/datasets/dataset1.json"]],
            "DatasetTuple"))
        expect_identical(datcat[["/api/datasets/dataset1.json"]]@body,
            index(datcat)[["/api/datasets/dataset1.json"]])
        expect_identical(index(datcat[2:3]), index(datcat)[2:3])
        expect_error(datcat[[500]], "subscript out of bounds")
    })

    test_that("names", {
        expect_identical(names(datcat),
            c("ECON.sav", "an archived dataset", "test ds"))
    })

    test_that("entity method for tuple", {
        expect_true(is.dataset(entity(datcat[["/api/datasets/dataset1.json"]])))
    })
})
