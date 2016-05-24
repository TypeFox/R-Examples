context("Batch catalog")

with_mock_HTTP({
    ds <- loadDataset("test ds")
    test_that("batches method", {
        expect_true(inherits(batches(ds), "BatchCatalog"))
        expect_identical(length(batches(ds)), 3L)
        expect_identical(urls(batches(ds)),
            c("/api/datasets/dataset1/batches/0.json",
            "/api/datasets/dataset1/batches/2.json",
            "/api/datasets/dataset1/batches/3.json"))
    })

    test_that("imported/pending", {
        expect_identical(urls(imported(batches(ds))),
            c("/api/datasets/dataset1/batches/0.json",
            "/api/datasets/dataset1/batches/2.json"))
        expect_identical(urls(pending(batches(ds))),
            "/api/datasets/dataset1/batches/3.json")
    })

    test_that("show method for batch catalog", {
        expect_identical(capture.output(print(batches(ds))),
            capture.output(print(data.frame(id=c(0, 2, 3),
            status=c("imported", "imported", "conflict")))))
    })
})
