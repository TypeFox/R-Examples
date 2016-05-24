context("Deleting rows of a dataset")

with_mock_HTTP({
    ds <- loadDataset("test ds")
    test_that("dropRows generates the right request", {
        expect_error(dropRows(ds, ds$gender == "Male"),
            paste0('POST /api/datasets/dataset1/table.json ',
            '{"command":"delete","filter":{"function":"==",',
            '"args":[{"variable":"/api/datasets/dataset1/variables/gender.json"},',
            '{"value":1,"type":{"function":"typeof","args":[',
            '{"variable":"/api/datasets/dataset1/variables/gender.json"}]}}]}}'),
            fixed=TRUE)
    })
})

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            test_that("dropRows really removes rows", {
                try(ds <- dropRows(ds, ds$v4 == "C"))
                expect_identical(dim(ds), c(10L, ncol(df)))
                expect_identical(as.vector(ds$v4, mode="id"), rep(1, 10))
                expect_identical(as.vector(ds$v3), seq(8, 26, 2))
            })
        })
        with(test.dataset(df), {
            exclusion(ds) <- ds$v4 == "B"
            test_that("dropRows correctly drops with an exclusion applied", {
                expect_identical(nrow(ds), 10L)
                try(ds <- dropRows(ds, ds$v3 > 10 & ds$v3 <= 15))
                expect_identical(dim(ds), c(7L, ncol(df)))
                exclusion(ds) <- NULL
                expect_identical(nrow(ds), 15L)
            })
        })
    })
}
