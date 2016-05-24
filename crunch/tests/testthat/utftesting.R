test_that("encoding + JSON reads correctly", {
    s <- iconv("aided_follow_grid:ElCorteInglés", to="UTF-8")
    expect_identical(Encoding(s), "UTF-8")
    expect_true(grepl("Inglés", s))
    sj <- toJSON(s)
    expect_true(grepl("Inglés", sj))
    s2 <- fromJSON(sj)
    expect_identical(s2, s)
    expect_identical(fromJSON("utf-test.json"), "Budějovický Budvar")
})

with_mock_HTTP({
    ds <- loadDataset("test ds")
    test_that("Reading UTF in tests", {
        expect_identical(description(ds$textVar), "Budějovický Budvar")
    })
})

with(test.authentication, {
    with(test.dataset(df), {
        test_that("Properly encoded UTF is sent and received", {
            s <- iconv("aided_follow_grid:ElCorteInglés", to="UTF-8")
            name(ds$v1) <- s
            expect_identical(name(ds$v1), s)
            expect_identical(name(refresh(ds)$v1), s)
            s2 <- "Budějovický Budvar"
            name(ds$v2) <- s2
            expect_identical(name(ds$v2), s2)
            expect_identical(name(refresh(ds)$v2), s2)
        })
    })
})
