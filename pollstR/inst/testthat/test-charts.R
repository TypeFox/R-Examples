context("charts")

chartsu <- pollstr_charts(convert=FALSE)

test_that("charts is in expected format if not converted", {
    expect_is(chartsu, "list")
    expect_true(length(chartsu) > 0)
})

charts <- pollstr_charts()

test_that("charts is in expected format", {
    expect_is(charts, "pollstr_charts")
    expect_equal(length(charts), 2)
    expect_equal(names(charts), c("charts", "estimates"))
    expect_equal(unname(sapply(charts, class)), rep("data.frame", 2))
})

test_that("query returns data", {
    expect_true(nrow(charts[["charts"]]) > 0)
    expect_true(nrow(charts[["estimates"]]) > 0)
})

