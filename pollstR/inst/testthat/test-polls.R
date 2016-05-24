context("polls")

pollsu <- pollstr_polls(max_pages = 1, convert = FALSE)

test_that("polls is in expected format if not converted", {
    expect_is(pollsu, "list")
    expect_true(length(pollsu) > 0)
})

polls <- pollstr_polls(max_pages = 1)

test_that("polls data is in the correct format", {
    expect_is(polls, "pollstr_polls")
    expect_equal(length(polls), 4)
    expect_equal(names(polls), c("polls", "questions", "survey_houses", "sponsors"))
    expect_equal(unname(sapply(polls, class)), rep("data.frame", 4))
})

test_that("query returns data", {
    expect_true(nrow(polls[["polls"]]) > 0)
    expect_true(nrow(polls[["questions"]]) > 0)
})

