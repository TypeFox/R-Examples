context("Get notice")

test_that("Retrieve A Notice", {
    (n <- ldnotice(1))
    expect_equal(class(n), "lumen_notice")
})

test_that("Print Notice Summary", {
    n <- ldnotice(1)
    summary(n)
})

test_that("Retrieve Entities", {
    (e <- ldentities(list(term = "joe")))
    expect_equal(class(e[[1]]), "lumen_entity")
})

test_that("Retrieve Topics", {
    e <- ldtopics()
    expect_equal(class(e[[1]]), "lumen_topic")
})

test_that("Search Works", {
    (s1 <- ldsearch(list(term = "youtube")))
    expect_equal(class(s1[[1]]), "lumen_notice")
    
    (s2 <- ldsearch(list(term = "youtube"), per_page = 2))
    expect_equal(length(s2), 2)
})
