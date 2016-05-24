library("papeR")
context("label functions")

data <- data.frame(a = 1:10, b = 10:1, c = rep(1:2, 5))

############################################################
## labels.default
############################################################

test_that("labels.default is dispatched correctly", {
    expect_equal(labels(matrix(1:10)),
                 list(as.character(1:10), as.character(1)))
})

############################################################
## labels / labels.data.frame
############################################################

library("foreign")
spss.data <- suppressWarnings(
    read.spss(system.file("SPSS/data.sav", package = "papeR"),
              to.data.frame = TRUE))
lbls <- c(x = "Predictor", y = "Outcome", Notes = "Additional Notes")

test_that("SPSS labels are correctly imported", {
    expect_equal(labels(spss.data), lbls)
    expect_equal(labels(spss.data, "x"), lbls[1])
    expect_false(is.ldf(spss.data))
})

test_that("abbreviate works as expected", {
    lbls <- c("This is a long label", "This is another long label",
              "This also")
    lbls_short <- c(a = "Thsisalngl", b = "Thsisantll", c = "This also")
    lbls_shorter <- c(a = "Thsisalnl", b = "Thsisanll", c = "Ta")
    labels(data) <- lbls
    names(lbls) <- c("a", "b", "c")
    expect_equal(labels(data), lbls)
    expect_equal(labels(data, abbreviate = TRUE, minlength = 10),
                 lbls_short)
    expect_equal(labels(data, abbreviate = TRUE, minlength = 2),
                 lbls_shorter)
})

############################################################
## setting labels
############################################################

test_that("labels can be set and reset", {
    labels(data) <- c("my_a", "my_b", "my_c")
    expect_equal(labels(data), c(a = "my_a", b = "my_b", c = "my_c"))
    expect_true(is.ldf(data))
    labels(data) <- NULL
    expect_equal(labels(data), c(a = "a", b = "b", c = "c"))
})


test_that("labels for subsets of the data can be set", {
    labels(data, which = c("a", "b")) <- c("x", "y")
    expect_equal(labels(data), c(a = "x", b = "y", c = "c"))
    expect_true(is.ldf(data))
    ## new variable
    data$z <- as.factor(rep(2:3, each = 5))
    expect_equal(labels(data), c(a = "x", b = "y", c = "c", z = "z"))
    labels(data, which = "z") <- "new_label"
    expect_equal(labels(data, "z"), c(z = "new_label"))
    ## subsets with [] operator
    labels(data)[1] <- "A"
    expect_equal(labels(data, "a"), c(a = "A"))
})

test_that("labels for subsets of a labeled data frame can be set", {
    data <- as.ldf(data)
    labels(data, which = c("a", "b")) <- c("x", "y")
    expect_equal(labels(data), c(a = "x", b = "y", c = "c"))
})

############################################################
## CLEAN_LABELS
############################################################

test_that("label cleaning works", {
    ## drop variable [note that this also drops all non-data.frame attributes]
    expect_equal(labels(spss.data[-1]), c(y = "y", Notes = "Notes"))
    ## add variable
    spss.data$z <- 4:1
    expect_equal(labels(spss.data), c(lbls, z = "z"))
    ## reorder data [note that this also drops all non-data.frame attributes]
    expect_equal(labels(spss.data[, c(4, 2, 1, 3)]),
                 c(z = "z", y = "y", x = "x", Notes = "Notes"))
    ## rename variable
    names(spss.data)[3] <- "comments"
    expect_equal(labels(spss.data),
                 c(x = "Predictor", y = "Outcome", comments = "comments", z = "z"))
})

############################################################
## as.ldf / convert.labels / is.ldf
############################################################

test_that("conversion of labels works (1)", {
    spss.data <- convert.labels(spss.data)
    expect_true(is.ldf(spss.data))
    expect_equivalent(labels(spss.data$x), labels(spss.data, "x"))
    expect_equivalent(labels(spss.data$x, abbreviate = TRUE, minlength = 5),
                      labels(spss.data, "x", abbreviate = TRUE, minlength = 5))
    expect_equal(labels(spss.data, "x"), lbls[1])
})

test_that("conversion of labels works (2)", {
    spss.data <- as.ldf(spss.data)
    expect_true(is.ldf(spss.data))
    expect_equivalent(labels(spss.data$x), labels(spss.data, "x"))
    expect_equivalent(labels(spss.data$x, abbreviate = TRUE, minlength = 5),
                      labels(spss.data, "x", abbreviate = TRUE, minlength = 5))
    expect_equal(labels(spss.data, "x"), lbls[1])
})
