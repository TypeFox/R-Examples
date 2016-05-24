#############################################################################
# test-style_default.R
# 
# Testing default style functions
#############################################################################

# style.default.character
#############################################################################

test_that("character takes character vectors",
{
    expect_identical(style.default.character("test"),
                     "test")
    expect_identical(style.default.character(c("test",
                                               "test")),
                     c("test", "test"))
})

test_that("character takes no empty input",
{
    expect_error(style.default.character())
})


# style.default.numeric
#############################################################################

test_that("numeric takes numeric vectors",
{
    expect_identical(style.default.numeric(12),
                     "12.00")
    expect_identical(style.default.numeric(c(12,
                                             13)),
                     c("12.00",
                       "13.00"))
})

test_that("numeric takes no empty input",
{
    expect_error(style.default.numeric())
})
