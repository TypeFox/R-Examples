#############################################################################
# test-out_all.R
# 
# Testing all output functions together
#############################################################################

# out.default.specialchar
# out.html.specialchar
# out.latex.specialchar
#############################################################################

test_that("specialchar takes character vectors",
{
    # arguments
    t1 <- "test"
    t2 <- c("test1", "<", "test2")

    # character vector (length = 1)
    expect_identical(out.default.specialchar(t1),
                     "test")
    expect_identical(out.html.specialchar(t1), 
                     "test") 
    expect_identical(out.latex.specialchar(t1), 
                     "test") 

    #  character vector (length > 1)
    expect_identical(out.default.specialchar(t2),
                     c("test1", 
                       "<",
                       "test2"))
    expect_identical(out.html.specialchar(t2), 
                     c("test1", 
                       "&lt;",
                       "test2"))
    expect_identical(out.latex.specialchar(t2), 
                     c("test1", 
                       "\\ifmmode<\\else\\textless\\fi",
                       "test2"))
})

test_that("specialchar takes no empty input",
{
    expect_error(out.default.specialchar())
    expect_error(out.html.specialchar())
    expect_error(out.latex.specialchar())
})

test_that("specialchar takes no multiple arguments",
{
    # arguments
    t1  <- "test1"
    t2  <- "test2"

    expect_error(out.default.specialchar(t1, t2))
    expect_error(out.html.specialchar(t1, t2))
    expect_error(out.latex.specialchar(t1, t2))
})


# out.default.math
# out.html.math
# out.latex.math
#############################################################################

test_that("math takes character vectors",
{
    # arguments
    t1 <- "test"
    t2 <- c("test", "test")
    tmmode <- TRUE

    # character vector (length = 1)
    expect_identical(out.default.math(t1, 
                                      mmode = tmmode), 
                     "test")
    expect_identical(out.html.math(t1, 
                                    mmode = tmmode), 
                     "<MATH>test</MATH>")
    expect_identical(out.latex.math(t1, 
                                    mmode = tmmode), 
                     "\\ensuremath{test}")

    # character vector (length > 1)
    expect_identical(out.default.math(t2, 
                                      mmode = tmmode), 
                     c("test",
                       "test"))
    expect_identical(out.html.math(t2, 
                                    mmode = tmmode), 
                     c("<MATH>test</MATH>",
                       "<MATH>test</MATH>"))
    expect_identical(out.latex.math(t2, 
                                    mmode = tmmode), 
                     c("\\ensuremath{test}",
                       "\\ensuremath{test}"))

    # arguments
    tmmode <- FALSE

    # character vector (length = 1)
    expect_identical(out.default.math(t1, 
                                      mmode = tmmode), 
                     "test")
    expect_identical(out.html.math(t1, 
                                    mmode = tmmode), 
                     "test")
    expect_identical(out.latex.math(t1, 
                                    mmode = tmmode), 
                     "test")

    # character vector (length > 1)
    expect_identical(out.default.math(t2, 
                                      mmode = tmmode), 
                     c("test",
                       "test"))
    expect_identical(out.html.math(t2, 
                                    mmode = tmmode), 
                     c("test",
                       "test"))
    expect_identical(out.latex.math(t2, 
                                    mmode = tmmode), 
                     c("test",
                       "test"))
})

test_that("math takes multiple arguments",
{
    # arguments
    t1 <- "test1"
    t2 <- "test2"
    t3 <- c("test1", "test2")
    t4 <- c("test3", "test4")
    tmmode <- TRUE

    # character vector (length = 1)
    expect_identical(out.default.math(t1, t2,
                                      mmode = tmmode), 
                     "test1test2")
    expect_identical(out.html.math(t1, t2,
                                    mmode = tmmode), 
                     "<MATH>test1test2</MATH>")
    expect_identical(out.latex.math(t1, t2,
                                    mmode = tmmode), 
                     "\\ensuremath{test1test2}")

    # character vector (length > 1)
    expect_identical(out.default.math(t3, t4,
                                      mmode = tmmode), 
                     c("test1test3",
                       "test2test4"))
    expect_identical(out.html.math(t3, t4,
                                    mmode = tmmode), 
                     c("<MATH>test1test3</MATH>",
                       "<MATH>test2test4</MATH>"))
    expect_identical(out.latex.math(t3, t4,
                                    mmode = tmmode), 
                     c("\\ensuremath{test1test3}",
                       "\\ensuremath{test2test4}"))

    # arguments
    tmmode <- FALSE

    # character vector (length = 1)
    expect_identical(out.default.math(t1, t2,
                                      mmode = tmmode), 
                     "test1test2")
    expect_identical(out.html.math(t1, t2,
                                    mmode = tmmode), 
                     "test1test2")
    expect_identical(out.latex.math(t1, t2,
                                    mmode = tmmode), 
                     "test1test2")

    # character vector (length > 1)
    expect_identical(out.default.math(t3, t4,
                                      mmode = tmmode), 
                     c("test1test3",
                       "test2test4"))
    expect_identical(out.html.math(t3, t4,
                                    mmode = tmmode), 
                     c("test1test3",
                       "test2test4"))
    expect_identical(out.latex.math(t3, t4,
                                    mmode = tmmode), 
                     c("test1test3",
                       "test2test4"))
})

test_that("math needs mmode argument",
{
    expect_error(out.default.math())
    expect_error(out.html.math())
    expect_error(out.latex.math())
})

test_that("math takes empty output",
{
    # arguments
    tmmode <- TRUE

    expect_identical(out.default.math(mmode = tmmode),
                     character(0))
    expect_identical(out.html.math(mmode = tmmode),
                     "<MATH></MATH>")
    expect_identical(out.latex.math(mmode = tmmode),
                     "\\ensuremath{}")

    # arguments
    tmmode <- FALSE

    expect_identical(out.default.math(mmode = tmmode),
                     character(0))
    expect_identical(out.html.math(mmode = tmmode),
                     character(0))
    expect_identical(out.latex.math(mmode = tmmode),
                     character(0))
})


# out.default.subscript
# out.html.subscript
# out.latex.subscript
#############################################################################

test_that("subscript takes character vectors",
{
    # arguments
    t1 <- "2"
    t2 <- c("2", "3")

    # character vector (length = 1)
    expect_identical(out.default.subscript(t1), 
                     "2") 
    expect_identical(out.html.subscript(t1), 
                     "<sub>2</sub>") 
    expect_identical(out.latex.subscript(t1), 
                     "\\ifmmode_{2}\\else\\textsubscript{2}\\fi") 

    # character vector (length > 1)
    expect_identical(out.default.subscript(t2), 
                     c("2",
                       "3")) 
    expect_identical(out.html.subscript(t2), 
                     c("<sub>2</sub>",
                       "<sub>3</sub>"))
    expect_identical(out.latex.subscript(t2), 
                     c("\\ifmmode_{2}\\else\\textsubscript{2}\\fi",
                       "\\ifmmode_{3}\\else\\textsubscript{3}\\fi"))
})

test_that("specialchar takes no multiple arguments",
{
    # arguments
    t1 <- "test"

    expect_error(out.default.subscript(t1, t1))
    expect_error(out.html.subscript(t1, t1))
    expect_error(out.latex.subscript(t1, t1))
})

test_that("specialchar takes no empty input",
{
    expect_error(out.default.subscript())
    expect_error(out.html.subscript())
    expect_error(out.latex.subscript())
})


# out.default.superscript
# out.html.superscript
# out.latex.superscript
#############################################################################

test_that("superscript takes character vectors",
{
    # arguments
    t1 <- "2"
    t2 <- c("2", "3")

    # character vector (length = 1)
    expect_identical(out.default.superscript(t1), 
                     "2") 
    expect_identical(out.html.superscript(t1), 
                     "<sup>2</sup>") 
    expect_identical(out.latex.superscript(t1), 
                     "\\ifmmode^{2}\\else\\textsuperscript{2}\\fi") 

    # character vector (length > 1)
    expect_identical(out.default.superscript(t2), 
                     c("2",
                       "3")) 
    expect_identical(out.html.superscript(t2), 
                     c("<sup>2</sup>",
                       "<sup>3</sup>"))
    expect_identical(out.latex.superscript(t2), 
                     c("\\ifmmode^{2}\\else\\textsuperscript{2}\\fi",
                       "\\ifmmode^{3}\\else\\textsuperscript{3}\\fi"))
})

test_that("specialchar takes no multiple arguments",
{
    # arguments
    t1 <- "test"

    expect_error(out.default.superscript(t1, t1))
    expect_error(out.html.superscript(t1, t1))
    expect_error(out.latex.superscript(t1, t1))
})

test_that("specialchar takes no empty input",
{
    expect_error(out.default.superscript())
    expect_error(out.html.superscript())
    expect_error(out.latex.superscript())
})

# out.default.bracket
#############################################################################

test_that("bracket takes character vectors",
{
    # arguments
    t1 <- "t(12) = 5"
    t2 <- "t(15) = 8"
    t3 <- paste(t1, t2)
    t4 <- "( t[12] = 5 t[15] = 8 )"
    t5 <- c(t1, t2)

    # character vector (length = 1)
    expect_identical(out.default.bracket(t1,
                                         brackets = c("(", ")", "[", "]")), 
                     "(t[12] = 5)") 
    expect_identical(out.default.bracket(t3,
                                         brackets = c("(", ")", "[", "]")), 
                     "(t[12] = 5 t[15] = 8)") 
    expect_identical(out.default.bracket(t4,
                                         brackets = c("(", ")", "[", "]")), 
                     "([ t(12) = 5 t(15) = 8 ])") 

    # character vector (length > 1)
    expect_identical(out.default.bracket(t5,
                                         brackets = c("(", ")", "[", "]")), 
                     c("(t[12] = 5)", "(t[15] = 8)"))
})

test_that("bracket takes no multiple arguments",
{
    # arguments
    t1 <- "test"

    expect_error(out.default.bracket(t1, t1, 
                                     brackets = c("(", ")", "[", "]")))
})

test_that("bracket takes no empty input",
{
    expect_error(out.default.bracket(brackets = c("(", ")", "[", "]")))
})

# out.default.above
# out.latex.above
#############################################################################

test_that("above takes character vectors",
{
    # arguments
    t1 <- "x"
    t2 <- "y"

    # character vector (length = 1)
    expect_identical(out.default.above(t1, t2),
                     "xy") 
    expect_identical(out.latex.above(t1, t2),
                     "\\ifmmode\\overset{y}{x}\\else x\\textsuperscript{y}\\fi")
})

## out.default.below
## out.latex.below
##############################################################################

test_that("below takes character vectors",
{
    # arguments
    t1 <- "x"
    t2 <- "y"

    # character vector (length = 1)
    expect_identical(out.default.below(t1, t2),
                     "xy") 
    expect_identical(out.latex.below(t1, t2),
                     "\\ifmmode\\underset{y}{x}\\else x\\textsubscript{y}\\fi")
})
