context("Checking print.regex")

test_that("print.regex prints a character vector",{

    m <- "foo" %:)% "bar"
    expect_equivalent(capture.output(regexr:::print.subcom(m)), "[1] \"foo\"")
    
})



