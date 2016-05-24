
# Test ecitmatch() -----------------------------------------------------------

context("Testing 'ecitmatch()'")

if (getOption('reutils.test.remote')) {
  citstrings <- c("proc natl acad sci u s a|1991|88|3248|mann bj|Art1|",
                  "science|1987|235|182|palmenber ac|Art2|")
  x <- ecitmatch(citstrings)
  
  test_that("ecitmatch() content gets parsed correctly", {
    expect_is(x, 'ecitmatch')
    expect_equal(content(x), "proc natl acad sci u s a|1991|88|3248|mann bj|Art1|2014248\nscience|1987|235|182||Art2|3026048\n")
    expect_equal(content(x, "parsed"), c("2014248", "3026048"))
  })

}
