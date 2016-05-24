context("UTF-8 Encoding")

## Meta test of encoding of the test files
test_that("All test- files are ASCII (for CHECK)", {
    testfiles <- dir(pattern="^test.*R$")
    for (i in testfiles) {
        expect_warning(scan(i, what=character(), fileEncoding="ascii", quiet=TRUE),
            NA, info=i)
    }
})

if (run.integration.tests) {
    ## Move the actual tests to a different file so that non-Unicode
    ## environments won't fail to parse (even if not running these tests)
    source("utftesting.R", local=TRUE)
}
