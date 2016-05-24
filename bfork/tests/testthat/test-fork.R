dest <- "./test.csv"
fn <- function() {
    write.csv(writedf, dest)
}

writedf <- data.frame(a=c(1,2,3,4), b=1,2,3,4)
cpid <- fork(fn)
rpid <- waitpid(cpid)

test_that("forked pid equals returned pid", {
          expect_equal(cpid, rpid)
})

test_that("forked process did create a file", {
          expect_true(file.exists(dest))
})

readdf <- read.csv(dest, header=TRUE)

test_that("written data fram is the same as the read data frame", {
          expect_equal(readdf$a, writedf$a)
          expect_equal(readdf$b, writedf$b)
})

file.remove(dest)

test_that("check to ensure cleanup", {
          expect_false(file.exists(dest))
})
