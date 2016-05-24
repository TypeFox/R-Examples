library(unittest, quietly = TRUE)

expect_equal <- function (expr, expected) {
    actual <- capture.output(expr)
    if (!identical(all.equal(actual, expected), TRUE)) {
        stop("ok_group output didn't match: ", actual)
    }
}

# Expression printed
expect_equal(ok_group("camels", TRUE), c(
    "# camels"))

# Can have multiple lines in a vector
expect_equal(ok_group(c("camels", "ostriches"), TRUE), c(
    "# camels",
    "# ostriches"))

# Can divide lines with newlines too
expect_equal(ok_group(c("camels\nbadgers\r\nhoney badgers", "ostriches"), FALSE), c(
    "# camels",
    "# badgers",
    "# honey badgers",
    "# ostriches"))

# Expression evaluated after printing section message
expect_equal(ok_group("camels", print("moo")), c(
    "# camels",
    '[1] "moo"'))

# Return NULL
expect_equal({
    if (!is.null(ok_group("camels", 6))) stop("Didn't return NULL")
}, c("# camels"))
