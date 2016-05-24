context("UNF: Version dispatch")

test_that("Version dispatch", {
    expect_equivalent(attr(unf(1:3, version=3), "version"), 3, label="Version 3")
    expect_equivalent(attr(unf(1:3, version=4), "version"), 4, label="Version 4")
    expect_equivalent(attr(unf(1:3, version=4.1), "version"), 4.1, label="Version 4.1")
    expect_equivalent(attr(unf(1:3, version=5), "version"), 5, label="Version 5")
    expect_equivalent(attr(unf(1:3, version=6), "version"), 6, label="Version 6")
    expect_error(attr(unf(1:3, version=2), "version"), label="Unrecognized numeric version")
    expect_error(attr(unf(1:3, version="a"), "version"), label="Unrecognized character version")
})
