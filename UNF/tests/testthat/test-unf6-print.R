context("UNFv6: Signature Printing")
test_that("Version is printed", {
    expect_equal(unf6(1)$formatted, "UNF6:tv3XYCv524AfmlFyVOhuZg==")
})
test_that("Digits are printed", {
    expect_equal(unf6(1, digits = 7)$formatted, "UNF6:tv3XYCv524AfmlFyVOhuZg==", label="digits == 7")
    expect_equal(unf6(1, digits = 1)$formatted, "UNF6:N1:tv3XYCv524AfmlFyVOhuZg==", label="digits < 7")
    expect_equal(unf6(1, digits = 10)$formatted, "UNF6:N10:tv3XYCv524AfmlFyVOhuZg==", label="digits > 7")
})
test_that("Characters are printed", {
    expect_equal(unf6(1, characters = 128)$formatted, "UNF6:tv3XYCv524AfmlFyVOhuZg==", label="characters == 128")
    expect_equal(unf6(1, characters = 1)$formatted, "UNF6:X1:BiVmCZpe+iuKPSR6rRl4nQ==", label="characters < 128")
})
test_that("Truncation is printed", {
    expect_equal(unf6(1, truncation = 128)$formatted, "UNF6:tv3XYCv524AfmlFyVOhuZg==", label="truncation == 128")
    expect_equal(unf6(1, truncation = 256)$formatted, "UNF6:H256:tv3XYCv524AfmlFyVOhuZo3W84VyoLXzTqyDU+3OHG8=", label="truncation > 128")
})

test_that("Multiple header features are printed", {
    expect_equal(unf6(1, digits = 1, characters = 1)$formatted, "UNF6:N1,X1:BiVmCZpe+iuKPSR6rRl4nQ==", label="digits < 7, characters < 128")
    expect_equal(unf6(1, digits = 1, truncation = 256)$formatted, "UNF6:N1,H256:tv3XYCv524AfmlFyVOhuZo3W84VyoLXzTqyDU+3OHG8=", label="digits > 7, truncation > 128")
    expect_equal(unf6(1, digits = 1, characters = 150, truncation=256)$formatted, "UNF6:N1,X150,H256:tv3XYCv524AfmlFyVOhuZo3W84VyoLXzTqyDU+3OHG8=", label="digits > 7, characters > 128, truncation > 128")
})
