context("Checking run_split")

test_that("run_split works on character strings",{

    expect_equal(
        run_split(c("122333444455555666666", NA, "abbcccddddeeeeeffffff")),
        list(c("1", "22", "333", "4444", "55555", "666666"), NA_character_, 
            c("a", "bb", "ccc", "dddd", "eeeee", "ffffff"))
    )

})

