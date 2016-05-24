context("check_list_class")

test_that("check_list_class behaves as expected", {
    expect_true(check_list_class(list(numeric(0), numeric(0), numeric(0)), "numeric"))
    expect_true(check_list_class(list(character(0), character(0), character(0)), "character"))

    ### Use named lists
    expect_true(check_list_class(list(x = numeric(0), y = numeric(0), z = numeric(0)), "numeric"))
    expect_true(check_list_class(list(x = character(0), y = character(0), z = character(0)), "character"))

    ### Empty list
    expect_warning(check_list_class(list(), "character"), "no elements")
    expect_true(suppressWarnings(check_list_class(list(), "character")))
})
