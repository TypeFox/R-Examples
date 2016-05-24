context("list_classes")

test_that("list_classes behaves as expected", {
    expect_equal(list_classes(list(numeric(0), integer(0), character(0))),
                 c("numeric", "integer", "character"))

    ### Use a named list
    expect_that(list_classes(list(x = numeric(0), y = integer(0), z = character(0))),
                is_equivalent_to(c("numeric", "integer", "character"))) # Need to use is_equivalent_to here to ignore attributes

    ### Empty list
    expect_equal(list_classes(list()), c())
})
