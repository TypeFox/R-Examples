library(postGIStools)
context("hstore operator")

# Create example data frame with empty hstore
contacts <- data.frame(name = c("Anne", "Bert", "Chris"))
contacts$phone <- new_hstore(3)

test_that("%->% works on empty hstore", {
    expect_equal(contacts$phone %->% "home", rep(NA, 3))
})

# Make assignments to whole hstore and subset, including assignment to NULL
contacts$phone %->% "home" <- c("555-123-4567", "555-923-9134", "555-276-1123")
contacts$phone[2] %->% "home" <- NULL
contacts$phone[2:3] %->% "cell" <- c("555-889-9134", "555-852-0137")
rownames(contacts) <- c("A", "B", "C")

test_that("%->% works with various subsetting", {
    expect_equal(contacts$phone %->% "home", c("555-123-4567", NA, "555-276-1123"))
    expect_equal(contacts[-1, 2] %->% "cell", c("555-889-9134", "555-852-0137"))
    expect_equal(contacts$phone[contacts$name == "Anne"] %->% "cell", NA)
    expect_equal(contacts["C", "phone"] %->% "cell", "555-852-0137")
})


test_that("assignment works between hstores", {
    new_contacts <- new_hstore(2)
    new_contacts %->% "home" <- contacts$phone[2:3] %->% "home"
    expect_equal(new_contacts %->% "home", contacts$phone[2:3] %->% "home")
})


test_that("new_hstore and `%->%` fail on bad inputs", {
    expect_error(new_hstore("1"))
    expect_error(new_hstore(1:2))
    expect_error(contacts$name %->% "name")
    expect_error(contacts$phone[[1]] %->% "home")
    expect_error(contacts$phone %->% 1)
    expect_error(contacts$phone %->% "home" <- "555-123-4567")
    expect_error(contacts$phone %->% c("home", "cell"))
})

