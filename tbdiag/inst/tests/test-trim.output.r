
################################################################################
# trim.output defaults to terse when an invalid verbosity argument is passed
test_that("trim.output defaults to terse and gives a warning when an invalid verbosity argument is passed", {

expect_that(trim.output("Positive dood!"), equals("Positive"))
expect_that(trim.output("Positive dood!", "loquacious"), equals("Positive"))
expect_that(trim.output("Positive dood!", "1"), equals("Positive"))
expect_that(trim.output("Positive dood!", NA), equals("Positive"))

# trim.output warns when an invalid verbosity argument is passed
expect_that(trim.output("Positive dood!", "loquacious"), gives_warning())
expect_that(trim.output("Positive dood!", "1"), gives_warning())
expect_that(trim.output("Positive dood!", NA), gives_warning())

})


################################################################################
# trim.output returns appropriate result permutations
test_that("trim.output returns appropriate result permutations", {

expect_that(trim.output("Positive dood!"), equals("Positive"))
expect_that(trim.output("Positive dood!", "onechar"), equals("P"))
expect_that(trim.output("Positive dood!", "terse"), equals("Positive"))
expect_that(trim.output("Positive dood!", "verbose"), 
            equals("Positive dood!"))
expect_that(trim.output("Positive dood!", "loquacious"), equals("Positive"))

})


