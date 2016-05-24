context("Test listOrganisations")
orgs <- listOrganisations()
test_that("List organisations", {
    expect_true(is.data.frame(orgs))
    expect_true(nrow(orgs)>200)
})
