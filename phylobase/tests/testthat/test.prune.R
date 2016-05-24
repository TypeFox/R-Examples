#
# --- Test prune.R ---
#

data(geospiza)
gtree <- extractTree(geospiza)

context("prune")

test_that("prune works on phylo4 objects", {
    # function(phy, tip, trim.internal = TRUE, subtree = FALSE, ...)
    expect_equal(gtree, prune(gtree, character(0)))
})

test_that("prune works on phylo4d objects", {
    # function(phy, tip, trim.internal = TRUE, subtree = FALSE, ...)
    expect_equal(geospiza, prune(geospiza, character(0)))
})
