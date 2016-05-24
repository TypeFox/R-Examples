context("Testing compareAttributes")




context("Test compareAList")

# some lists
l1 <- list(b=1:5, c=5:1)
l2 <- list(a=2, b=1:5, c=1:5, d=3)

          

test_that("compareAlist returns NA when intersection of names sets is empty", {
          l1 <- list(a=2, b=1:5, c=1:5, d=3)
          l2 <- list(B=1:5, C=5:1)
          expect_identical(compareAlist(l1, l2), as.character(NA))
} )
