########################################################################
# ---------------------------------------------------------------------#
# Tests for objects from the imbalance class and associated functions  #
# ---------------------------------------------------------------------#
########################################################################
context("Imbalance")

test_that("imbalance returns valid object", {
	
    type <- sample(c("imb", "absImb", "loss", "maxImb"), 1)
	  expect_is(imbal(type = type), "imbal")
	  expect_is(imbal(type = type), "issue")
	
	  type <- sample(c(2, "wrongtype"), 1)
	  expect_error(imbal(type = type))
  }
)


test_that("imbalance is working correctly", {
    CR <- getAllSeq(crPar(4))
    is1 <- imbal("imb")
    is2 <- imbal("absImb")
    is3 <- imbal("loss")
    is4 <- imbal("maxImb")
    A <- assess(CR, is1, is2, is3, is4)@D  
    
    expect_equal(A$imb, c(-4, -2, -2, 0, -2, 0, 0, 2, -2, 0, 0, 2, 0, 2, 2, 4))
    expect_equal(A$absImb, c(4, 2, 2, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0, 2, 2, 4))
    expect_equal(A$loss, c(4, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 4))
    expect_equal(A$maxImb, c(4, 2, 2, 2, 2, 1, 1, 3, 3, 1, 1, 2, 2, 2, 2, 4))
    
  }
)