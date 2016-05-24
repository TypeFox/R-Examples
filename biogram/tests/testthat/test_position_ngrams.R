context("Converting n-grams to positions")

test_that("Extract ngrams for different distances",{
  seqs <- structure(c(4L, 2L, 1L, 1L, 4L, 3L, 1L, 3L, 1L, 1L), .Dim = c(2L, 
                                                                        5L))
  expect_equal(position_ngrams("2_1.1.2_0.1"), 
               structure(list(`2` = structure(1L, .Label = c("1_0", "2_0"), class = "factor"), 
                              `3` = structure(1L, .Label = c("1_0", "2_0"), class = "factor"), 
                              `5` = structure(2L, .Label = c("1_0", "2_0"), class = "factor")), 
                         .Names = c("2", "3", "5")))
  
  #multiple n-grams
  expect_equal(position_ngrams(c("2_1.1.2_0.1", "3_1.1.2_0.0", "3_2.2.2_0.0")),
               structure(list(`2` = structure(1L, .Label = c("1_0", "2_0"), class = "factor"), 
                              `3` = structure(c(1L, 1L, 2L), .Label = c("1_0", "2_0"), class = "factor"), 
                              `4` = structure(1:2, .Label = c("1_0", "2_0"), class = "factor"), 
                              `5` = structure(c(2L, 2L, 2L), .Label = c("1_0", "2_0"), class = "factor")), 
                         .Names = c("2", "3", "4", "5")))
  
  #no overlap in positions
  expect_equal(position_ngrams(c("2_1.1.2_0.1", "7_1.1.2_0.0", "10_2.2.2_0.0")),
               structure(list(`2` = structure(1L, .Label = c("1_0", "2_0"), class = "factor"), 
                              `3` = structure(1L, .Label = c("1_0", "2_0"), class = "factor"), 
                              `5` = structure(2L, .Label = c("1_0", "2_0"), class = "factor"), 
                              `7` = structure(1L, .Label = c("1_0", "2_0"), class = "factor"), 
                              `8` = structure(1L, .Label = c("1_0", "2_0"), class = "factor"), 
                              `9` = structure(2L, .Label = c("1_0", "2_0"), class = "factor"), 
                              `10` = structure(2L, .Label = c("1_0", "2_0"), class = "factor"), 
                              `11` = structure(2L, .Label = c("1_0", "2_0"), class = "factor"), 
                              `12` = structure(2L, .Label = c("1_0", "2_0"), class = "factor")), 
                         .Names = c("2", "3", "5", "7", "8", "9", "10", "11", "12")))
})
