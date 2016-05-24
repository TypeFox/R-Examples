context("Converting sequence to n-grams")

test_that("Extract ngrams for different distances",{
  seqs <- structure(c(4L, 2L, 1L, 1L, 4L, 3L, 1L, 3L, 1L, 1L), .Dim = c(2L, 
                                                                        5L))
  expect_equal(seq2ngrams(seqs, 3, 1L:4), 
               structure(c("4.1.4_0.0", "2.1.3_0.0", "1.4.1_0.0", "1.3.3_0.0", 
                           "4.1.1_0.0", "3.3.1_0.0"), .Dim = 2:3))
  
  #the same with distance
  expect_equal(seq2ngrams(seqs, 3, 1L:4, d = c(0, 1)), 
               structure(c("4.1.1_0.1", "2.1.3_0.1", "1.4.1_0.1", "1.3.1_0.1"), 
                         .Dim = c(2L, 2L)))
  
  #only one n-gram possible
  expect_equal(seq2ngrams(seqs, 3, 1L:4, d = 1), 
               structure(c("4.4.1_1.1", "2.3.1_1.1"), .Dim = c(2L, 1L)))
  
  #are all extracted n-gram names correct? 
  expect_true(all(sapply(seq2ngrams(sample(1L:4, 20, replace = TRUE), 3, 1L:4), 
                         is_ngram)))
  
  #are all extracted n-gram names correct? same as above, but for letters
  expect_true(all(sapply(seq2ngrams(sample(c("a", "c", "g", "t"), 20, replace = TRUE), 3, c("a", "c", "g", "t")), 
                         is_ngram)))
  
})
