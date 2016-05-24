require(tm)
test_that("Integration with tm", {
              nirvana <- c("hello hello hello how low", "hello hello hello how low",
                           "hello hello hello how low", "hello hello hello",
                           "with the lights out", "it's less dangerous",
                           "here we are now", "entertain us",
                           "i feel stupid", "and contagious", "here we are now", "entertain us",
                           "a mulatto", "an albino", "a mosquito", "my libido", "yeah", "hey yay")
              expect_true("TermDocumentMatrix" %in% class(TermDocumentMatrix(Corpus(VectorSource(nirvana)), control = list(tokenize = function(x) ngramrr(x, ngmax = 2)))))
              expect_true("TermDocumentMatrix" %in% class(TermDocumentMatrix(Corpus(VectorSource(nirvana)), control = list(tokenize = function(x) ngramrr(x, ngmax = 2, char = TRUE)))))
          })
         
