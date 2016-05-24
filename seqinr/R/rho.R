rho <- function (sequence, wordsize = 2, alphabet = s2c("acgt"))
{
  wordcount <- count(sequence, wordsize, freq = FALSE, alphabet = alphabet)
  uni <- count(sequence, 1, freq = TRUE, alphabet = alphabet)

  expected_wordfreq <- function (wordsize, uni) 
  {
    if (wordsize == 1) 
        return(uni)
    else kronecker(uni, expected_wordfreq(wordsize - 1, uni))
  }

  expected_wordcount <- sum(wordcount)*expected_wordfreq(wordsize, uni)
   
  return(wordcount/expected_wordcount)
}
