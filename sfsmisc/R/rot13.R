##' Generalized  rot13   --> ../man/rot13.Rd
rotn <- function (ch, n = 13)
{
    ch <- as.character(ch) # or error
    stopifnot(0 <= n, n <= 26)
    i <- c(if(n < 26) (n+1):26, seq_len(n))
    chartr(old = paste(c(letters,   LETTERS   ), collapse=""),
	   new = paste(c(letters[i],LETTERS[i]), collapse=""),
	   x = ch)
}
