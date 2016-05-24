rot13 <- function(string){
  if(!is.character(string)) stop("character string expected")
  old <- c2s(c(letters, LETTERS))
  new <- c2s(c(letters[14:26], letters[1:13], LETTERS[14:26], LETTERS[1:13]))
  chartr(old = old, new = new, x = string)
}
