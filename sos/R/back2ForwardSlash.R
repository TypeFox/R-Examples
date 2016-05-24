back2ForwardSlash <- function(nmax=1, what=character(), sep='\n', ...){
  x2 <- scan(what=what, nmax=nmax, sep=sep, ...)
  x. <- gsub('\\', '/', x2, fixed = TRUE)
}
