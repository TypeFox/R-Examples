stresc <- function(strings)
{
  #
  # For some reason I don't remember, the objects "fromchar" and "tochar" were
  # originally arguments of the function to allow for more flexibillity but
  # are now hard-encoded within the function.
  #
  fromchar <- s2c("\\{}$^_%#&~[]|")
  tochar <- c("$\\backslash$", "\\{", "\\}", "\\$", "\\^{}", "\\_", "\\%", "\\#", "\\&", "\\~{}",
  "\\lbrack{}", "\\rbrack{}","\\textbar{}")
  #
  # Definition of the function to escape LaTeX character in one string:
  #
  f <- function(string){
    c2s( sapply(s2c(string), function(x) ifelse(x %in% fromchar, tochar[which(x == fromchar)], x)))
  }
  #
  # Now apply it to all elements of vector string:
  #
  sapply(strings, f, USE.NAMES = FALSE)
}
