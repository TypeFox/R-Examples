# tests if the argument is a start-of-code-chunk marker
startOfCodeChunk <- function(x)
{
   if (!is.null(x))
   {
      s <- strsplit(trim(x), '')[[1]]
      n <- length(s)
      n >= 5 && all(s[c(1, 2)] == '<' & s[c(n-2, n-1)] == '>') && s[n] == '='
   } else {
      FALSE
   }
}

# tests if the argument is an end-of-code-chunk marker
endOfCodeChunk <- function(x)
{
   ! is.null(x) && trim(x) == '@'
}

# tests if the argument contains \Sexpr{...} or \SweaveOpts{...}
startOfCodeExpr <- function(x)
{
   if (is.null(x))
      x <- ''
   regexpr('\\\\(Sexpr|SweaveOpts)\\{', x)
}

endOfCodeExpr <- function(x)
{
   if (is.null(x))
      x <- ''
   regexpr('\\}', x)
}
