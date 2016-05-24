quotewords <- function(x, delim = "\\s+")  {
    .PerlPackage("Text::ParseWords")
    .Perl("quotewords", .args = c(list(delim, 1), x))
}
