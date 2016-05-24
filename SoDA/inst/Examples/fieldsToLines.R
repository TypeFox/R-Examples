fieldsToLines <- function(input, nfield, sep = "",
                          outputSep = if(nchar(sep)) sep else " ",
                            addQuotes = FALSE, ...)
  {
    text <- scan(input, what = "", sep = sep, ...)
    if(addQuotes)
      text <- shQuote(text)
    text <- matrix(text, nrow = nfield)
    output <- apply(text, 2, paste, collapse = outputSep)
    textConnection(output)
  }
