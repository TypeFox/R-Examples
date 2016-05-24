capitalize <- function(strings){
  ## capitalize characters that begin strings, or that follow a character that
  ## is not a digit, underscore, or letter
  upstrings <- toupper(strings)
  downstrings <- tolower(strings)
  marks <- gsub("^[[:alpha:]]", "X", gsub("[^[:alnum:]_][[:alnum:]_]", " X", downstrings))
  for(i in 1:length(strings)){
    n <- nchar(strings[i])
    pos <- (1:n)[substring(marks[i], 1:n, 1:n) == "X"]
    chars <- substring(strings[i], 1:n, 1:n)
    chars[pos] <- substring(upstrings[i], 1:n, 1:n)[pos]
    strings[i] <- paste(chars, collapse = "")
  }
  strings
}

