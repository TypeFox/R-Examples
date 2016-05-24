naturalorder <- function(text, decreasing=FALSE, na.last=TRUE) {  # different with base::order in order or arguments
   if (!is.logical(decreasing) || length(decreasing) != 1) {
      decreasing <- as.logical(decreasing)[1]
   }
   if (is.na(decreasing)) {
      stop("'decreasing' must be either TRUE or FALSE")
   }
   if (!is.logical(na.last) || length(na.last) != 1) {
      na.last <- as.logical(na.last)[1]
   }
   if (!is.character(text)) {
      text <- as.character(text)
   }
   if (length(text) == 0L) {
      return(integer(0L))
   }
   
   sign <- (-1L) ^ decreasing
   removingNA <- is.na(na.last)  # used at last to remove NAs
   ## na.last | (is.na(na.last) || na.last)
   ## --------+----------------------------
   ## NA      | TRUE
   ## TRUE    | TRUE
   ## FALSE   | FALSE
   na.last <- xor(is.na(na.last) || na.last, decreasing)
   
   ## If strsplit is applied to an empty character, an empty character vector is returned.
   ## Therefore, if all elements in 'text' are empty, 'maxLength' will be 0.
   ## Otherwise, when there is at least one ordinal value or NA in 'text', 'maxLength' will be greater than 0.
   tokenList <- strsplit(text, "(?<=\\d)(?=\\D)|(?<=\\D)(?=\\d)", perl=TRUE)
   maxLength <- max(sapply(tokenList, length))
   if (maxLength == 0L) {  # all elements are empty ("").
      return(seq_along(text))
   }
   tokenList <- lapply(tokenList, function(tokens) c(tokens, rep("", maxLength - length(tokens))))
   tokenList <- Reduce(rbind, tokenList, matrix(, 0, maxLength))
   tokenList <- as.data.frame(tokenList, stringsAsFactors=FALSE)
   
   ranks <- lapply(tokenList, function(tokens) {
      isInteger <- grepl("^\\d+$", tokens, useBytes=TRUE)
      ## stability for zero-padding equivalent values
      ## (e.g. 1, 01, 001, ...)
      integers <- ifelse(isInteger, tokens, "")
      zeroPaddingRank <- sign * rank(integers, na.last=TRUE)
      ## sort as integer values
      ## string values sholud be ranked identically
      integers <- rep(-1, length(text))
      integers[isInteger] <- as.integer(tokens[isInteger])
      integerRank <- sign * rank(integers, na.last=TRUE)
      ## sort as string values
      ## integer values sholud be ranked identically
      strings <- ifelse(isInteger, "0", tokens)
      stringRank <- sign * rank(strings, na.last=na.last)
      list(stringRank, integerRank, zeroPaddingRank)
   })
   ranks <- unlist(ranks, recursive=FALSE)
   orderFunction <- sprintf("order(%s)", paste(names(ranks), collapse=","))
   result <- with(ranks, eval(parse(text=orderFunction)))
   if (removingNA) {
      result <- result[!(result %in% which(is.na(text)))]
   }
   result
}

naturalsort <- function(text, decreasing=FALSE, na.last=NA) {
   text[naturalorder(text, decreasing=decreasing, na.last=na.last)]
}
