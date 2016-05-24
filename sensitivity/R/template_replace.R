template.replace <- function(text, replacement, eval = FALSE,
                             key.pattern = NULL, code.pattern = NULL) {
  if (is.null(key.pattern)) {
    key.pattern <- "\\$\\(KEY\\)"
  }
  if (is.null(code.pattern)) {
    code.pattern = "@\\{CODE\\}"
  }
  code.pattern <- sub("CODE", ".+?", code.pattern)
  
  # loop on the template text lines
      
  for (i in 1 : length(text)) {
      
    # replacement of the keys

    for (keyname in names(replacement)) {
      text[i] <- gsub(sub("KEY", keyname, key.pattern, perl = TRUE),
                      paste(replacement[[keyname]]), text[i], perl = TRUE)
    }

    if (eval) {
      
      # code evaluation

      reg <- regexpr(code.pattern, text[i], perl = TRUE)
      while (reg > -1) {
        matched.first <- as.numeric(reg)
        matched.last <- matched.first + attr(reg, "match.length") - 1
        matched.text <- substr(text[i], matched.first + 2, matched.last - 1)
        
        val.matched.text <- eval(parse(text = matched.text))

        line.begin <- substr(text[i], 1, matched.first - 1)
        line.middle <- paste(val.matched.text)
        line.end <- substr(text[i], matched.last + 1, nchar(text[i]))
        text[i] <- paste(line.begin, line.middle, line.end, sep = "")
        
        reg <- regexpr(code.pattern, text[i], perl = TRUE)
      }
    }
  }

  return(text)
}
