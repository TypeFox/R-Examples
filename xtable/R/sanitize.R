sanitize <- function(str, type = "latex") {
  if(type == "latex"){
    result <- str
    result <- gsub("\\\\", "SANITIZE.BACKSLASH", result)
    result <- gsub("$", "\\$", result, fixed = TRUE)
    result <- gsub(">", "$>$", result, fixed = TRUE)
    result <- gsub("<", "$<$", result, fixed = TRUE)
    result <- gsub("|", "$|$", result, fixed = TRUE)
    result <- gsub("{", "\\{", result, fixed = TRUE)
    result <- gsub("}", "\\}", result, fixed = TRUE)
    result <- gsub("%", "\\%", result, fixed = TRUE)
    result <- gsub("&", "\\&", result, fixed = TRUE)
    result <- gsub("_", "\\_", result, fixed = TRUE)
    result <- gsub("#", "\\#", result, fixed = TRUE)
    result <- gsub("^", "\\verb|^|", result, fixed = TRUE)
    result <- gsub("~", "\\~{}", result, fixed = TRUE)
    result <- gsub("SANITIZE.BACKSLASH", "$\\backslash$", result, fixed = TRUE)
    return(result)
  } else {
    result <- str
    result <- gsub("&", "&amp;", result, fixed = TRUE)
    result <- gsub(">", "&gt;", result, fixed = TRUE)
    result <- gsub("<", "&lt;", result, fixed = TRUE)
    return(result)
  }
}


sanitize.numbers <- function(str, type,
                             math.style.negative = FALSE,
                             math.style.exponents = FALSE){
  if (type == "latex"){
    result <- str
    if ( math.style.negative ) {
      for(i in 1:length(str)) {
        result[i] <- gsub("-", "$-$", result[i], fixed = TRUE)
      }
    }
    if ( math.style.exponents ) {
      if (is.logical(math.style.exponents) && ! math.style.exponents ) {
      } else if (is.logical(math.style.exponents) && math.style.exponents ||
                 math.style.exponents == "$$"
                 ) {
        for(i in 1:length(str)) {
          result[i] <-
            gsub("^\\$?(-?)\\$?([0-9.]+)[eE]\\$?(-?)\\+?\\$?0*(\\d+)$",
                 "$\\1\\2 \\\\times 10^{\\3\\4}$", result[i])
        }
      } else if (math.style.exponents == "ensuremath") {
        for(i in 1:length(str)) {
          result[i] <-
            gsub("^\\$?(-?)\\$?([0-9.]+)[eE]\\$?(-?)\\+?\\$?0*(\\d+)$",
                 "\\\\ensuremath{\\1\\2 \\\\times 10^{\\3\\4}}",
                 result[i])
        }
      } else if (math.style.exponents == "UTF8" ||
                 math.style.exponents == "UTF-8") {
        for(i in 1:length(str)) {
          ## this code turns 1e5 into a UTF-8 representation of 1\times10^5
          if (all(grepl("^\\$?(-?)\\$?([0-9.]+)[eE]\\$?(-?)\\+?\\$?0*(\\d+)$",
                        result[i]))) {
            temp <- strsplit(result[i],"eE",result[i])
            result[i] <-
              paste0(temp[1],
                     "\u00d710",
                     chartr("-1234567890",
                            "\u207b\u00b9\u00b2\u00b3\u2074\u2075\u20746\u20747\u20748\u20749\u2070",
                            temp[2]))
          }
        }
      }
    }
    return(result)
  } else {
    return(str)
  }
}


sanitize.final <- function(str, type){
  if (type == "latex"){
    return(str)
  } else {
    str$text <- gsub("  *", " ",  str$text, fixed = TRUE)
    str$text <- gsub(' align="left"',  "", str$text,
                     fixed = TRUE)
    return(str)
  }
}

### Some trivial helper functions
### Suggested by Stefan Edwards, sme@iysik.com
### Helper function for disabling sanitizing
as.is <- function(str) {str}

### Helper function for embedding names in a math environment
as.math <- function(str, ...) { paste0('$',str,'$', ...) }
