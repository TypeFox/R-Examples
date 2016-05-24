# Read XML while fishing-out some errors

# Function `tryCatch` looks at conditions thrown by `read_xml`. If the condition
# is of class "Rcpp::exception" a handler function `handle_rcpp_exception` is
# executed. This function looks at the message of the original condition. If it
# looks like "PCDATA invalid Char value" it extracts the value of the invalid
# character and the error number (in square brackets in the original message). A
# new condition object is created of class `"invalid_char_value"` which extends
# the original class `"Rcpp::exception"`. Extracted offending character value
# and error number are added as attributes to the new condition. Finally, the
# condition is signalled.
#
# Conditions not caught by `tryCatch`, so all others apart from
# `"Rcpp::exception"`s, are signaled as usual.
#
# This allows for special handling of different conditions somewhere else

read_xml_with_errors <- function(x, ...) {
  # Handle Rcpp exceptions,
  # which include libxml errors
  handle_rcpp_exception <- function(cond) {
    # better error message
    cond$message <- paste("xml2::read_xml says:", cond$message)
    # enhance for selected errors
    if(grepl("PCDATA invalid Char value", cond$message)) {
      cond <- condition(
        subclass=c("invalid_char_value", class(cond)),
        message = cond$message,
        call=cond$call,
        error_no = as.numeric(stringr::str_extract(cond$message, "(?<=\\[)[0-9]+(?=\\]$)" )),
        char_value = as.numeric(stringr::str_extract(cond$message, "(?<=value )[0-9]+(?= \\[)" ))
      )
    }
    stop(cond)
  }

  # Catch!
  tryCatch( xml2::read_xml(x, ...),
            "Rcpp::exception" = handle_rcpp_exception
  )
}




# Read XML safely
#
# Removes illegal characters
# @param invalid_as NA or character. If not NA, replace invalid characters with `invalid_as`
read_xml_safely <- function(x, ..., xml_invalid_as=getOption("oai.xml_invalid_as", "") ) {
  repeat {
    tryCatch( return(read_xml_with_errors(x, ...)),
              # Removing offending characters
              invalid_char_value = function(er, invalid_as=xml_invalid_as) {
                charint <- attr(er, "char_value")
                if(is.na(invalid_as)) {
                  stop(er)
                } else {
                  stopifnot(is.character(xml_invalid_as))
                  x <<- gsub(intToUtf8(charint), invalid_as, x)
                  msg <- paste0(er$message, ", replacing offending characters with ", dQuote(invalid_as))
                  warning(msg)
                }
              }
    )
  }
}
