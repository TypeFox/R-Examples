# Look for OAI-PMH exceptions
# https://www.openarchives.org/OAI/openarchivesprotocol.html#ErrorConditions
# xml = parsed xml
# Return TRUE if OK or stop
# handle_errors <- function(xml) {
#   nodeset <- xml2::xml_find_all(xml, ".//*[local-name()='error']")
#   if( length(nodeset) > 0 ) {
#     msgs <- sapply(nodeset, function(n)
#       paste0( xml2::xml_attr(n, "code"), ": ", xml2::xml_text(n) ) )
#     stop( paste0("OAI-PMH exceptions: ", paste(msgs, collapse="\n")) )
#   }
#
# }
#




handle_errors <- function(xml) {
  # find error tags
  req <- xml2::xml_text(xml2::xml_find_one(xml, ".//*[local-name()='request']" ))
  nodeset <- xml2::xml_find_all(xml, ".//*[local-name()='error']")
  # collect error information, if any
  if( length(nodeset) == 0 ) {
    return(TRUE)
  } else {
    errors <- lapply(nodeset, function(n)
      c(code=xml2::xml_attr(n, "code"),
        message=xml2::xml_text(n) ) )
    cond <- condition(c("oai-pmh_error", "error"),
                      message = paste0("OAI-PMH errors: ",
                                       paste( sapply(errors, paste, collapse=": "), collapse=",\n")),
                      request=req,
                      error_codes = sapply(errors, "[", "code") )
    stop(cond)
  }
}