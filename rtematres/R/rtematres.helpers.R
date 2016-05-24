# clean html strings
clean_html_string <- function(string) {
  remove_tags = gsub("<.*?>", "", string)
  remove_newlines = gsub("\\n"," ",remove_tags)
  replace_escaped_quotes = gsub("\"", "'", remove_newlines)
  remove_square_brackets = gsub("\\[.?\\]", "", replace_escaped_quotes)
  remove_empty_brackets = gsub("\\(\\s*\\)", "", remove_square_brackets)
  remove_trailing_ws = gsub("\\s*$", "", remove_empty_brackets)
  remove_leading_ws = gsub("^\\s*", "", remove_trailing_ws)
  remove_unnecessary_ws = gsub("[ ]{2,}", " ", remove_leading_ws)
  cleaned_string = remove_unnecessary_ws
  return(cleaned_string)
}

# a wrapper function of xpathSApply. which will return NA instead of zero-lenght list when nodes not found
xmlNodesValue <- function(doc, path){
  out = xpathSApply(doc, path, xmlValue, trim=T, ignoreComments=T)
  out = Filter(function(x) x!="", out)
  if (length(out) == 0) return(NA)
  out
}
