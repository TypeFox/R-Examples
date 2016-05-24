#----------------------------------------------------------------------------
# localsolver
# Copyright (c) 2014, WLOG Solutions
#----------------------------------------------------------------------------


remove.comments <- function(modelText) {
  modelText <- gsub('//[^\n]*\n', '\n', modelText)  # mono-line comment
  modelText <- gsub('//[^\n]*$', '', modelText)   # last mono-line comment
  modelText <- gsub('/\\*.*?\\*/', '', modelText) # multi-line comment
  return(modelText)  
}

detect.functions <- function(modelText) {
  modelText <- gsub('\r\n', '\n', modelText) # fix eol
  modelText <- remove.comments(modelText)
  result <- list()

  functions <- unlist(strsplit(modelText, "\\s*function\\s+"))
  if (functions[[1]] == "") {
    functions <- functions[-1]
  }
  for(fun in functions) {
#    m <- regexec("^([_[:alpha:]][_[:alnum:]]*)[[:space:]]*(\\([^{]*\\))[[:space:]]*\\{(.*)}[[:space:]]*$", fun)
    m <- regexec("^([_a-zA-Z][_a-zA-Z0-9]*)[[:space:]]*(\\([^{]*\\))[[:space:]]*\\{(.*)}[[:space:]]*$", fun)
    
    rm <- regmatches(fun, m)[[1]]
    if (length(rm) == 0) {
      msg <- sprintf("Function declaration expected near '%s...'", substr(fun, 1, 25))
      stop(msg)
    }
    
    name <- rm[[2]]
    decl <- sprintf('function %s%s', name, rm[[3]])
    result[[name]] <- list(decl = decl, body = rm[[4]])
  }
  
  return(result)
}
