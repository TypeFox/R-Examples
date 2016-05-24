#----------------------------------------------------------------------------
# localsolver
# Copyright (c) 2014, WLOG Solutions
#----------------------------------------------------------------------------

parse.ls.error <- function(lsp, func.locations, err.file) {
  error.str <- readLines(err.file, warn=FALSE)[[1]]

  # syntax error parsing
  m <- regmatches(error.str, regexec("^Error:  At line ([[:digit:]]+): (.+)$", error.str))[[1]]
  if (length(m) != 0) {
    line.no <- as.integer(m[[2]])
    err.text <- m[[3]]
    return(sprintf("%s\n%s", err.text, find.ls.context(lsp, func.locations, line.no)))
  }

  # runtime error parsing
  m <- regmatches(error.str, regexec("^^Error \\[[^]]+ line: ([[:digit:]]+) ]: (.+)$", error.str))[[1]]
  if (length(m) != 0) {
    line.no <- as.integer(m[[2]])
    err.text <- m[[3]]
    return(sprintf("%s\n%s", err.text, find.ls.context(lsp, func.locations, line.no)))
  }
  
  return(error.str)
}

find.ls.context <- function(lsp, func.locations, line.no) {
  for(funcName in names(func.locations)) {
    loc <- func.locations[[funcName]]
    if (loc$start.line <= line.no && line.no < loc$end.line) {
      inFuncLine <- line.no - loc$start.line + 1
      funcDesc <- lsp$functions[[funcName]]

      lines <- paste(".", strsplit(funcDesc$body, '\n')[[1]])
      lines[[inFuncLine]] <- sprintf('>%s', substring(lines[[inFuncLine]], 2))
      
      ctx <- lines[max(1, inFuncLine-1):min(inFuncLine+1, length(lines))]
      if (inFuncLine <= 2) {
        ctx[[1]] <- sprintf('%s {%s', funcDesc$decl, ctx[[1]])        
      }
      
      return(sprintf("in '%s' function, arround\n%s", funcName, paste(ctx, collapse='\n')))
    }
  }
  return("At unexpected position (line: %d)", line.no)
}
