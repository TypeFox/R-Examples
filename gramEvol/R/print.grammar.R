print.grammar <- function (x, ..., max.line.len = 60) {
  
  # generic function to print n spaces
  gen.spaces <- function(n) do.call(paste0, as.list(rep(" ", n)))
  
  # max len of all rulenames
  maxLen = max(sapply (x$def, function(rule) nchar(rule[[1]])))
  
  for (rule in x$def) {
    spaces = maxLen - nchar(rule[[1]])  # use spaces as tabs usually don't work
    expressions = lapply(rule[[2]], function(x) unescape.gt.lt(as.character(x))) # extract character versions of rules
    
    # if current string is too long, add to the next expression, a '\r' and enough spaces
    for (i in 1:length(expressions)) {
      if (i + 1 > length(expressions)) break
      
      if (i == 1) {
        leadspace = maxLen + 4
      } else {
        leadspace = 0
      }
      
      if (leadspace + nchar(expressions[[i]]) > max.line.len) {
        expressions[[i + 1]] = paste0('\n', gen.spaces(maxLen + 7), expressions[[i + 1]])
      }
    }
    
    # print them with appropriate style
    cat(paste0("<", rule[[1]], ">" ), 
        gen.spaces(spaces),
        " ::= ",
        do.call(function(...) {paste(..., sep=" | ")}, expressions), 
        "\n", ..., sep = "")
  }
}

