normVarNames <- function(vars, sep="_")
{
  if (sep == ".") sep <- "\\."
    
  # Replace all _ and . and ' ' with the nominated separator. Note
  # that I used [] originally but that fails so use |.
  
  pat  <- '_|\u00a0|\u2022| |,|-|:|/|&|\\.|\\?|\\[|\\]|\\{|\\}|\\(|\\)'
  rep  <- sep
  vars <- gsub(pat, rep, vars)

  # Replace any all capitals words with Initial capitals. This uses an
  # extended perl regular expression. The ?<! is a zero-width negative
  # look-behind assertion that matches any occurrence of the following
  # pattern that does not follow a Unicode property (the \p) of a
  # letter (L) limited to uppercase (u). Not quite sure of the
  # use-case for the look-behind.
  
  pat  <- '(?<!\\p{Lu})(\\p{Lu})(\\p{Lu}*)'
  rep  <- '\\1\\L\\2'
  vars <- gsub(pat, rep, vars, perl=TRUE)
  
  # Replace any capitals not at the beginning of the string with _ 
  # and then the lowercase letter.
  
  pat  <- '(?<!^)(\\p{Lu})'
  rep  <- paste0(sep, '\\L\\1')
  vars <- gsub(pat, rep, vars, perl=TRUE)
  
  # WHY DO THIS? Replace any number sequences not preceded by an
  # underscore, with it preceded by an underscore. The (?<!...) is a
  # lookbehind operator.
  
  pat  <- paste0('(?<![', sep, '\\p{N}])(\\p{N}+)')
  rep  <- paste0(sep, '\\1')
  vars <- gsub(pat, rep, vars, perl=TRUE)

  # Remove any resulting initial or trailing underscore or multiples:
  #
  # _2level -> 2level

  vars <- gsub("^_+", "", vars)
  vars <- gsub("_+$", "", vars)
  vars <- gsub("__+", "_", vars)
  
  # Convert to lowercase
  
  vars <- tolower(vars)
    
  # Remove repeated separators.
  
  pat  <- paste0(sep, "+")
  rep  <- sep
  vars <- gsub(pat, rep, vars)
  
  return(vars)
}
