evalText <- function(text, where = .GlobalEnv) {
   for(expr in parse(text = text))
     value <- eval(expr, envir = where)
   value
}
