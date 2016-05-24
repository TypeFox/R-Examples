sql.wrap <- function(x, term = "spec", operator = "=", BOOL = "OR", regex = TRUE, by = 500) {
  
  operator <- match.arg(operator, c("=", "~"))
  if ( !is.numeric(x) ) {
    if ( regex ){
      x <- gsub("'", ".", x) # e.g. "Gigantochloa_sp._'daluoensis'"
      x <- gsub("([(]|[)]|[+])", "[\\1]", x) # e.g. "Clivia_sp._RHA+CA_7b", "Oreobolus_sp._1_(Laegaard_70382)
    }
    if ( operator == "=" ) x <- paste("'", x, "'", sep = "")
    else x <- paste("'^", x, "$'", sep = "")
  }
  if ( !is.null(term) ) x <- paste(term, operator, x, sep = "")
  if ( !is.null(BOOL) ){
    id <- seq(from = 1, to = length(x), by = 500)
    id <- paste(id, c(id[-1] - 1, length(x)), sep = ":")
    id <- lapply(id, function(obj) eval(parse(text = obj)))
    x <- sapply(id, function(obj, id) paste(obj[id], collapse = paste("", BOOL, "")), obj = x)
  }
  x
}